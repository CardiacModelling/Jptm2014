/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTMAKELOOKUPTABLE_HPP_
#define TESTMAKELOOKUPTABLE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>

#include "CheckpointArchiveTypes.hpp"

#include "SetupModel.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "LookupTableGenerator.hpp"
#include "LookupTableReader.hpp"

/**
 * This class precomputes APD for a hyper-cube in parameter space, with the corners
 * [0,1] in each of 'DIM' dimensions.
 *
 * The resulting data table is then used, with interpolation, to estimate
 * APD for any combination of the DIM parameters by the main simulation code.
 */
class TestMakeLookupTable : public CxxTest::TestSuite
{
private:
    std::string mFileName;
    std::string mArchiveFilename;

public:

    // Change this to decide how many parameters (ion channel blocks) to sweep over.
    static const unsigned DIM = 5u;

    void TestGenerateAndArchiveLookupTable() throw (Exception)
    {
        /*
         * Parameters to change here.
         */
        const unsigned model_index = 2u; // Ten tusscher '06
        const double hertz = 1.0;
        const unsigned max_num_evaluations = 1000000u;

        // You shouldn't need to change any settings under here.

        SetupModel setup_model(hertz, model_index);
        const std::string model_name = setup_model.GetModel()->GetSystemName();

        std::stringstream file_name;
        OutputFileHandler handler("LookupTables", false); // Don't wipe the output folder!
        if (DIM==5)
        {
            file_name << "5d_hERG_IKs_INa_ICaL_Ito";
        }
        else if (DIM==4)
        {
            file_name << "4d_hERG_IKs_INa_ICaL";
        }
        else if(DIM==2)
        {
            file_name << "2d_hERG_IKs";
        }
        else
        {
            EXCEPTION("Not set up for DIMs other than 2, 4 or 5.");
        }
        file_name << "_" << hertz << "Hz";
        mFileName = file_name.str();

        mArchiveFilename = handler.GetOutputDirectoryFullPath() + model_name + "/" + mFileName + "_generator.arch";

        FileFinder archive_file(mArchiveFilename, RelativeTo::Absolute);

        LookupTableGenerator<DIM>* p_generator;
        if (archive_file.IsFile())
        {
            // Create a pointer to the input archive
            std::ifstream ifs(mArchiveFilename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_generator;
        }
        else
        {
            // Set up a new generator
            p_generator = new LookupTableGenerator<DIM>(model_index, mFileName, "LookupTables/" + model_name);
            p_generator->SetPacingFrequency(hertz);

            p_generator->SetParameterToScale("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0 , 1.0);
            p_generator->SetParameterToScale("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0 , 1.0);
            if (DIM >= 4u)
            {
                p_generator->SetParameterToScale("membrane_fast_sodium_current_conductance", 0.0 , 1.0);
                p_generator->SetParameterToScale("membrane_L_type_calcium_current_conductance", 0.0 , 1.0);
            }
            if (DIM >= 5u)
            {
                // At the moment we are tending to use TT06 for speed and robustness, this should really 
                // be fast Ito in models that have it separately.
                p_generator->SetParameterToScale("membrane_transient_outward_current_conductance", 0.0, 1.0);
            }
            p_generator->AddQuantityOfInterest(Apd90, 1 /*ms*/); // QoI and tolerance
            p_generator->AddQuantityOfInterest(Apd50, 1 /*ms*/); // QoI and tolerance
            p_generator->SetMaxNumPaces(30u*60u); // This is 30 minutes of 1Hz pacing
            p_generator->SetMaxVariationInRefinement(5u); // This prevents over-refining in one area.
        }

        const unsigned start_evaluations = p_generator->GetNumEvaluations();
        std::cout << "Started with " << start_evaluations << " evaluations.\n";
        const unsigned evaluations_per_checkpoint = 500;

        for (unsigned i=start_evaluations;
             i < max_num_evaluations;
             i += evaluations_per_checkpoint)
        {
            // Run some evaluations
            p_generator->SetMaxNumEvaluations(i+evaluations_per_checkpoint);
            p_generator->GenerateLookupTable();

            // Overwrite archive entry
            {
                LookupTableGenerator<DIM>* const p_arch_generator = p_generator;

                std::ofstream ofs(mArchiveFilename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                output_arch << p_arch_generator;
            }
        }
        delete p_generator;
    }

    void TestInterpolationViaGenerator() throw(Exception)
    {
        LookupTableGenerator<DIM>* p_generator;

        // Create a pointer to the input archive
        std::ifstream ifs(mArchiveFilename.c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        // restore from the archive
        input_arch >> p_generator;

        // Do the interpolation directly with the generator and the
        // parameter box methods.
        // Some sampling points we would like to get an estimate at.
        std::vector<c_vector<double,DIM> > sampling_points;

        c_vector<double,DIM> sample;
        for (double u = 0; u<=1.0; u+=0.01)
        {
            sample[0] = u;
            for (double v = 0; v<=1.0; v+=0.01)
            {
                sample[1] = v;
                if (DIM==4)
                {
                    for (double w = 0; w<=1.0; w+=0.1)
                    {
                        sample[2] = w;
                        for (double x = 0; x<=1.0; x+=0.1)
                        {
                            sample[3] = x;
                            sampling_points.push_back(sample);
                        }
                    }
                }
                else if (DIM==2)
                {
                    sampling_points.push_back(sample);
                }
                else
                {
                    EXCEPTION("Not set up to deal with DIM other than 2 or 4.");
                }
            }
        }

        std::cout << "\nInterpolated QoIs\n";
        std::vector<std::vector<double> > predictions = p_generator->Interpolate(sampling_points);

        assert(predictions.size()==sampling_points.size());

        for (unsigned j=0; j<predictions.size(); j++) // Loop over parameter space points
        {
            std::cout << sampling_points[j][0] << "\t" << sampling_points[j][1] ;
            if (DIM==4)
            {
                std::cout << "\t" << sampling_points[j][2] << "\t" << sampling_points[j][3];
            }
            for (unsigned i=0; i< predictions[j].size(); i++) // Loop over QoIs.
            {
                std::cout << "\t" << predictions[j][i];
            }
            std::cout << std::endl;
        }

        delete p_generator;
    }
};

#endif // TESTMAKELOOKUPTABLE_HPP_
