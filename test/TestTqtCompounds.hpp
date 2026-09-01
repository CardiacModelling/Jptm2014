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

#ifndef TESTTQTCOMPOUNDS_HPP_
#define TESTTQTCOMPOUNDS_HPP_

#include <boost/shared_ptr.hpp>
#include <cxxtest/TestSuite.h>

// Chaste source includes
#include "CheckpointArchiveTypes.hpp"

#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "SteadyStateRunner.hpp"
#include "OutputFileHandler.hpp"
#include "CellProperties.hpp"

// ApPredict includes
#include "SetupModel.hpp"
#include "AbstractDataStructure.hpp"
#include "LogisticDistribution.hpp"
#include "LogLogisticDistribution.hpp"
#include "BayesianInferer.hpp"
#include "SingleActionPotentialPrediction.hpp"
#include "LookupTableGenerator.hpp"

// This project includes
#include "HtsDataReader.hpp"

// Should always be last, can do drugs in parallel with this...
#include "PetscSetupAndFinalize.hpp"

/**
 * This test will run the Bayesian inference scheme as per the Elkins paper,
 * and then use the lookup table generator to get credible intervals instantly
 * rather than running hundreds of simulations as we did before.
 */
class TestTqtCompounds : public CxxTest::TestSuite
{
private:
    /** Frequency of pacing to use */
    double mFreq;

    std::vector<double> GetIc50Samples(const double& ic50,
                                       const double& rSigma,
                                       const unsigned& rNumSamples)
    {
        // We will say that the underlying distribution really is
        // a logistic with this pIc50 and sigma.
        LogisticDistribution logistic;
        BayesianInferer inferer(PIC50);

        std::vector<double> observed_pic50s;
        observed_pic50s.push_back(AbstractDataStructure::ConvertIc50ToPic50(ic50));

        inferer.SetObservedData(observed_pic50s);
        inferer.SetSpreadOfUnderlyingDistribution(rSigma);
        inferer.PerformInference();

        // We now take samples from the inferred PDF
        std::vector<double> samples = inferer.GetSampleMedianValue(rNumSamples);

        // Convert them from pIC50 to IC50 for dose-response calculations.
        double temp_mean = 0u;
        for (unsigned i=0; i<samples.size(); i++)
        {
            temp_mean += samples[i];
            samples[i] = AbstractDataStructure::ConvertPic50ToIc50(samples[i]);
        }
        temp_mean /= samples.size();
        std::cout << "Original iC50 = " << ic50 << ", mean of inferred samples = " << AbstractDataStructure::ConvertPic50ToIc50(temp_mean) << std::endl << std::flush;

        return samples;
    }

    LookupTableGenerator<5u>* LoadLookupTable(const std::string& rModelName)
    {
        LookupTableGenerator<5u>* p_generator;

        OutputFileHandler handler("LookupTables", false); // Don't wipe the output folder!

        std::stringstream filename;
        filename << handler.GetOutputDirectoryFullPath() << rModelName << "/5d_hERG_IKs_INa_ICaL_Ito_" << mFreq << "Hz_generator.arch";

        FileFinder file_finder(filename.str(), RelativeTo::Absolute);
        if (!file_finder.IsFile())
        {
            EXCEPTION("Archive file not found, looking for: " << file_finder.GetAbsolutePath());
        }

        // Create a pointer to the input archive
        std::ifstream ifs((filename.str()).c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        // restore from the archive
        input_arch >> p_generator;

        return p_generator;
    }

public:
    /**
     * This test will wipe $CHASTE_TEST_OUTPUT/TQTStudy/<Model>/<options>
     *
     * Parameters are defined as command line arguments.
     *
     */
    void TestDrugAffectByModifyingConductances(void) throw (Exception)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments.
        std::cout << "# " << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc == 1)
        {
            std::cerr << "TestTqtCompounds::Please input arguments\n"
                         "* EITHER --herg, OR --herg-cav-nav, OR --herg-cav-nav-iks-ito\n"
                         "         Whether to add experimental noise to herg, or 3 or 5 channels.\n"
                         "* --number-runs\n"
                         "    start with this number of simulations (must be >0). Apart from the\n"
                         "    first one, the rest are with randomly inferred parameters.\n"
                         "\n"
                         "* --model\n"
                         "    The index of the model to use (see SetupModel.hpp for options).\n"
                         "* [--drug-file]\n"
                         "    a drug filename (default is input_data/AZ_HTS_data.txt).\n"
                         "* [--drugs]\n"
                         "    the name of particular drugs to re-run the analysis on\n"
                         "    (these must be in the correct capitalisation as per the above file).\n"
                         "* [--set-hills-to-one]\n"
                         "    set all the hill coefficients equal to one (defaults to false).\n"
                         "* [--frequency]\n"
                         "    Frequency of pacing to use in Hz (defaults to 1Hz).\n"
                         "\n"
                         "In all the above '[argument]' denotes optional.\n\n";
            return;
        }

        mFreq = 1.0; // Default is 1Hz.
        if (CommandLineArguments::Instance()->OptionExists("--frequency"))
        {
            mFreq = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("--frequency");
        }

        // Pick the model, set the correct pacing frequency etc...
        // This class reads the command line for the model index with argument "--model"
        SetupModel setup_model(mFreq);
        boost::shared_ptr<AbstractCvodeCell> p_model = setup_model.GetModel();
        p_model->UseCellMLDefaultStimulus();

        bool set_hills_to_one = false;
        if (CommandLineArguments::Instance()->OptionExists("--set-hills-to-one"))
        {
            set_hills_to_one = true;
        }

        unsigned channels_to_vary; // This will take 0 for hERG only,
                                    // 1 for hERG, CaV1.2, NaV1.5
                                    // 2 for hERG, Cav1.2, NaV1.5, IKS, Ito,f

        if (CommandLineArguments::Instance()->OptionExists("--herg"))
        {
            channels_to_vary = 0u;
            if (   (CommandLineArguments::Instance()->OptionExists("--herg-cav-nav"))
                || (CommandLineArguments::Instance()->OptionExists("--herg-cav-nav-iks-ito")) )
            {
                EXCEPTION("You can have EITHER --herg or --herg-cav-nav or --herg-cav-nav-iks-ito, not more than one");
            }
        }
        else if (CommandLineArguments::Instance()->OptionExists("--herg-cav-nav"))
        {
            channels_to_vary = 1u;
        }
        else if (CommandLineArguments::Instance()->OptionExists("--herg-cav-nav-iks-ito"))
        {
            channels_to_vary = 2u;
            if  (CommandLineArguments::Instance()->OptionExists("--herg-cav-nav"))
            {
                EXCEPTION("You can have EITHER --herg or --herg-cav-nav or --herg-cav-nav-iks-ito, not more than one");
            }
        }
        else
        {
            EXCEPTION("Please specify EITHER --herg or --herg-cav-nav or --herg-cav-nav-iks-ito");
        }

        std::string drug_file = "projects/Jptm2014/input_data/AZ_HTS_data.txt";
        if (CommandLineArguments::Instance()->OptionExists("--drug-file"))
        {
            drug_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--drug-file");
        }

        std::vector<std::string> drugs_to_run;
        bool run_for_all_drugs = true;
        if (CommandLineArguments::Instance()->OptionExists("--drugs"))
        {
            drugs_to_run = CommandLineArguments::Instance()->GetStringsCorrespondingToOption("--drugs");
            run_for_all_drugs = false;
        }

        unsigned number_runs = 0u;
        if (!CommandLineArguments::Instance()->OptionExists("--number-runs"))
        {
            EXCEPTION("Please also specify number of runs to perform with --number-runs");
        }
        else
        {
            number_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--number-runs");
            if (number_runs==0u)
            {
                EXCEPTION("The number of runs must be positive.");
            }
        }

        // To prevent each process doing the same random exploration of lookup tables.
        RandomNumberGenerator::Instance()->Reseed(PetscTools::GetMyRank());
        std::cout << "Beginning " << number_runs-1 << " random runs for each drug\n" << std::flush;

        ///////// END DEFINE PARAMETERS ////////////////////////

        // Load drug data
        FileFinder file(drug_file, RelativeTo::ChasteSourceRoot);

        if(!file.IsFile())
        {
            EXCEPTION("Drug data file: " <<  file.GetAbsolutePath() << " does not exist.");
        }

        //
        // Based on the data file, specify appropriate levels of control variability.
        //
        double sigma_na_pic50;
        double sigma_cal_pic50;
        double sigma_herg_pic50;
        double sigma_iks_pic50;
        double sigma_ito_pic50;
        if (   file.GetLeafNameNoExtension()=="AZ_HTS_data"
            || file.GetLeafNameNoExtension()=="AZ_HTS_data_Gintant_hERG")
        {
            // The following are hardcoded from a matlab script which analyses how much variation
            // IonWorks pIC50s have based upon repeated controls in KATE or AZ's equivalent database...
            // See the Elkins 2013 Jptm paper for more details on this.
            sigma_na_pic50 =  0.0760348; // From AZ IonWorks Quattro Nav1.5 control- n=2307
            sigma_cal_pic50 =  0.159676; // From AZ IonWorks Quattro Cav control
            sigma_herg_pic50 = 0.103456; // From AZ data: very similar to within 4 dp of GSK data
            sigma_iks_pic50 = 0.139736283; // Data averaged from all GSK (n=79, 451) and AZ data(n=525)- Average weighted on n's
            sigma_ito_pic50 = 0.0859756; // From AZ IonWorks Quattro Ito,f control
            //        const double beta_na_hill =  11.97; // From AZ IonWorks Quattro Nav1.5 control- n=2307
            //        const double beta_cal_hill =  8.295; // From AZ IonWorks Quattro Cav control
            //        const double beta_herg_hill = 5.605; // From AZ IonWorks hERG data
            //        const double beta_iks_hill = 6.0325; // Average of AZ and GSK IKs data (similar)
            //        const double beta_ito_hill = 11.6; // From AZ IonWorks Quattro Ito,f control
        }
        else if (file.GetLeafNameNoExtension()=="GSK_HTS_data")
        {
            // The following are hardcoded from a matlab script which analyses how much variation
            // IonWorks pIC50s have based upon repeated controls in KATE or AZ's equivalent database...
            // it is in RyanElkins/data/graphs/fit_IC50_histograms.m
            sigma_na_pic50 =  0.1388; // From GSK IonWorks Quattro Nav1.5 control
            sigma_cal_pic50 =  0.12; // From GSK Barracuda data
            sigma_herg_pic50 = 0.146; // From GSK Barracuda data
            sigma_iks_pic50 = 0.139736283; // Data averaged from all GSK (n=79, 451) and AZ data(n=525)- Average weighted on n's
            sigma_ito_pic50 = 0.0859756; // From AZ IonWorks Quattro Ito,f control (Not actually used as all data is pIC50 = 0).
        }
        else
        {
            if (number_runs > 1u)
            {
                EXCEPTION("Undefined control variability for this data file.");
            }
        }

        // Load the drug data from file.
        HtsDataReader drug_data(file);

        // The following names are fixed and correspond to metadata tags.
        const double default_g_na = p_model->GetParameter("membrane_fast_sodium_current_conductance");
        const double default_g_cal= p_model->GetParameter("membrane_L_type_calcium_current_conductance");
        const double default_g_kr = p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
        const double default_g_ks = p_model->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");

        double g_to;
        if (p_model->HasParameter("membrane_fast_transient_outward_current_conductance"))
        {
            g_to = p_model->GetParameter("membrane_fast_transient_outward_current_conductance");
        }
        else if (p_model->HasParameter("membrane_transient_outward_current_conductance"))
        {
            WARNING(p_model->GetSystemName() << " model does not have fast Ito tagged, so general Ito current (fast + slow) was drug blocked.");
            g_to = p_model->GetParameter("membrane_transient_outward_current_conductance");
        }
        else
        {
            EXCEPTION(p_model->GetSystemName() << " model has neither fast nor 'general' Ito currents, fast preferred for tagging.");
        }
        const double default_g_to = g_to;

        boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = p_model->GetSystemInformation();
        std::string model_name = p_model->GetSystemName();

        // Set up foldernames for each model and protocol set.
        std::string foldername = "TQTStudy/" + model_name + "/" + file.GetLeafNameNoExtension();
        std::stringstream steady_state_foldername;

        steady_state_foldername << foldername << "/ic50_varying_";
        if (set_hills_to_one)
        {
            steady_state_foldername << "hill_fixed_to_one_";
        }

        steady_state_foldername << 2*channels_to_vary+1 << "_channel"; // Hardcoded to 1,3,5 channels!

        // Make the above directories
        // don't clean them
        OutputFileHandler steady_handler(steady_state_foldername.str(),false);

        // Speed things up by recording the control steady state vars and using where appropriate.
        bool control_recorded = false;
        N_Vector solution_at_control;

        /**
         * DECIDE WHAT CONCENTRATIONS TO TEST AT, in uM (micro M)
         */
        double divisions = 4.0; // per log unit.
        std::vector<double> drug_conc;
        {
            // This range of drug concentrations corresponds to those shown in Figure 1 of
            // Redfern-et-al (2003) doi:10.1016/S0008-6363(02)00846-5
            drug_conc.push_back(0);
            unsigned i=0;
            while (drug_conc.back()<100)
            {
                drug_conc.push_back(pow(10.0,((double)(i)-3.0*divisions)/divisions));
                i++;
            }
        }

        LookupTableGenerator<5u>* p_lookup_table;
        unsigned max_runs = UNSIGNED_UNSET;
        if (number_runs > 1u)
        {
            p_lookup_table = LoadLookupTable(model_name);
            max_runs = p_lookup_table->GetMaxNumPaces();
        }

        // This object actually runs the steady state action potential predictions.
        SingleActionPotentialPrediction ap_predictor(p_model);
        ap_predictor.SetLackOfOneToOneCorrespondenceIsError();

        // NB O'Hara-Rudy seemed to give some fairly different results with and without 1800 pace cap.
        // Particularly for strong hERG block where O'Hara fails to repolarise.
        // Ten Tusscher '06 didn't have this problem.
        ap_predictor.SetMaxNumPaces(max_runs);

        /**
         * START LOOP OVER EACH DRUG
         *
         */
        std::cout << "Running simulations..." << std::endl;
        for(unsigned drug_index = 0; drug_index < drug_data.GetNumDrugs(); drug_index++)
        {
            if (PetscTools::IsParallel() && drug_index % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
            {
                // Let another processor do this drug
                continue;
            }
            // Loop over each drug
            const std::string this_drug_name = drug_data.GetDrugName(drug_index);
            std::cout << "DRUG = " << this_drug_name;
            if (run_for_all_drugs==false)
            {
                // if this drug (in loop over file) is not in the list of ones we want to re-run, then...
                if (find(drugs_to_run.begin(), drugs_to_run.end(), this_drug_name)==drugs_to_run.end())
                {
                    std::cout << " - skipping" << std::endl << std::flush;
                    // ...skip to the next one
                    continue;
                }
            }
            std::cout << std::endl << std::flush;

            // If we ever end up running this for multiple experimental repeats
            // the reader will need changing to return a vector of IC50s,
            // which will then get passed to the inference scheme.
            const double mean_na_ic50 = drug_data.GetIC50Value(this_drug_name,"NaV1.5");
            const double mean_cal_ic50 = drug_data.GetIC50Value(this_drug_name,"CaV1.2");
            const double mean_herg_ic50 = drug_data.GetIC50Value(this_drug_name,"hERG");
            const double mean_iks_ic50 = drug_data.GetIC50Value(this_drug_name,"IKs");
            const double mean_ito_ic50 = drug_data.GetIC50Value(this_drug_name,"Ito,f");

            const double mean_na_hill = drug_data.GetHillCoefficient(this_drug_name,"NaV1.5");
            const double mean_cal_hill = drug_data.GetHillCoefficient(this_drug_name,"CaV1.2");
            const double mean_herg_hill = drug_data.GetHillCoefficient(this_drug_name,"hERG");
            const double mean_iks_hill = drug_data.GetHillCoefficient(this_drug_name,"IKs");
            const double mean_ito_hill = drug_data.GetHillCoefficient(this_drug_name,"Ito,f");

            // Open files and write headers
            out_stream steady_voltage_results_file = steady_handler.OpenOutputFile("voltage_results_" + this_drug_name + ".dat");
            //out_stream steady_calcium_results_file = steady_handler.OpenOutputFile("calcium_results_" + drug_data.GetDrugName(drug_index) + ".dat");
            *steady_voltage_results_file << "RandomRun\tNaV1.5_IC50(uM)\tNaV1.5_Hill\tCaV1.2_IC50(uM)\tCaV1.2_Hill\thERG_IC50(uM)\thERG_Hill\tIKs_IC50(uM)\tIKs_Hill\tItof_IC50(uM)\tItof_Hill\tConcentration(uM)\tAPD50(ms)\tAPD90(ms)\n";
            //*steady_voltage_results_file << "RandomRun\tNaV1.5_IC50(uM)\tNaV1.5_Hill\tCaV1.2_IC50(uM)\tCaV1.2_Hill\thERG_IC50(uM)\thERG_Hill\tIKs_IC50(uM)\tIKs_Hill\tItof_IC50(uM)\tItof_Hill\tConcentration(uM)\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\n";
            //*steady_calcium_results_file << "RandomRun\tNaV1.5_IC50(uM)\tNaV1.5_Hill\tCaV1.2_IC50(uM)\tCaV1.2_Hill\thERG_IC50(uM)\thERG_Hill\tIKs_IC50(uM)\tIKs_Hill\tItof_IC50(uM)\tItof_Hill\tConcentration(uM)\tPeakChangeInCa(mM/ms)\tPeakCa(mM)\tCaD50(ms)\tCaD90(ms)\n";

            // Get the Bayesian samples for IC50 values
            std::vector<double> herg_ic50_samples = GetIc50Samples(mean_herg_ic50, sigma_herg_pic50, number_runs-1);
            std::vector<double> na_ic50_samples = GetIc50Samples(mean_na_ic50, sigma_na_pic50,       number_runs-1);
            std::vector<double> cal_ic50_samples = GetIc50Samples(mean_cal_ic50, sigma_cal_pic50,    number_runs-1);
            std::vector<double> iks_ic50_samples = GetIc50Samples(mean_iks_ic50, sigma_iks_pic50,    number_runs-1);
            std::vector<double> ito_ic50_samples = GetIc50Samples(mean_ito_ic50, sigma_ito_pic50,    number_runs-1);

            // Choose IC50 values
            double na_ic50 = mean_na_ic50;
            double cal_ic50 = mean_cal_ic50;
            double herg_ic50 = mean_herg_ic50;
            double iks_ic50 = mean_iks_ic50;
            double ito_ic50 = mean_ito_ic50;

            // Choose Hill coefficients
            double na_hill = mean_na_hill;
            double cal_hill = mean_cal_hill;
            double herg_hill = mean_herg_hill;
            double iks_hill = mean_iks_hill;
            double ito_hill = mean_ito_hill;

            // Perform Monte-Carlo style simulations...
            for (unsigned random_idx = 0u; random_idx < number_runs; random_idx++)
            {
                // First run uses the 'mean' or experimentally measured IC50 value.
                if (random_idx > 0)
                {
                    // Next runs use ones we have got from the Bayesian inference procedure.

                    // Pick random numbers from the distributions we associate with these IC50s...
                    herg_ic50 = herg_ic50_samples[random_idx-1];
                    if (channels_to_vary > 0u) // If not the vary herg only case
                    {
                        na_ic50 = na_ic50_samples[random_idx-1];
                        cal_ic50 = cal_ic50_samples[random_idx-1];
                    }
                    if (channels_to_vary >= 2) // If the vary everything case
                    {
                        iks_ic50 = iks_ic50_samples[random_idx-1];
                        ito_ic50 = ito_ic50_samples[random_idx-1];
                    }
                }

                // If we are setting Hill coefficients to one, then just overwrite whatever else has been calculated.
                if (set_hills_to_one)
                {
                    herg_hill = 1.0;
                    na_hill = 1.0;
                    cal_hill = 1.0;
                    iks_hill = 1.0;
                    ito_hill = 1.0;
                }

                std::cout << "Random run #" << random_idx << ":\n"
                        "NaV1.5 IC50 = " << na_ic50 << " uM,\tHill = " << na_hill << std::endl <<
                        "CaV1.2 IC50 = " << cal_ic50 << " uM,\tHill = " << cal_hill << std::endl <<
                        "hERG   IC50 = " << herg_ic50<< " uM,\tHill = " << herg_hill << std::endl <<
                        "Iks    IC50 = " << iks_ic50 << " uM,\tHill = " << iks_hill << std::endl <<
                        "Ito,f  IC50 = " << ito_ic50 << " uM,\tHill = " << ito_hill << std::endl << std::flush;

                /**
                 * START LOOP OVER EACH CONCENTRATION TO TEST WITH
                 */
                for (unsigned conc_index=0; conc_index<drug_conc.size(); conc_index++)
                {
                    // Because of the way the Lookup table is made, it is safest to compare with
                    // a max number of runs from the control case each time.
                    if (control_recorded)
                    {
                        p_model->SetStateVariables(solution_at_control);
                    }

                    std::cout << this_drug_name << ": " << drug_conc[conc_index] << " uM\n" << std::flush;

                    // Here we calculate the proportion of the different channels which are still active
                    // (at this concentration of this drug)
                    double gNa_factor = AbstractDataStructure::CalculateConductanceFactor(drug_conc[conc_index],na_ic50,na_hill);
                    double gCaL_factor= AbstractDataStructure::CalculateConductanceFactor(drug_conc[conc_index],cal_ic50,cal_hill);
                    double gKr_factor = AbstractDataStructure::CalculateConductanceFactor(drug_conc[conc_index],herg_ic50,herg_hill);
                    double gKs_factor = AbstractDataStructure::CalculateConductanceFactor(drug_conc[conc_index],iks_ic50,iks_hill);
                    double gTo_factor = AbstractDataStructure::CalculateConductanceFactor(drug_conc[conc_index],ito_ic50,ito_hill);

                    std::cout << "gNa factor = " << gNa_factor << "\n" << std::flush;
                    std::cout << "gCaL factor = "<< gCaL_factor<< "\n" << std::flush;
                    std::cout << "gKr factor = " << gKr_factor << "\n" << std::flush;
                    std::cout << "gKs factor = " << gKs_factor << "\n" << std::flush;
                    std::cout << "gTo factor = " << gTo_factor << "\n" << std::flush;

                    // The following names are fixed and correspond to metadata tags.
                    p_model->SetParameter("membrane_fast_sodium_current_conductance",
                                          default_g_na*gNa_factor);
                    p_model->SetParameter("membrane_L_type_calcium_current_conductance",
                                          default_g_cal*gCaL_factor);
                    p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance",
                                          default_g_kr*gKr_factor);
                    p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance",
                                          default_g_ks*gKs_factor);

                    if (p_model->HasParameter("membrane_fast_transient_outward_current_conductance"))
                    {
                        p_model->SetParameter("membrane_fast_transient_outward_current_conductance",
                                              default_g_to*gTo_factor);
                    }
                    else
                    {
                        p_model->SetParameter("membrane_transient_outward_current_conductance",
                                              default_g_to*gTo_factor);
                    }

                    /**
                     * STEADY STATE PACING EXPERIMENT
                     */
                    double apd90 = 0;
                    double apd50 = 0;

                    //
                    // We only run a real simulation for the first case (no 'noise' on experimental inputs).
                    //
                    if (random_idx==0)
                    {
                        try
                        {
                            ap_predictor.RunSteadyPacingExperiment(drug_conc[conc_index]);
                            apd90 = ap_predictor.GetApd90();
                            apd50 = ap_predictor.GetApd50();
                        }
                        catch (Exception& e)
                        {
                            // Here we handle evaluation errors in the same way as the Lookup table generator does.
                            /// \todo combine these bits of code inside SingleActionPotentialPrediction ??
                            std::string error_code = ap_predictor.GetErrorMessage();
                            std::cout << "Action potential evaluation error code = " << error_code << "\n" << std::flush;
                            if (  (error_code=="NoActionPotential_2" || error_code =="NoActionPotential_3") )
                            {
                                // For an APD calculation failure on repolarisation put in the stimulus period.
                                apd50 = 1000.0/mFreq;
                                apd90 = 1000.0/mFreq;
                            }
                            else
                            {
                                // For everything else (failure to depolarize "NoActionPotential_1") just put in zero for now.
                                apd50 = 0.0;
                                apd90 = 0.0;
                            }
                        }

                        //upstroke = ap_predictor.GetUpstrokeVelocity();
                        //peak = ap_predictor.GetPeakVoltage();

                        ap_predictor.Reset();

                        // The first ever time round at control record the state variables for using again later.
                        if (abs(drug_conc[conc_index])<1e-10 && !control_recorded)
                        {
                            control_recorded = true;
                            solution_at_control = p_model->GetStateVariables(); // This needs changing to a std::vector when using Chaste cells.
                        }
                    }

                    //
                    // Otherwise the action potential predictions come from Lookup table.
                    //
                    if (number_runs > 1u)
                    {
                        // In the lookup table the order of parameters is given in the filename:
                        // "5d_hERG_IKs_INa_ICaL_Ito_generator.arch"
                        c_vector<double,5u> sample_required_at;
                        sample_required_at[0] = gKr_factor;
                        sample_required_at[1] = gKs_factor;
                        sample_required_at[2] = gNa_factor;
                        sample_required_at[3] = gCaL_factor;
                        sample_required_at[4] = gTo_factor;

                        std::vector<c_vector<double,5u> > sampling_points;
                        sampling_points.push_back(sample_required_at);

                        std::vector<std::vector<double> > predictions = p_lookup_table->Interpolate(sampling_points);
                        assert(predictions.size()==1);

                        if (random_idx==0)
                        {
                            // Check that the interpolated values make sense
                            // (don't alter them)
                            std::cout << "Evaluated APD50 = " <<  apd50 << ", interpolated = " << predictions[0][1] << "\n";
                            std::cout << "Evaluated APD90 = " <<  apd90 << ", interpolated = " << predictions[0][0] << "\n";
                        }
                        else
                        {
                            // Take the interpolated values from the lookup table.
                            apd50 = predictions[0][1];
                            apd90 = predictions[0][0];
                        }
                    }

                    std::cout << "1Hz APD50 = " << apd50 << ", APD90 = " << apd90 << "\n" << std::flush;
                    //std::cout << "1Hz Upstroke velocity = " << upstroke << ", Peak mV = " << peak << ", APD50 = " << apd50 << ", APD90 = " << apd90 << "\n" << std::flush;
                    *steady_voltage_results_file << random_idx << "\t" << na_ic50 << "\t" << na_hill <<
                            "\t" << cal_ic50 <<  "\t" <<cal_hill <<
                            "\t" << herg_ic50  << "\t" <<herg_hill <<
                            "\t" << iks_ic50  << "\t" <<iks_hill <<
                            "\t" << ito_ic50  << "\t" <<ito_hill <<
                            "\t" << drug_conc[conc_index] << "\t" << apd50  << "\t" << apd90 << "\n";
                            //"\t" << drug_conc[conc_index] << "\t" << upstroke << "\t" << peak << "\t" << apd50  << "\t" << apd90 << "\n";

                    // Write out results...
                    *steady_voltage_results_file << std::flush;
                    //*steady_calcium_results_file << std::flush;
                }// Conc
            } // Random loop

            // Tidy up
            steady_voltage_results_file->close();
            //steady_calcium_results_file->close();
        } // Drug

        // Tidy up CVODE vector.
        DeleteVector(solution_at_control);

        if (number_runs > 1u) // Otherwise 'new' was never called...
        {
            delete p_lookup_table;
        }
    } // Test
};

#endif // TESTTQTCOMPOUNDS_HPP_
