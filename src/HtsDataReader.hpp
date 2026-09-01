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

#ifndef HTSDATAREADER_HPP_
#define HTSDATAREADER_HPP_

#include <fstream>
#include <vector>

#include "UblasVectorInclude.hpp"
#include "FileFinder.hpp"
#include "AbstractDrugDataStructure.hpp"

/**
 * A little class which is just designed to read in data
 * from the input_data folder and look it up nicely.
 */
class HtsDataReader: public AbstractDrugDataStructure<5>
{
private:
    /** Set up the #mChannelNameToIdMap */
    void SetUpMap();

protected:

    /**
     * A method which hard-codes the format of this file
     * @param rLine  a line in stringsteam format.
     */
    void LoadALine(std::stringstream& rLine);

    /** A map between the name of a channel and its index in the vectors this class returns */
    std::map<std::string, unsigned> mChannelNameToIdMap;

    /**
     * Default Constructor (empty)
     */
    HtsDataReader(){};

    /**
     * Helper method to get the index of the channel as stored in a vector
     * @param rChannelName  The name of the channel
     */
    unsigned GetChannelIndex(const std::string& rChannelName);

    /**
     * Read a header line if present.
     *
     * We just say we have read it and then move on.
     */
    bool LoadHeaderLine(std::stringstream& rLine)
    {
        return true;
    }

public:

    /**
     * Constructor. Data will be immediately loaded into memory by the constructor.
     *
     * @param fileName  the name of a file to load (relative or absolute).
     */
    HtsDataReader(std::string fileName)
      : AbstractDrugDataStructure<5>()
    {
        LoadDataFromFile(fileName);
        SetUpMap();
    };

    /**
     * Constructor. Data will be immediately loaded into memory by the constructor.
     *
     * @param rFileFinder  a file finder pointing to the data file to load.
     */
    HtsDataReader(FileFinder& rFileFinder)
      : AbstractDrugDataStructure<5>()
    {
        LoadDataFromFile(rFileFinder.GetAbsolutePath());
        SetUpMap();
    };

    /**
     * Destructor (empty)
     */
    virtual ~HtsDataReader(){};

    /**
     * @return The number of drugs in the data file.
     */
    unsigned GetNumDrugs(void);

    /**
     * Get the name of a drug from its index in the data file.
     *
     * @param drugIndex  The index of the drug (row of the data file on which it appears)
     * @return  The name of the drug
     */
    std::string GetDrugName(unsigned drugIndex);

    /**
     * Get the index of a drug (row of data file) from its name.
     *
     * @param rName  The name of the drug
     * @return  The index in the current drug list.
     */
    unsigned GetDrugIndex(const std::string& rName);

    /**
     * Return the IC50 value associated with a particular channel.
     *
     * @param drugIndex  index of the drug as listed in rows of the data file.
     * @param channelName  The channel name, NaV1.5, hERG, CaV1.2, IKs or Ito,f
     * @return the IC50 value for a certain drug on a certain channel.
     */
    double GetIC50Value(const std::string& rDrugName, const std::string& rChannelName);

    /**
     * Return the hill coefficient associated with this drug and this channel's dose-reponse curve
     *
     * @param drugIndex  The index of the drug (in drug_data.dat file)
     * @param channelName  The channel name, NaV1.5, hERG, CaV1.2, IKs or Ito,f
     * @return the hill coefficient
     */
    double GetHillCoefficient(const std::string& rDrugName, const std::string& rChannelName);

};
#endif // HTSDATAREADER_HPP_
