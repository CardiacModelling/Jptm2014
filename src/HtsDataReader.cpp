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

#include <map>

#include "HtsDataReader.hpp"
#include "Exception.hpp"


void HtsDataReader::SetUpMap()
{
    mChannelNameToIdMap.clear();
    mChannelNameToIdMap["hERG"]  = 0u;
    mChannelNameToIdMap["NaV1.5"]= 1u;
    mChannelNameToIdMap["IKs"]   = 2u;
    mChannelNameToIdMap["Ito,f"] = 3u;
    mChannelNameToIdMap["CaV1.2"]= 4u;
}

void HtsDataReader::LoadALine(std::stringstream& rLine)
{
    std::string name;
    c_vector<double, 5> ic50s;
    c_vector<double, 5> hills;
    ic50s.clear();
    hills.clear();

    rLine >> name;
    // herg
    rLine >> ic50s(0);
    rLine >> hills(0);
    // nav
    rLine >> ic50s(1);
    rLine >> hills(1);
    // iks
    rLine >> ic50s(2);
    rLine >> hills(2);
    // ito
    rLine >> ic50s(3);
    rLine >> hills(3);
    // cav
    rLine >> ic50s(4);
    rLine >> hills(4);

    this->mDrugNames.push_back(name);
    this->mIc50values.push_back(ic50s);
    this->mHillCoefficients.push_back(hills);
}

unsigned HtsDataReader::GetChannelIndex(const std::string& rChannelName)
{
    std::map<std::string,unsigned>::iterator it = mChannelNameToIdMap.find(rChannelName);

    if (it==mChannelNameToIdMap.end())
    {
        EXCEPTION("No channel '" << rChannelName << "' in data file.");
    }

    return it->second;
}

unsigned HtsDataReader::GetNumDrugs(void)
{
    return mDrugNames.size();
}


std::string HtsDataReader::GetDrugName(unsigned drugIndex)
{
    assert(drugIndex < GetNumDrugs());
    return mDrugNames[drugIndex];
}

unsigned HtsDataReader::GetDrugIndex(const std::string& rName)
{
    unsigned idx = UINT_MAX;
    for (unsigned i=0; i<mDrugNames.size(); ++i)
    {
        if (mDrugNames[i] == rName)
        {
            idx = i;
            break;
        }
    }
    if (idx==UINT_MAX)
    {
        EXCEPTION("Drug " << rName << " not found.");
    }
    return idx;
}

double HtsDataReader::GetIC50Value(const std::string& rDrugName, const std::string& rChannelName)
{
    unsigned channel_index = GetChannelIndex(rChannelName);
    unsigned drug_index = GetDrugIndex(rDrugName);

    assert(drug_index < GetNumDrugs());
    assert(channel_index < 5u);

    return mIc50values[drug_index](channel_index);
}

double HtsDataReader::GetHillCoefficient(const std::string& rDrugName, const std::string& rChannelName)
{
    unsigned channel_index = GetChannelIndex(rChannelName);
    unsigned drug_index = GetDrugIndex(rDrugName);

    assert(drug_index < GetNumDrugs());
    assert(channel_index < 5u);

    double hill = mHillCoefficients[drug_index](channel_index);
    if (hill < 0)
    {   // Should default to 1
        hill = 1.0;
    }
    return hill;
}
