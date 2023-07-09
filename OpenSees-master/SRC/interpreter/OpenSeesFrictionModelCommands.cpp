/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Description: command to create friction models

#include <FrictionModel.h>
#include <elementAPI.h>
#include <map>

void* OPS_Coulomb();
void* OPS_VelDependent();
void* OPS_VelPressureDep();
void* OPS_VelDepMultiLinear();
void* OPS_VelNormalFrcDep();

namespace {

    struct char_cmp {
        bool operator () (const char *a, const char *b) const
        {
            return strcmp(a, b) < 0;
        }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap frictionModelsMap;


    static int setUpFrictionModels(void)
    {
        frictionModelsMap.insert(std::make_pair("Coulomb", &OPS_Coulomb));
        frictionModelsMap.insert(std::make_pair("VelDependent", &OPS_VelDependent));
        frictionModelsMap.insert(std::make_pair("VelPressureDep", &OPS_VelPressureDep));
        frictionModelsMap.insert(std::make_pair("VelDepMultiLinear", &OPS_VelDepMultiLinear));
        frictionModelsMap.insert(std::make_pair("VelNormalFrcDep", &OPS_VelNormalFrcDep));

        return 0;
    }
}

int
OPS_FrictionModel()
{
    static bool initDone = false;
    if (initDone == false) {
        setUpFrictionModels();
        initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING too few arguments: frictionModel type? tag? ...\n";
        return -1;
    }

    const char* frnType = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = frictionModelsMap.find(frnType);
    if (iter == frictionModelsMap.end()) {
        opserr << "WARNING friction model type " << frnType << " is unknown\n";
        return -1;
    }

    FrictionModel *theFrnMdl = (FrictionModel*)(*iter->second)();
    if (theFrnMdl == 0) {
        return -1;
    }

    // Now add the friction model to the modelBuilder
    if (OPS_addFrictionModel(theFrnMdl) == false) {
        opserr << "ERROR could not add friction model.\n";
        delete theFrnMdl;
        return -1;
    }

    return 0;
}
