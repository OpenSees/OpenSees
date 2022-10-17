/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

#include <Channel.h>
#include <CoulombDamperMaterial.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <Parameter.h>
#include <Vector.h>
#include <elementAPI.h>
#include <string.h>

#include <cmath>

void *OPS_CoulombDamperMaterial(void) {
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr
            << "Invalid #args,  want: uniaxialMaterial CoulombDamper "
               "tag? Tangent? FrictionForce? flipDamping? factor? \n";
        return 0;
    }

    int iData[1];
    double dData[2];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid tag for uniaxialMaterial "
                  "CoulombDamper\n";
        return 0;
    }

    numData = 2;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for uniaxial CoulombDamper "
               << iData[0] << "\n";
        return 0;
    }

    double flipDamping = 1.0, factor = -1.0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &flipDamping) < 0) {
            opserr << "WARNING: invalid tol for CoulombDamper\n";
            return 0;
        }
    }
    if (flipDamping <= 0) {
        flipDamping = 1.0;
    }
    if (OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &factor) < 0) {
            opserr << "WARNING: invalid tol for CoulombDamper\n";
            return 0;
        }
    }
    if (factor >= 1.0) {
        factor = -1.0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new CoulombDamperMaterial(
        iData[0], dData[0], dData[1], flipDamping, factor);
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type "
                  "CoulombDamperMaterial"
               << endln;
        return 0;
    }

    return theMaterial;
}

CoulombDamperMaterial::CoulombDamperMaterial(int tag, double k,
                                             double fc, double d,
                                             double f)
    : UniaxialMaterial(tag, MAT_TAG_CoulombDamperMaterial),
      trialStrain(0.0),
      trialStrainRate(0.0),
      tangent(k),
      friction(fc),
      prevTrialStrainRate(0),
      flipDamping(d),
      factor(f),
      curFactor(1.0),
      parameterID(0) {}

CoulombDamperMaterial::CoulombDamperMaterial()
    : UniaxialMaterial(0, MAT_TAG_CoulombDamperMaterial),
      trialStrain(0.0),
      trialStrainRate(0.0),
      tangent(0.0),
      friction(0.0),
      prevTrialStrainRate(0),
      flipDamping(1.0),
      factor(-1.0),
      curFactor(1.0),
      parameterID(0) {}

CoulombDamperMaterial::~CoulombDamperMaterial() {
    // does nothing
}

int CoulombDamperMaterial::setTrialStrain(double strain,
                                          double strainRate) {
    if (factor > 0) {
        // prev is the last iteration
        // when choosing to reduce friction
        prevTrialStrainRate = trialStrainRate;
    }
    trialStrain = strain;
    trialStrainRate = strainRate;

    return 0;
}

int CoulombDamperMaterial::setTrial(double strain, double &stress,
                                    double &tan, double strainRate) {
    if (factor > 0) {
        // prev is the last iteration
        // when choosing to reduce friction
        prevTrialStrainRate = trialStrainRate;
    }
    trialStrain = strain;
    trialStrainRate = strainRate;

    stress = tangent * strain + sign();
    tangent = tan;

    return 0;
}

double CoulombDamperMaterial::getStress(void) {
    return tangent * trialStrain + sign();
}

double CoulombDamperMaterial::getTangent(void) { return tangent; }

double CoulombDamperMaterial::getDampTangent(void) { return dsign(); }

double CoulombDamperMaterial::getInitialTangent(void) {
    return tangent;
}

int CoulombDamperMaterial::commitState(void) {
    prevTrialStrainRate = trialStrainRate;
    curFactor = 1.0;
    return 0;
}

int CoulombDamperMaterial::revertToLastCommit(void) {
    trialStrainRate = prevTrialStrainRate;
    curFactor = 1.0;
    return 0;
}

int CoulombDamperMaterial::revertToStart(void) {
    trialStrain = 0.0;
    trialStrainRate = 0.0;
    prevTrialStrainRate = 0.0;
    curFactor = 1.0;
    return 0;
}

UniaxialMaterial *CoulombDamperMaterial::getCopy(void) {
    CoulombDamperMaterial *theCopy = new CoulombDamperMaterial(
        this->getTag(), tangent, friction, flipDamping, factor);
    theCopy->trialStrain = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    theCopy->prevTrialStrainRate = prevTrialStrainRate;
    theCopy->parameterID = parameterID;
    return theCopy;
}

int CoulombDamperMaterial::sendSelf(int cTag, Channel &theChannel) {
    int res = 0;
    static Vector data(7);
    data(0) = this->getTag();
    data(1) = tangent;
    data(2) = friction;
    data(3) = parameterID;
    data(4) = prevTrialStrainRate;
    data(5) = flipDamping;
    data(6) = factor;
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "CoulombDamperMaterial::sendSelf() - failed to "
                  "send data"
               << endln;

    return res;
}

int CoulombDamperMaterial::recvSelf(int cTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker) {
    int res = 0;
    static Vector data(7);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);

    if (res < 0) {
        opserr << "CoulombDamperMaterial::recvSelf() - failed to "
                  "receive data"
               << endln;
        tangent = 0;
        friction = 0;
        this->setTag(0);
        friction = 0;
    } else {
        this->setTag(int(data(0)));
        tangent = data(1);
        friction = data(2);
        parameterID = (int)data(3);
        prevTrialStrainRate = data(4);
        flipDamping = data(5);
        factor = data(6);
    }

    return res;
}

void CoulombDamperMaterial::Print(OPS_Stream &s, int flag) {
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "CoulombDamperMaterial tag: " << this->getTag() << endln;
        s << "  Tangent: " << tangent
          << ", Friciton force: " << friction << "\n";
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"CoulombDamperMaterial\", ";
        s << "\"Tangent\": " << tangent << ", ";
        s << "\"Friction\": " << friction << "}";
    }
}

int CoulombDamperMaterial::setParameter(const char **argv, int argc,
                                        Parameter &param) {
    if (strcmp(argv[0], "Tangent") == 0) {
        param.setValue(tangent);
        return param.addObject(1, this);
    }
    if (strcmp(argv[0], "Friction") == 0) {
        param.setValue(friction);
        return param.addObject(2, this);
    }
    return -1;
}

int CoulombDamperMaterial::updateParameter(int parameterID,
                                           Information &info) {
    switch (parameterID) {
        case 1:
            tangent = info.theDouble;
            return 0;
        case 2:
            friction = info.theDouble;
            return 0;
        default:
            return -1;
    }
}

int CoulombDamperMaterial::activateParameter(int paramID) {
    parameterID = paramID;

    return 0;
}

double CoulombDamperMaterial::getStressSensitivity(int gradIndex,
                                                   bool conditional) {
    if (parameterID == 1) return trialStrain;
    if (parameterID == 2) return sign();

    return 0.0;
}

double CoulombDamperMaterial::getTangentSensitivity(int gradIndex) {
    if (parameterID == 1) return 1.0;
    return 0.0;
}

double CoulombDamperMaterial::getDampTangentSensitivity(
    int gradIndex) {
    if (parameterID == 2) {
        return dsign();
    }
    return 0.0;
}

double CoulombDamperMaterial::getInitialTangentSensitivity(
    int gradIndex) {
    if (parameterID == 1) return 1.0;
    return 0.0;
}

int CoulombDamperMaterial::commitSensitivity(double strainGradient,
                                             int gradIndex,
                                             int numGrads) {
    // Nothing to commit ... path independent
    return 0;
}

double CoulombDamperMaterial::sign() {
    double res = 0.0;

    if (factor > 0) {
        // this is to reduce friction to reach to converge
        double tol = 1e-5;
        if (trialStrainRate > tol) {
            res = factor * friction;
        } else if (trialStrainRate < -tol) {
            res = -factor * friction;
        } else {
            res = friction / tol * trialStrainRate;
        }

        if ((prevTrialStrainRate > 0 && trialStrainRate < 0) ||
            (prevTrialStrainRate < 0 && trialStrainRate > 0)) {
            curFactor *= factor;
        }
    } else {
        // this is to use damping for flipped velocity
        if (prevTrialStrainRate > 0) {
            // prev is > 0, only positive friction
            if (trialStrainRate > 0) {
                // constant positive friction
                res = friction;
            } else {
                // damping for negative velocity
                res = trialStrainRate * flipDamping + friction;
            }
        }

        if (prevTrialStrainRate < 0) {
            // prev is < 0, only negative friction
            if (trialStrainRate < 0) {
                // constant negative friction
                res = -friction;
            } else {
                // damping for positive velocity
                res = trialStrainRate * flipDamping - friction;
            }
        }
    }

    return res;
}

double CoulombDamperMaterial::dsign() {
    double res = 0.0;

    if (factor > 0) {
        // this is to reduce friction to reach to converge
        double tol = 1e-5;
        if (trialStrainRate > tol) {
            res = 0.0;
        } else if (trialStrainRate < -tol) {
            res = 0.0;
        } else {
            res = friction / tol;
        }
    } else {
        // this is to use damping for flipped velocity
        if (prevTrialStrainRate > 0) {
            // prev is > 0, only positive friction
            if (trialStrainRate > 0) {
                // constant positive friction
                res = 0.0;
            } else {
                // damping for negative velocity
                res = flipDamping;
            }
        }

        if (prevTrialStrainRate < 0) {
            // prev is < 0, only negative friction
            if (trialStrainRate < 0) {
                // constant negative friction
                res = 0.0;
            } else {
                // damping for positive velocity
                res = flipDamping;
            }
        }
    }

    return res;
}