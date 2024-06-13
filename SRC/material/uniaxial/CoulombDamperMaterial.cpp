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
               "tag? Tangent? FrictionForce? -tol tol? -numFlipped "
               "numFlipped? -reduceFc? -dampOutTangent? "
               "dampOutTangent\n";
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

    double tol = 1e-6;
    double dampOutTangent = 1.0;
    int method = 1;
    int numFlipped = 2;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-tol") == 0 &&
            OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numData, &tol) < 0) {
                opserr << "WARNING: failed to get tol\n";
                return 0;
            }
        } else if (strcmp(opt, "-numFlipped") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numData, &numFlipped) < 0) {
                opserr << "WARNING: failed to get numFlipped\n";
                return 0;
            }
        } else if (strcmp(opt, "-dampOutTangent") == 0 &&
                   OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetDoubleInput(&numData, &dampOutTangent) < 0) {
                opserr << "WARNING: failed to get dampOutTangent\n";
                return 0;
            }
            if (dampOutTangent > 0) {
                method = 3;
            }
        } else if (strcmp(opt, "-reduceFc") == 0) {
            method = 2;
        }
    }

    // Parsing was successful, allocate the material
    theMaterial =
        new CoulombDamperMaterial(iData[0], dData[0], dData[1], tol,
                                  dampOutTangent, method, numFlipped);
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type "
                  "CoulombDamperMaterial"
               << endln;
        return 0;
    }

    return theMaterial;
}

CoulombDamperMaterial::CoulombDamperMaterial(int tag, double k,
                                             double fc, double t,
                                             double damp, int m,
                                             int n)
    : UniaxialMaterial(tag, MAT_TAG_CoulombDamperMaterial),
      trialStrain(0.0),
      trialStrainRate(0.0),
      tangent(k),
      friction(fc),
      commitTrialStrainRate(0),
      flipped(0),
      tol(t),
      dampOutTangent(damp),
      method(m),
      numFlipped(n),
      parameterID(0) {}

CoulombDamperMaterial::CoulombDamperMaterial()
    : UniaxialMaterial(0, MAT_TAG_CoulombDamperMaterial),
      trialStrain(0.0),
      trialStrainRate(0.0),
      tangent(0.0),
      friction(0.0),
      commitTrialStrainRate(0),
      flipped(0),
      tol(1e-6),
      dampOutTangent(1.0),
      method(1),
      numFlipped(2),
      parameterID(0) {}

CoulombDamperMaterial::~CoulombDamperMaterial() {
    // does nothing
}

int CoulombDamperMaterial::setTrialStrain(double strain,
                                          double strainRate) {
    trialStrain = strain;
    trialStrainRate = strainRate;

    // flipped
    if ((commitTrialStrainRate > 0 && trialStrainRate < -tol) ||
        (commitTrialStrainRate < 0 && trialStrainRate > tol)) {
        ++flipped;
    }

    return 0;
}

int CoulombDamperMaterial::setTrial(double strain, double &stress,
                                    double &tan, double strainRate) {
    trialStrain = strain;
    trialStrainRate = strainRate;

    stress = tangent * strain + sign();
    tangent = tan;

    // flipped
    if ((commitTrialStrainRate > 0 && trialStrainRate < -tol) ||
        (commitTrialStrainRate < 0 && trialStrainRate > tol)) {
        ++flipped;
    }

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
    commitTrialStrainRate = trialStrainRate;
    flipped = 0;
    return 0;
}

int CoulombDamperMaterial::revertToLastCommit(void) {
    trialStrainRate = commitTrialStrainRate;
    flipped = 0;
    return 0;
}

int CoulombDamperMaterial::revertToStart(void) {
    trialStrain = 0.0;
    trialStrainRate = 0.0;
    commitTrialStrainRate = 0.0;
    flipped = 0;
    return 0;
}

UniaxialMaterial *CoulombDamperMaterial::getCopy(void) {
    CoulombDamperMaterial *theCopy = new CoulombDamperMaterial(
        this->getTag(), tangent, friction, tol, dampOutTangent,
        method, numFlipped);
    theCopy->trialStrain = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    theCopy->commitTrialStrainRate = commitTrialStrainRate;
    theCopy->parameterID = parameterID;
    return theCopy;
}

int CoulombDamperMaterial::sendSelf(int cTag, Channel &theChannel) {
    int res = 0;
    static Vector data(9);
    data(0) = this->getTag();
    data(1) = tangent;
    data(2) = friction;
    data(3) = parameterID;
    data(4) = commitTrialStrainRate;
    data(5) = tol;
    data(6) = dampOutTangent;
    data(7) = method;
    data(8) = numFlipped;
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
    static Vector data(9);
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
        commitTrialStrainRate = data(4);
        tol = data(5);
        dampOutTangent = data(6);
        method = (int)data(7);
        numFlipped = (int)data(8);
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
    double dampTangent = friction / tol;

    if (flipped > numFlipped) {
        // rate flipped n times
        if (method == 1) {
            res = factor() * dampTangent * trialStrainRate;
        } else if (method == 2) {
            if (trialStrainRate > tol) {
                res = factor() * friction;
            } else if (trialStrainRate < -tol) {
                res = -factor() * friction;
            }
        } else if (dampOutTangent == 3) {
            res = dampOutTangent * trialStrainRate;
        }

    } else if (trialStrainRate > -tol && trialStrainRate < tol) {
        // linear rate
        res = dampTangent * trialStrainRate;

    } else {
        // rate not flipped
        if (trialStrainRate > tol) {
            // rate > 0
            res = friction;
        } else if (trialStrainRate < -tol) {
            // rate < 0
            res = -friction;
        }
    }

    return res;
}

double CoulombDamperMaterial::dsign() {
    double res = 0.0;
    double dampTangent = friction / tol;

    if (flipped > numFlipped) {
        // rate flipped n times
        if (method == 1) {
            res = factor() * dampTangent;
        } else if (method == 2) {
            res = 0.0;
        } else if (method == 3) {
            res = dampOutTangent;
        }
    } else if (trialStrainRate > -tol && trialStrainRate < tol) {
        // linear rate
        res = dampTangent;

    } else {
        res = 0.0;
    }

    return res;
}

double CoulombDamperMaterial::factor() {
    double res = 1.0;
    for (int i = 0; i < flipped; i += 2) {
        res *= 0.5;
    }
    return res;
}