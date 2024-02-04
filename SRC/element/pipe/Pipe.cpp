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

// Minjie
#include <CrdTransf.h>
#include <Node.h>
#include <Pipe.h>

#include <cmath>

void *OPS_PipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "Invalid #args,  want: element pipe "
                  "tag? nd1? nd2? transfTag? pipeMatTag? pipeSecTag?"
                  "<-T0 T0? -p p? -cMass? -releasey? -releasez?>\n";
        return 0;
    }

    // get tag
    int iData[6];
    int numData = 6;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid integer input for pipe element\n";
        return 0;
    }

    // get data
    double T0 = 0.0, pressure = 0.0;
    int cMass = 0;
    int releasez = 0;
    int releasey = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *theType = OPS_GetString();
        if (strcmp(theType, "-T0") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &T0) < 0) {
                    opserr << "WARNING: failed to read T0\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-p") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &pressure) < 0) {
                    opserr << "WARNING: failed to read internal "
                              "pressure\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-cMass") == 0) {
            cMass = 1;
        } else if (strcmp(theType, "-releasez") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &releasez) < 0) {
                    opserr << "WARNING: failed to get releasez";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-releasey") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &releasey) < 0) {
                    opserr << "WARNING: failed to get releasey";
                    return 0;
                }
            }
        }
    }

    auto *theSect = dynamic_cast<PipeSection *>(
        OPS_getSectionForceDeformation(iData[5]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[5]
               << " is not found or not a pipe section\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[4]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[4]
               << " is not found or not a pipe material\n";
        return 0;
    }

    auto *theTrans = OPS_getCrdTransf(iData[3]);
    if (theTrans == 0) {
        opserr << "WARNING: CrdTransf " << iData[3]
               << " is not found\n";
        return 0;
    }

    auto *ele =
        new Pipe(iData[0], iData[1], iData[2], *theTrans, *theMat,
                 *theSect, T0, pressure, cMass, releasez, releasey);

    return ele;
}

Pipe::Pipe()
    : ElasticBeam3d(),
      theMat(0),
      theSect(0),
      T0(0.0),
      pressure(0.0) {}

Pipe::Pipe(int tag, int nd1, int nd2, CrdTransf &theTransf,
           PipeMaterial &mat, PipeSection &sect, double t0,
           double pre, int cm, int rz, int ry)
    : ElasticBeam3d(ELE_TAG_Pipe),
      theMat(0),
      theSect(0),
      T0(t0),
      pressure(pre) {
    // nodes
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    // section
    theSect = dynamic_cast<PipeSection *>(sect.getCopy());
    if (theSect == 0) {
        opserr << "Pipe element - failed to "
                  "get a copy of section "
               << sect.getTag() << "\n";
        exit(-1);
    }

    // material
    theMat = dynamic_cast<PipeMaterial *>(mat.getCopy());
    if (theMat == 0) {
        opserr << "Pipe element - failed to get a copy of "
                  "material with tag "
               << mat.getTag() << "\n";
        exit(-1);
    }

    // transf
    theCoordTransf = theTransf.getCopy3d();

    if (!theCoordTransf) {
        opserr << "Pipe element -- failed to get "
                  "copy of coordinate transformation\n";
        exit(-1);
    }

    // Make no release if input not 0, 1, 2, or 3
    releasez = rz;
    releasey = ry;
    if (releasez < 0 || releasez > 3) {
        releasez = 0;
    }
    if (releasey < 0 || releasey > 3) {
        releasey = 0;
    }

    cMass = cm;
}

Pipe::~Pipe() {
    if (theMat != 0) {
        delete theMat;
    }
    if (theSect != 0) {
        delete theSect;
    }
}

const char *Pipe::getClassType(void) const { return "Pipe"; };

void Pipe::setDomain(Domain *theDomain) {
    // check domain
    if (theDomain == 0) {
        opserr << "Pipe::setDomain -- Domain is null\n";
        exit(-1);
    }

    int ndm = OPS_GetNDM();
    if (ndm != 3) {
        opserr << "WARNING: pipe element must be 3D\n";
        exit(-1);
    }

    // get nodes
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

    if (theNodes[0] == 0) {
        opserr << "Pipe::setDomain  tag: " << this->getTag()
               << " -- Node 1: " << connectedExternalNodes(0)
               << " does not exist\n";
        exit(-1);
    }

    if (theNodes[1] == 0) {
        opserr << "Pipe::setDomain  tag: " << this->getTag()
               << " -- Node 2: " << connectedExternalNodes(1)
               << " does not exist\n";
        exit(-1);
    }

    // check dofs
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != getNumDOF() / 2) {
        opserr << "Pipe::setDomain  tag: " << this->getTag()
               << " -- Node 1: " << connectedExternalNodes(0)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    if (dofNd2 != getNumDOF() / 2) {
        opserr << "Pipe::setDomain  tag: " << this->getTag()
               << " -- Node 2: " << connectedExternalNodes(1)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    // update section data
    if (updateSectionData() < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return;
    }

    // update material data
    if (updateMaterialData() < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return;
    }

    // set domain
    this->ElasticBeam3d::setDomain(theDomain);
}

double Pipe::aveTemp() {
    // get average element temperature
    double Ti = theNodes[0]->getTemp();
    double Tj = theNodes[1]->getTemp();
    double Tavg = 0.5 * (Ti + Tj);
    return Tavg - T0;
}

int Pipe::updateMaterialData() {
    // select point based on temperature
    int retVal = 0;
    auto Tpt = theMat->selectPoint(aveTemp(), retVal);
    if (retVal < 0) {
        return retVal;
    }

    // set ElasticBeam3d data
    ElasticBeam3d::E = Tpt.E;
    double nu = Tpt.xnu;
    ElasticBeam3d::G = ElasticBeam3d::E / (2 * (1.0 + nu));

    return retVal;
}

int Pipe::updateSectionData() {
    ElasticBeam3d::A = theSect->AREA();
    ElasticBeam3d::Jx = theSect->JX();
    ElasticBeam3d::Iy = theSect->IY();
    ElasticBeam3d::Iz = theSect->IZ();
    ElasticBeam3d::rho = theSect->RHO();
    return 0;
}

void Pipe::zeroLoad(void) {
    this->ElasticBeam3d::zeroLoad();

    double alp = theSect->ALFAV();
    double temp = aveTemp();
    if (temp > 0) {
        ElasticBeam3d::q0[0] -= E * A * alp * temp;
    }
}