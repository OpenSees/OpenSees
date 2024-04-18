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
#include <ElementalLoad.h>
#include <LinearCrdTransf3d.h>
#include <Node.h>
#include <Pipe.h>

#include <cmath>

void *OPS_PipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "Invalid #args,  want: element pipe "
                  "tag? nd1? nd2? pipeMatTag? pipeSecTag?"
                  "<-T0 T0? -p p? -cMass? -releasey releasey? "
                  "-releasez releasez?>\n";
        return 0;
    }

    // get tag
    int iData[5];
    int numData = 5;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid integer input for pipe element\n";
        return 0;
    }

    // get data
    double T0 = 0.0, pressure = 0.0;
    int cMass = 0;
    int releasez = 0;
    int releasey = 0;
    numData = 1;
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
        OPS_getSectionForceDeformation(iData[4]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[4]
               << " is not found or not a pipe section\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[3]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[3]
               << " is not found or not a pipe material\n";
        return 0;
    }

    auto *ele =
        new Pipe(iData[0], iData[1], iData[2], *theMat, *theSect, T0,
                 pressure, cMass, releasez, releasey);

    return ele;
}

Pipe::Pipe()
    : ElasticBeam3d(),
      theMat(0),
      theSect(0),
      alp(0.0),
      nu(0.0),
      T0(0.0),
      pressure(0.0) {}

Pipe::Pipe(int tag, int classTag)
    : ElasticBeam3d(tag, classTag),
      theMat(0),
      theSect(0),
      alp(0.0),
      nu(0.0),
      T0(0.0),
      pressure(0.0) {}

Pipe::Pipe(int tag, int nd1, int nd2, PipeMaterial &mat,
           PipeSection &sect, double t0, double pre, int cm, int rz,
           int ry)
    : ElasticBeam3d(tag, ELE_TAG_Pipe),
      theMat(0),
      theSect(0),
      alp(0.0),
      nu(0.0),
      T0(t0),
      pressure(pre) {
    if (createPipe(nd1, nd2, mat, sect, cm, rz, ry, pre) < 0) {
        opserr << "WARNING: failed to create pipe element\n";
        exit(-1);
    }
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

    const auto &crdsI = theNodes[0]->getCrds();
    const auto &crdsJ = theNodes[1]->getCrds();

    Vector IJ = crdsJ;
    IJ -= crdsI;
    IJ.Normalize();

    Vector gzAxis(ndm);
    gzAxis(2) = 1.0;

    Vector dir;
    if (crossProduct(IJ, gzAxis, dir) < 0) {
        exit(-1);
    }
    dir.Normalize();

    Vector zAxis(ndm);
    if (dir.Norm() < 0.1) {
        // parallel to global z axis
        zAxis(0) = 1.0;
    } else {
        zAxis(2) = 1.0;
    }

    if (ElasticBeam3d::theCoordTransf != 0) {
        delete ElasticBeam3d::theCoordTransf;
        ElasticBeam3d::theCoordTransf = 0;
    }
    ElasticBeam3d::theCoordTransf =
        new LinearCrdTransf3d(nextTransfTag(), zAxis);
    if (ElasticBeam3d::theCoordTransf == 0) {
        opserr << "WARNING: failed to crete Transformation object -- "
                  "CurvedPipe\n";
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
    nu = Tpt.xnu;
    ElasticBeam3d::G = ElasticBeam3d::E / (2 * (1.0 + nu));
    alp = Tpt.alp;
    if (E <= 0) {
        opserr << "E <= 0\n";
        return -1;
    }
    if (G <= 0) {
        opserr << "G <= 0\n";
        return -1;
    }
    if (alp <= 0) {
        opserr << "alp <= 0\n";
        return -1;
    }

    return retVal;
}

int Pipe::updateSectionData() {
    ElasticBeam3d::A = theSect->AREA();
    ElasticBeam3d::Jx = theSect->JX();
    ElasticBeam3d::Iy = theSect->IY();
    ElasticBeam3d::Iz = theSect->IZ();
    ElasticBeam3d::rho = theSect->RHO();
    double alphaV = theSect->ALFAV();
    if (alphaV > 99) {
        alphaV = 0.0;
    }
    ElasticBeam3d::alphaVz = alphaV;
    ElasticBeam3d::alphaVy = alphaV;
    return 0;
}

void Pipe::zeroLoad(void) {
    // update section data
    if (updateSectionData() < 0) {
        opserr << "Pipe::zeroLoad failed to update section data\n";
        return;
    }

    // update material data
    if (updateMaterialData() < 0) {
        opserr << "Pipe::zeroLoad failed to update material data\n";
        return;
    }

    this->ElasticBeam3d::zeroLoad();

    // due to thermal
    double temp = aveTemp();
    if (temp > 0) {
        ElasticBeam3d::q0[0] -= E * A * alp * temp;
    }

    // due to internal pressure
    if (pressure != 0) {
        double dout = theSect->DOUT();
        double thk = theSect->WALL();

        ElasticBeam3d::q0[0] -=
            0.25 * pressure * (dout - thk) * (1. - 2 * nu) * A / thk;
    }
}

int Pipe::addLoad(ElementalLoad *theLoad, double loadFactor) {
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = theCoordTransf->getInitialLength();
    double B1, B2, C1, C2;  // shear coefficients
    shearCoefficients(B1, B2, C1, C2);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
        // transformation matrix
        Matrix T(3, 3);
        Vector xAxis(3), yAxis(3), zAxis(3);
        Vector global_w(3);
        if (theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis) < 0) {
            return -1;
        }
        xAxis.Normalize();
        yAxis.Normalize();
        zAxis.Normalize();
        for (int j = 0; j < 3; ++j) {
            T(0, j) = xAxis(j);
            T(1, j) = yAxis(j);
            T(2, j) = zAxis(j);
        }

        // global wy, wz, wx
        global_w(1) = data(0) * loadFactor;
        global_w(2) = data(1) * loadFactor;
        global_w(0) = data(2) * loadFactor;

        // local
        Vector local_w(3);
        local_w.addMatrixVector(0.0, T, global_w, 1.0);

        double wx = local_w(0);  // Transverse
        double wy = local_w(1);  // Transverse
        double wz = local_w(2);  // Axial (+ve from node I to J)

        this->wx += wx;
        this->wy += wy;
        this->wz += wz;

        double Vy = 0.5 * wy * L;
        double Mz = Vy * L / 6.0;  // wy*L*L/12
        double Vz = 0.5 * wz * L;
        double My = Vz * L / 6.0;  // wz*L*L/12
        double P = wx * L;

        // Reactions in basic system
        p0[0] -= P;
        p0[1] -= Vy;
        p0[2] -= Vy;
        p0[3] -= Vz;
        p0[4] -= Vz;

        // Fixed end forces in basic system
        q0[0] -= 0.5 * P;
        if (releasez == 0) {
            q0[1] -= Mz * (2 * B1 - C1);
            q0[2] += Mz * (2 * B1 - C1);
        }
        if (releasez == 1) {
            q0[2] += wy * L * L / 8;
        }
        if (releasez == 2) {
            q0[1] -= wy * L * L / 8;
        }

        if (releasey == 0) {
            q0[3] += My * (2 * B2 - C2);
            q0[4] -= My * (2 * B2 - C2);
        }
        if (releasey == 1) {
            q0[4] -= wz * L * L / 8;
        }
        if (releasey == 2) {
            q0[3] += wz * L * L / 8;
        }

    } else {
        opserr << "Pipe::addLoad()  -- load type unknown "
                  "for element with tag: "
               << this->getTag() << endln;
        return -1;
    }

    return 0;
}

int Pipe::createPipe(int nd1, int nd2, PipeMaterial &mat,
                     PipeSection &sect, int cm, int rz, int ry,
                     double pre) {
    // nodes
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    cMass = cm;

    // section
    theSect = dynamic_cast<PipeSection *>(sect.getCopy());
    if (theSect == 0) {
        opserr << "Pipe element - failed to "
                  "get a copy of section "
               << sect.getTag() << "\n";
        return -1;
    }

    // material
    theMat = dynamic_cast<PipeMaterial *>(mat.getCopy());
    if (theMat == 0) {
        opserr << "Pipe element - failed to get a copy of "
                  "material with tag "
               << mat.getTag() << "\n";
        return -1;
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

    pressure = pre;

    return 0;
}

int Pipe::crossProduct(const Vector &A, const Vector &B,
                       Vector &res) {
    if (A.Size() != 3 || B.Size() != 3) {
        opserr << "WARNING: vector A and B's size must be 3 -- "
                  "CurvedPipe::crossProduct\n";
        return -1;
    }

    res.resize(3);
    res.Zero();

    res(0) = A(1) * B(2) - A(2) * B(1);
    res(1) = A(2) * B(0) - A(0) * B(2);
    res(2) = A(0) * B(1) - A(1) * B(0);

    return 0;
}

int Pipe::nextTransfTag() {
    ID tags = OPS_getAllCrdTransfTags();
    int gap = 10;
    if (tags.Size() == 0) {
        return gap;
    }

    int maxTag = tags(0);
    for (int i = 1; i < tags.Size(); ++i) {
        if (tags(i) < maxTag) {
            maxTag = tags(i);
        }
    }

    return maxTag + gap;
}