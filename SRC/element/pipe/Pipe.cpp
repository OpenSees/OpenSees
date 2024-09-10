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
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <LinearCrdTransf3d.h>
#include <Node.h>
#include <Pipe.h>

#include <cmath>

void *OPS_PipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr
            << "Invalid #args,  want: element pipe "
               "tag? nd1? nd2? pipeMatTag? pipeSecTag?"
               "<-T0 T0? -p p? -noThermalLoad? -noPressureLoad?>\n";
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
    bool thermalLoad = true, pressureLoad = true;
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
        } else if (strcmp(theType, "-noThermalLoad") == 0) {
            thermalLoad = false;
        } else if (strcmp(theType, "-noPressureLoad") == 0) {
            pressureLoad = false;
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

    auto *ele = new Pipe(iData[0], iData[1], iData[2], *theMat,
                         *theSect, T0, pressure, thermalLoad, pressureLoad);

    return ele;
}

Pipe::Pipe()
    : Element(0, ELE_TAG_Pipe),
      theMat(0),
      theSect(0),
      T0(0.0),
      pressure(0.0),
      thermalLoad(true),
      pressureLoad(true),
      K(12, 12),
      P(12),
      Q(12),
      kb(12, 12),
      q(6),
      wx(0.0),
      wy(0.0),
      wz(0.0),
      connectedExternalNodes(2),
      theCoordTransf(0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    for (int i = 0; i < 5; ++i) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }
}

Pipe::Pipe(int tag, int nd1, int nd2, PipeMaterial &mat,
           PipeSection &sect, double t0, double pre, bool tload,
           bool pload)
    : Element(tag, ELE_TAG_Pipe),
      theMat(0),
      theSect(0),
      T0(t0),
      pressure(pre),
      thermalLoad(tload),
      pressureLoad(pload),
      K(12, 12),
      P(12),
      Q(12),
      kb(12, 12),
      q(6),
      wx(0.0),
      wy(0.0),
      wz(0.0),
      connectedExternalNodes(2),
      theCoordTransf(0) {
    // nodes
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    theNodes[0] = 0;
    theNodes[1] = 0;

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

    // q0 - pb0, q - pb, p0 - plw
    for (int i = 0; i < 5; ++i) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }
}

Pipe::~Pipe() {
    if (theMat != 0) {
        delete theMat;
    }
    if (theSect != 0) {
        delete theSect;
    }
    if (theCoordTransf != 0) {
        delete theCoordTransf;
    }
}

const char *Pipe::getClassType(void) const { return "Pipe"; };

int Pipe::getNumExternalNodes() const { return 2; }

const ID &Pipe::getExternalNodes() { return connectedExternalNodes; }

Node **Pipe::getNodePtrs() { return theNodes; }

int Pipe::getNumDOF() { return 12; }

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

    // check DOFs
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 6) {
        opserr << "ElasticBeam3d::setDomain  tag: " << this->getTag()
               << " -- Node 1: " << connectedExternalNodes(0)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    if (dofNd2 != 6) {
        opserr << "ElasticBeam3d::setDomain  tag: " << this->getTag()
               << " -- Node 2: " << connectedExternalNodes(1)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    // copy domain
    this->DomainComponent::setDomain(theDomain);

    // create transf
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

    if (theCoordTransf != 0) {
        delete theCoordTransf;
        theCoordTransf = 0;
    }
    theCoordTransf = new LinearCrdTransf3d(nextTransfTag(), zAxis);
    if (theCoordTransf == 0) {
        opserr << "WARNING: failed to crete Transformation object -- "
                  "Pipe\n";
        exit(-1);
    }

    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
        opserr
            << "WARNING: Pipe::setDomain  tag: " << this->getTag()
            << " -- Error initializing coordinate transformation\n";
        exit(-1);
    }

    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
        opserr << "WARNING: Pipe::setDomain  tag: " << this->getTag()
               << " -- Element has zero length\n";
        exit(-1);
    }

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return;
    }
}

int Pipe::commitState() {
    int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr
            << "ElasticBeam3d::commitState () - failed in base class";
    }
    if (theCoordTransf != 0) {
        retVal += theCoordTransf->commitState();
    }
    return retVal;
}

int Pipe::revertToLastCommit() {
    int retVal = 0;
    if (theCoordTransf != 0) {
        retVal += theCoordTransf->revertToLastCommit();
    }
    return retVal;
}

int Pipe::revertToStart() {
    int retVal = 0;
    if (theCoordTransf != 0) {
        retVal += theCoordTransf->revertToStart();
    }
    return retVal;
}

int Pipe::update() {
    if (theCoordTransf != 0) {
        return theCoordTransf->update();
    }
    return 0;
}

const Matrix &Pipe::getTangentStiff() {
    K.Zero();
    if (theCoordTransf == 0) {
        return K;
    }

    const Vector &v = theCoordTransf->getBasicTrialDisp();

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::getInitialStiff failed to update section "
                  "data\n";
        return K;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::getInitialStiff failed to update material "
                  "data\n";
        return K;
    }

    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double EoverL = E * oneOverL;
    double EAoverL = A * EoverL;         // EA/L
    double GJoverL = G * Jx * oneOverL;  // GJ/L
    double B1, B2, C1, C2;               // shear coefficients
    shearCoefficients(B1, B2, C1, C2);

    q(0) = EAoverL * v(0);
    q(5) = GJoverL * v(5);
    kb.Zero();
    kb(0, 0) = EAoverL;
    kb(5, 5) = GJoverL;

    double EIzoverL2 = 2.0 * Iz * EoverL;  // 2EIz/L
    double EIzoverL4 = 2.0 * EIzoverL2;    // 4EIz/L
    q(1) = B1 * EIzoverL4 * v(1) + C1 * EIzoverL2 * v(2);
    q(2) = C1 * EIzoverL2 * v(1) + B1 * EIzoverL4 * v(2);
    kb(1, 1) = kb(2, 2) = B1 * EIzoverL4;
    kb(2, 1) = kb(1, 2) = C1 * EIzoverL2;

    double EIyoverL2 = 2.0 * Iy * EoverL;  // 2EIy/L
    double EIyoverL4 = 2.0 * EIyoverL2;    // 4EIy/L
    q(3) = B2 * EIyoverL4 * v(3) + C2 * EIyoverL2 * v(4);
    q(4) = C2 * EIyoverL2 * v(3) + B2 * EIyoverL4 * v(4);
    kb(3, 3) = kb(4, 4) = B2 * EIyoverL4;
    kb(4, 3) = kb(3, 4) = C2 * EIyoverL2;

    q(0) += q0[0];
    q(1) += q0[1];
    q(2) += q0[2];
    q(3) += q0[3];
    q(4) += q0[4];

    return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &Pipe::getInitialStiff() {
    K.Zero();
    if (theCoordTransf == 0) {
        return K;
    }

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::getInitialStiff failed to update section "
                  "data\n";
        return K;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::getInitialStiff failed to update material "
                  "data\n";
        return K;
    }

    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double EoverL = E * oneOverL;
    double EAoverL = A * EoverL;         // EA/L
    double GJoverL = G * Jx * oneOverL;  // GJ/L
    double B1, B2, C1, C2;               // shear coefficients
    shearCoefficients(B1, B2, C1, C2);

    kb.Zero();
    kb(0, 0) = EAoverL;
    kb(5, 5) = GJoverL;
    double EIzoverL2 = 2.0 * Iz * EoverL;  // 2EIz/L
    double EIzoverL4 = 2.0 * EIzoverL2;    // 4EIz/L
    kb(1, 1) = kb(2, 2) = B1 * EIzoverL4;
    kb(2, 1) = kb(1, 2) = C1 * EIzoverL2;

    double EIyoverL2 = 2.0 * Iy * EoverL;  // 2EIy/L
    double EIyoverL4 = 2.0 * EIyoverL2;    // 4EIy/L
    kb(3, 3) = kb(4, 4) = B2 * EIyoverL4;
    kb(4, 3) = kb(3, 4) = C2 * EIyoverL2;

    return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}
const Matrix &Pipe::getMass() {
    K.Zero();
    if (theCoordTransf == 0) {
        return K;
    }

    double rho = theSect->RHO();

    if (rho <= 0.0) {
        return K;
    }

    // get initial element length
    double L = theCoordTransf->getInitialLength();

    // lumped mass matrix
    double m = 0.5 * rho * L;
    K(0, 0) = m;
    K(1, 1) = m;
    K(2, 2) = m;
    K(6, 6) = m;
    K(7, 7) = m;
    K(8, 8) = m;

    return K;
}

void Pipe::initLoad() {
    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::zeroLoad failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::zeroLoad failed to update material data\n";
        return;
    }

    // zero pb0(q0) and plw(p0)
    Q.Zero();
    for (int i = 0; i < 5; ++i) {
        q0[i] = 0.0;
        p0[i] = 0.0;
    }

    wx = 0.0;
    wy = 0.0;
    wz = 0.0;
}

void Pipe::zeroLoad(void) {
    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::zeroLoad failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::zeroLoad failed to update material data\n";
        return;
    }

    // init load
    initLoad();

    // q0 due to thermal
    double temp = aveTemp();
    if (temp > 0 && thermalLoad) {
        q0[0] -= E * A * alp * temp;
    }
    // q0 due to internal pressure
    if (pressure != 0 && pressureLoad) {
        double dout = theSect->DOUT();
        double thk = theSect->WALL();

        q0[0] -=
            0.25 * pressure * (dout - thk) * (1. - 2 * nu) * A / thk;
    }
}

int Pipe::addLoad(ElementalLoad *theLoad, double loadFactor) {
    if (theCoordTransf == 0) {
        return 0;
    }

    int type;
    const Vector &data = theLoad->getData(type, loadFactor);
    double L = theCoordTransf->getInitialLength();

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
        double Pw = wx * L;

        // Reactions in basic system
        p0[0] -= Pw;  // plw1
        p0[1] -= Vy;  // plw2
        p0[2] -= Vy;  // plw8
        p0[3] -= Vz;  // plw3
        p0[4] -= Vz;  // plw9

        // Fixed end forces in basic system
        q0[0] -= 0.5 * Pw;  // pb0_1
        q0[1] -= Mz;        // pb0_2
        q0[2] += Mz;        // pb0_3
        q0[3] += My;        // pb0_4
        q0[4] -= My;        // pb0_5

    } else {
        opserr << "Pipe::addLoad()  -- load type unknown "
                  "for element with tag: "
               << this->getTag() << "\n";
        return -1;
    }

    return 0;
}

int Pipe::addInertiaLoadToUnbalance(const Vector &accel) {
    // rho
    double rho = theSect->RHO();
    if (rho <= 0.0) {
        return 0;
    }
    if (theCoordTransf == 0) {
        return 0;
    }

    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr << "Pipe::addInertiaLoadToUnbalance matrix "
                  "and vector sizes are incompatible\n";
        return -1;
    }

    // want to add ( - fact * M R * accel ) to unbalance
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * rho * L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);

    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);

    return 0;
}

const Vector &Pipe::getResistingForce() {
    P.Zero();
    if (theCoordTransf == 0) {
        return P;
    }

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::getResistingForce failed to update section "
                  "data\n";
        return P;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::getResistingForce failed to update material "
                  "data\n";
        return P;
    }

    const Vector &v = theCoordTransf->getBasicTrialDisp();

    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double EoverL = E * oneOverL;
    double EAoverL = A * EoverL;         // EA/L
    double GJoverL = G * Jx * oneOverL;  // GJ/L
    double B1, B2, C1, C2;               // shear coefficients
    shearCoefficients(B1, B2, C1, C2);

    q(0) = EAoverL * v(0);
    q(5) = GJoverL * v(5);

    double EIzoverL2 = 2.0 * Iz * EoverL;  // 2EIz/L
    double EIzoverL4 = 2.0 * EIzoverL2;    // 4EIz/L
    q(1) = B1 * EIzoverL4 * v(1) + C1 * EIzoverL2 * v(2);
    q(2) = C1 * EIzoverL2 * v(1) + B1 * EIzoverL4 * v(2);

    double EIyoverL2 = 2.0 * Iy * EoverL;  // 2EIy/L
    double EIyoverL4 = 2.0 * EIyoverL2;    // 4EIy/L
    q(3) = B2 * EIyoverL4 * v(3) + C2 * EIyoverL2 * v(4);
    q(4) = C2 * EIyoverL2 * v(3) + B2 * EIyoverL4 * v(4);

    q(0) += q0[0];
    q(1) += q0[1];
    q(2) += q0[2];
    q(3) += q0[3];
    q(4) += q0[4];

    Vector p0Vec(p0, 5);
    P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

    // subtract external load P = P - Q
    if (rho > 0) {
        P.addVector(1.0, Q, -1.0);
    }

    return P;
}
const Vector &Pipe::getResistingForceIncInertia() {
    P.Zero();
    if (theCoordTransf == 0) {
        return P;
    }

    // static force
    P = this->getResistingForce();

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 ||
        betaKc != 0.0) {
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }

    // add inertia forces from element mass
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    double L = theCoordTransf->getInitialLength();
    double rho = theSect->RHO();
    if (rho > 0) {
        double m = 0.5 * rho * L;

        P(0) += m * accel1(0);
        P(1) += m * accel1(1);
        P(2) += m * accel1(2);

        P(6) += m * accel2(0);
        P(7) += m * accel2(1);
        P(8) += m * accel2(2);
    }

    return P;
}

int Pipe::sendSelf(int commitTag, Channel &theChannel) { return 0; }

int Pipe::recvSelf(int commitTag, Channel &theChannel,
                   FEM_ObjectBroker &theBroker) {
    return 0;
}

void Pipe::Print(OPS_Stream &s, int flag) {}

int Pipe::displaySelf(Renderer &theViewer, int displayMode,
                      float fact, const char **modes, int numModes) {
    return 0;
}

Response *Pipe::setResponse(const char **argv, int argc,
                            OPS_Stream &output) {
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "Pipe");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    // global forces
    if (strcmp(argv[0], "force") == 0 ||
        strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 ||
        strcmp(argv[0], "globalForces") == 0) {
        output.tag("ResponseType", "Px_1");
        output.tag("ResponseType", "Py_1");
        output.tag("ResponseType", "Pz_1");
        output.tag("ResponseType", "Mx_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Px_2");
        output.tag("ResponseType", "Py_2");
        output.tag("ResponseType", "Pz_2");
        output.tag("ResponseType", "Mx_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 2, P);

    } else if (strcmp(argv[0], "localForce") == 0 ||
               strcmp(argv[0], "localForces") == 0) {
        output.tag("ResponseType", "N_1");
        output.tag("ResponseType", "Vy_1");
        output.tag("ResponseType", "Vz_1");
        output.tag("ResponseType", "T_1");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "N_2");
        output.tag("ResponseType", "Vy_2");
        output.tag("ResponseType", "Vz_2");
        output.tag("ResponseType", "T_2");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "Mz_2");

        theResponse = new ElementResponse(this, 3, P);

    } else if (strcmp(argv[0], "basicForce") == 0 ||
               strcmp(argv[0], "basicForces") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Mz_1");
        output.tag("ResponseType", "Mz_2");
        output.tag("ResponseType", "My_1");
        output.tag("ResponseType", "My_2");
        output.tag("ResponseType", "T");

        theResponse = new ElementResponse(this, 4, Vector(6));

    } else if (strcmp(argv[0], "basicStiffness") == 0) {
        theResponse = new ElementResponse(this, 19, Matrix(6, 6));

    } else if (strcmp(argv[0], "stiffness") == 0) {
        theResponse = new ElementResponse(this, 1, Matrix(12, 12));

    } else if (strcmp(argv[0], "deformations") == 0 ||
               strcmp(argv[0], "basicDeformations") == 0) {
        output.tag("ResponseType", "eps");
        output.tag("ResponseType", "theta11");
        output.tag("ResponseType", "theta12");
        output.tag("ResponseType", "theta21");
        output.tag("ResponseType", "theta22");
        output.tag("ResponseType", "phi");
        theResponse = new ElementResponse(this, 5, Vector(6));

    } else if (strcmp(argv[0], "sectionX") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        if (argc > 1) {
            float xL = atof(argv[1]);
            if (xL < 0.0) xL = 0.0;
            if (xL > 1.0) xL = 1.0;
            theResponse = new ElementResponse(this, 6, Vector(6));
            Information &info = theResponse->getInformation();
            info.theDouble = xL;
        }
    } else if (strcmp(argv[0], "sectionI") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        theResponse = new ElementResponse(this, 7, Vector(6));

    } else if (strcmp(argv[0], "sectionJ") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        theResponse = new ElementResponse(this, 8, Vector(6));

    } else if (strcmp(argv[0], "sectionC") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        theResponse = new ElementResponse(this, 9, Vector(6));
    }

    // ElementOutput
    output.endTag();

    // no response found, send to transformation
    if (theResponse == 0 && theCoordTransf != 0) {
        theResponse = theCoordTransf->setResponse(argv, argc, output);
    }

    return theResponse;
}

int Pipe::getResponse(int responseID, Information &info) {
    if (theCoordTransf == 0) {
        return 0;
    }

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::getResistingForce failed to update section "
                  "data\n";
        return -1;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::getResistingForce failed to update material "
                  "data\n";
        return -1;
    }

    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    switch (responseID) {
        case 1:
            // stiffness
            return info.setMatrix(this->getTangentStiff());
        case 2:
            // global forces
            return info.setVector(this->getResistingForce());
        case 3:
            // local forces
            P(0) = -q(0) + p0[0];  // -pb1+plw1
            P(1) =
                (q(1) + q(2)) * oneOverL + p0[1];  // (pb2+pb3)/L+plw2
            P(2) = -(q(3) + q(4)) * oneOverL +
                   p0[3];  // -(pb4+pb5)/L+plw3
            P(3) = -q(5);  // -pb6
            P(4) = q(3);   // pb4
            P(5) = q(1);   // pb2
            P(6) = q(0);   // pb1
            P(7) = -(q(1) + q(2)) * oneOverL +
                   p0[2];  // -(pb2+pb3)/L+plw8
            P(8) =
                (q(3) + q(4)) * oneOverL + p0[4];  // (pb4+pb5)/L+plw9
            P(9) = q(5);                           // pb6
            P(10) = q(4);                          // pb5
            P(11) = q(2);                          // pb3
            return info.setVector(P);
        case 4:
            // basic forces
            return info.setVector(q);
        case 5:
            // basic deformations
            return info.setVector(
                theCoordTransf->getBasicTrialDisp());
        case 6: {
            // section forces
            double xL = info.theDouble;
            Vector s(6);
            getSectionForce(xL, s);
            return info.setVector(s);
        }

        case 7: {
            // section forces at I
            Vector s(6);
            getSectionForce(0.0, s);
            return info.setVector(s);
        }

        case 8: {
            // section forces at J
            Vector s(6);
            getSectionForce(1.0, s);
            return info.setVector(s);
        }

        case 9: {
            // section forces at Center
            Vector s(6);
            getSectionForce(0.5, s);
            return info.setVector(s);
        }

        case 19:
            // basic stiffness
            double B1, B2, C1, C2;
            shearCoefficients(B1, B2, C1, C2);
            kb.Zero();
            kb(0, 0) = E * A / L;
            kb(1, 1) = kb(2, 2) = B1 * 4 * E * Iz / L;
            kb(1, 2) = kb(2, 1) = C1 * 2 * E * Iz / L;
            kb(3, 3) = kb(4, 4) = B2 * 4 * E * Iy / L;
            kb(3, 4) = kb(4, 3) = C2 * 2 * E * Iy / L;
            kb(5, 5) = G * Jx / L;
            return info.setMatrix(kb);

        default:
            break;
    }
}

int Pipe::setParameter(const char **argv, int argc,
                       Parameter &param) {
    if (argc < 1) {
        return -1;
    }
    return 0;
}
int Pipe::updateParameter(int parameterID, Information &info) {
    return 0;
}

double Pipe::aveTemp() {
    // get average element temperature
    double Ti = theNodes[0]->getTemp();
    double Tj = theNodes[1]->getTemp();
    double Tavg = 0.5 * (Ti + Tj);
    return Tavg - T0;
}

int Pipe::updateMaterialData(double &E, double &nu, double &G,
                             double &alp) {
    // select point based on temperature
    int retVal = 0;
    auto Tpt = theMat->selectPoint(aveTemp(), retVal);
    if (retVal < 0) {
        return retVal;
    }

    // set data
    E = Tpt.E;
    nu = Tpt.xnu;
    G = E / (2 * (1.0 + nu));
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

int Pipe::updateSectionData(double &A, double &Jx, double &Iy,
                            double &Iz, double &rho, double &alphaV) {
    A = theSect->AREA();
    Jx = theSect->JX();
    Iy = theSect->IY();
    Iz = theSect->IZ();
    rho = theSect->RHO();
    alphaV = theSect->ALFAV();
    if (alphaV > 99) {
        alphaV = 0.0;
    }

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

void Pipe::shearCoefficients(double &B1, double &B2, double &C1,
                             double &C2) {
    if (theCoordTransf == 0) {
        return;
    }
    double a1 = 1.0;
    double a2 = 1.0;
    double b1 = 1.0;
    double b2 = 1.0;

    double L = theCoordTransf->getInitialLength();

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return;
    }

    if (alphaV > 0.0 && G > 0 && E > 0 && Iz > 0) {
        double Avy = A / alphaV;
        a1 += 3.0 * E * Iz / (G * Avy * L * L);
        b1 -= 6.0 * E * Iz / (G * Avy * L * L);
    }

    if (alphaV > 0.0 && G > 0 && E > 0 && Iy > 0) {
        double Avz = A / alphaV;
        a2 += 3.0 * E * Iy / (G * Avz * L * L);
        b2 -= 6.0 * E * Iy / (G * Avz * L * L);
    }

    B1 = 3 * a1 / (4 * a1 * a1 - b1 * b1);
    C1 = 3 * b1 / (4 * a1 * a1 - b1 * b1);
    B2 = 3 * a2 / (4 * a2 * a2 - b2 * b2);
    C2 = 3 * b2 / (4 * a2 * a2 - b2 * b2);
}

void Pipe::getSectionForce(double xL, Vector &s) {
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;
    double x = xL * L;

    // section forces
    s.resize(6);
    s(0) = q(0) + wx * (L - x);                            // N(x)
    s(1) = (q(1) + q(2)) * oneOverL + wy * (x - 0.5 * L);  // Vy
    s(2) = (q(3) + q(4)) * oneOverL + wz * (0.5 * L - x);  // Vz
    s(3) = q(5);                                           // T(x)
    s(4) = (xL - 1.0) * q(3) + xL * q(4) +
           0.5 * wz * x * (L - x);  // My(x)
    s(5) = (xL - 1.0) * q(1) + xL * q(2) +
           0.5 * wy * x * (x - L);  // Mz(x)
}