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
#include "CurvedPipe.h"

#include <Domain.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Node.h>

#include <cmath>

// 20 gauss points in (wi, xi)
std::vector<double> CurvedPipe::gaussPts = {
    0.1527533871307258,  -0.0765265211334973, 0.1527533871307258,
    0.0765265211334973,  0.1491729864726037,  -0.2277858511416451,
    0.1491729864726037,  0.2277858511416451,  0.1420961093183820,
    -0.3737060887154195, 0.1420961093183820,  0.3737060887154195,
    0.1316886384491766,  -0.5108670019508271, 0.1316886384491766,
    0.5108670019508271,  0.1181945319615184,  -0.6360536807265150,
    0.1181945319615184,  0.6360536807265150,  0.1019301198172404,
    -0.7463319064601508, 0.1019301198172404,  0.7463319064601508,
    0.0832767415767048,  -0.8391169718222188, 0.0832767415767048,
    0.8391169718222188,  0.0626720483341091,  -0.9122344282513259,
    0.0626720483341091,  0.9122344282513259,  0.0406014298003869,
    -0.9639719272779138, 0.0406014298003869,  0.9639719272779138,
    0.0176140071391521,  -0.9931285991850949, 0.0176140071391521,
    0.9931285991850949};

// function object for

void *OPS_CurvedPipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "Invalid #args,  want: element CurvedPipe "
                  "tag? nd1? nd2? pipeMatTag? pipeSecTag?"
                  "xC? yC? zC?"
                  "<-T0 T0? -p p? -tolWall? tolWall? -TI? "
                  "-noThermalLoad? -noPressureLoad?>\n";
        return 0;
    }

    // get tag
    int iData[5];
    int numData = 5;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid integer input for curved pipe "
                  "element\n";
        return 0;
    }

    // get center
    Vector center(3);
    numData = 3;
    if (OPS_GetDoubleInput(&numData, &center(0)) < 0) {
        opserr << "WARNING invalid center or radius input for curved "
                  "pipe element\n";
        return 0;
    }

    // get data
    double T0 = 0.0, pressure = 0.0;
    bool thermalLoad = true, pressureLoad = true;
    double tolWall = 0.1;
    numData = 1;
    bool inter = false;
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
        } else if (strcmp(theType, "-tolWall") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &tolWall) < 0) {
                    opserr << "WARNING: failed to read fraction of "
                              "wall thickness "
                              "to be used for dimensional tolerance "
                              "tests\n";
                    return 0;
                }
            }
            if (tolWall < 0 || tolWall > 1) {
                opserr << "WARNING: tolWall < 0 or > 1\n";
                return 0;
            }
        } else if (strcmp(theType, "-TI") == 0) {
            inter = true;
        }
    }

    auto *theSect = dynamic_cast<PipeSection *>(
        OPS_getSectionForceDeformation(iData[4]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[4]
               << " is not found or not a curved pipe section\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[3]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[3]
               << " is not found or not a curved pipe material\n";
        return 0;
    }

    auto *ele = new CurvedPipe(
        iData[0], iData[1], iData[2], *theMat, *theSect, center, T0,
        pressure, tolWall, inter, thermalLoad, pressureLoad);

    return ele;
}

CurvedPipe::CurvedPipe()
    : Element(0, ELE_TAG_CurvedPipe),
      pipeEle(0),
      center(3),
      radius(0.0),
      theta0(0.0),
      tolWall(0.1),
      kp(1.0),
      Length(0.0),
      intersection(false),
      alg(),
      abl() {
    pipeEle = new Pipe();
}

CurvedPipe::CurvedPipe(int tag, int nd1, int nd2, PipeMaterial &mat,
                       PipeSection &sect, const Vector &c, double to,
                       double pre, double tol, bool inter, bool tload,
                       bool pload)
    : Element(tag, ELE_TAG_CurvedPipe),
      pipeEle(0),
      center(3),
      radius(0.0),
      theta0(0.0),
      tolWall(tol),
      kp(1.0),
      Length(0.0),
      intersection(inter),
      alg(),
      abl() {
    // pipe element
    pipeEle =
        new Pipe(tag, nd1, nd2, mat, sect, to, pre, tload, pload);
    if (pipeEle == 0) {
        opserr << "WARNING: failed to create an internal pipe "
                  "element -- CurvedPipe\n";
        exit(-1);
    }

    // center
    if (center.Size() != 3) {
        center.resize(3);
    }
    for (int i = 0; i < c.Size(); ++i) {
        if (i >= 3) {
            break;
        }
        center(i) = c(i);
    }
}

CurvedPipe::~CurvedPipe() {
    if (pipeEle != 0) {
        delete pipeEle;
    }
}

const char *CurvedPipe::getClassType(void) const {
    return "CurvedPipe";
};

int CurvedPipe::getNumExternalNodes() const {
    return pipeEle->getNumExternalNodes();
}

const ID &CurvedPipe::getExternalNodes() {
    return pipeEle->getExternalNodes();
}

Node **CurvedPipe::getNodePtrs() { return pipeEle->getNodePtrs(); }

int CurvedPipe::getNumDOF() { return pipeEle->getNumDOF(); }

void CurvedPipe::setDomain(Domain *theDomain) {
    if (theDomain == 0) {
        opserr << "CurvedPipe::setDomain -- Domain is null\n";
        exit(-1);
    }

    int ndm = OPS_GetNDM();
    if (ndm != 3) {
        opserr << "WARNING: pipe element must be 3D\n";
        exit(-1);
    }

    auto theNodes = pipeEle->getNodePtrs();
    const auto &connectedExternalNodes = pipeEle->getExternalNodes();
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

    if (theNodes[0] == 0) {
        opserr << "CurvedPipe::setDomain  tag: " << this->getTag()
               << " -- Node 1: " << connectedExternalNodes(0)
               << " does not exist\n";
        exit(-1);
    }

    if (theNodes[1] == 0) {
        opserr << "CurvedPipe::setDomain  tag: " << this->getTag()
               << " -- Node 2: " << connectedExternalNodes(1)
               << " does not exist\n";
        exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

    if (dofNd1 != 6) {
        opserr << "CurvedPipe::setDomain  tag: " << this->getTag()
               << " -- Node 1: " << connectedExternalNodes(0)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    if (dofNd2 != 6) {
        opserr << "CurvedPipe::setDomain  tag: " << this->getTag()
               << " -- Node 2: " << connectedExternalNodes(1)
               << " has incorrect number of DOF\n";
        exit(-1);
    }

    this->DomainComponent::setDomain(theDomain);

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (pipeEle->updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (pipeEle->updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return;
    }

    // compute theta0
    if (this->getTheta0() < 0) {
        opserr << "WARNING: failed to compute theta0\n";
        exit(-1);
    }
}

int CurvedPipe::commitState() { return pipeEle->commitState(); }

int CurvedPipe::revertToLastCommit() {
    return pipeEle->revertToLastCommit();
}

int CurvedPipe::revertToStart() { return pipeEle->revertToStart(); }

int CurvedPipe::update() { return pipeEle->update(); }

const Matrix &CurvedPipe::getTangentStiff() {
    return getInitialStiff();
}

const Matrix &CurvedPipe::getInitialStiff() {
    // kb and pb0
    Matrix kbm;
    Vector pb0;
    auto &K = pipeEle->refK();
    K.Zero();
    if (kb(kbm, pb0) < 0) {
        opserr
            << "WARNING: failed to compute kb -- getInitialStiff\n";
        return K;
    }

    Matrix kl(12, 12);
    kl.addMatrixTripleProduct(0.0, abl, kbm, 1.0);
    K.addMatrixTripleProduct(0.0, alg, kl, 1.0);

    return K;
}

const Matrix &CurvedPipe::getMass() {
    auto &K = pipeEle->refK();
    K.Zero();

    double rho = pipeEle->getSection()->RHO();
    if (rho <= 0) {
        return K;
    }

    // get half curve length
    double s = theta0 * radius;

    // lumped mass matrix
    double m = rho * s;
    K(0, 0) = m;
    K(1, 1) = m;
    K(2, 2) = m;
    K(6, 6) = m;
    K(7, 7) = m;
    K(8, 8) = m;

    return K;
}

void CurvedPipe::zeroLoad(void) { pipeEle->initLoad(); }

int CurvedPipe::addLoad(ElementalLoad *theLoad, double loadFactor) {
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
        // transform global member load to local
        Matrix T(3, 3);
        Vector global_w(3);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                T(i, j) = alg(i, j);
            }
        }

        // wy, wz, wx
        global_w(1) = data(0) * loadFactor;
        global_w(2) = data(1) * loadFactor;
        global_w(0) = data(2) * loadFactor;

        Vector local_w(3);
        local_w.addMatrixVector(0.0, T, global_w, 1.0);

        pipeEle->refwx() +=
            local_w(0);  // Axial (+ve from node I to J)
        pipeEle->refwy() += local_w(1);  // Transverse
        pipeEle->refwz() += local_w(2);  // Transverse

    } else {
        opserr << "CurvedPipe::addLoad()  -- load type unknown for "
                  "element with tag: "
               << this->getTag() << "\n";
        return -1;
    }

    return 0;
}

int CurvedPipe::addInertiaLoadToUnbalance(const Vector &accel) {
    double rho = pipeEle->getSection()->RHO();
    if (rho == 0.0) return 0;

    // get R * accel from the nodes
    auto theNodes = pipeEle->getNodePtrs();
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);

    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
        opserr << "CurvedPipe::addInertiaLoadToUnbalance matrix and "
                  "vector sizes are incompatible\n";
        return -1;
    }

    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix

    // get half curve length
    double s = theta0 * radius;

    // lumped mass matrix
    double m = rho * s;

    auto &Q = pipeEle->refQ();
    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);

    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);

    return 0;
}

const Vector &CurvedPipe::getResistingForceIncInertia() {
    auto &P = pipeEle->refP();
    P = this->getResistingForce();

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 ||
        betaKc != 0.0) {
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }

    double rho = pipeEle->getSection()->RHO();
    if (rho == 0.0) return P;

    // add inertia forces from element mass
    auto theNodes = pipeEle->getNodePtrs();
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    // take advantage of lumped mass matrix
    // get half curve length
    double s = theta0 * radius;

    // lumped mass matrix
    double m = rho * s;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);

    P(6) += m * accel2(0);
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);

    return P;
}

const Vector &CurvedPipe::getResistingForce() {
    auto &P = pipeEle->refP();

    // kb and pb0
    Matrix kbm;
    Vector pb0;
    P.Zero();
    if (kb(kbm, pb0) < 0) {
        opserr << "WARNING: failed to compute kb -- "
                  "getResistingForce\n ";
        return P;
    }

    // ug
    Vector ug(12);
    auto theNodes = pipeEle->getNodePtrs();
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();
    for (int i = 0; i < 6; i++) {
        ug[i] = disp1(i);
        ug[i + 6] = disp2(i);
    }

    // ul
    Vector ul(12);
    ul.addMatrixVector(0.0, alg, ug, 1.0);

    // ub
    Vector ub(6);
    ub.addMatrixVector(0.0, abl, ul, 1.0);

    // vector pb
    Vector pb(6);
    pb.addMatrixVector(0.0, kbm, ub, 1.0);
    pb += pb0;

    // pl
    Vector pl(12);
    pl.addMatrixTransposeVector(0.0, abl, pb, 1.0);

    // plw
    Vector plwm(12);
    plw(plwm);
    pl += plwm;

    // P
    P.addMatrixTransposeVector(0.0, alg, pl, 1.0);

    // subtract external load P = P - Q
    double rho = pipeEle->getSection()->RHO();
    if (rho != 0) {
        auto &Q = pipeEle->refQ();
        P.addVector(1.0, Q, -1.0);
    }

    return P;
}

int CurvedPipe::sendSelf(int cTag, Channel &theChannel) { return 0; }

int CurvedPipe::recvSelf(int cTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker) {
    return 0;
}

void CurvedPipe::Print(OPS_Stream &s, int flag) {}

int CurvedPipe::displaySelf(Renderer &theViewer, int displayMode,
                            float fact, const char **modes,
                            int numMode) {
    return 0;
}

Response *CurvedPipe::setResponse(const char **argv, int argc,
                                  OPS_Stream &output) {
    Response *theResponse = 0;
    output.tag("ElementOutput");
    output.attr("eleType", "CurvedPipe");
    output.attr("eleTag", this->getTag());
    auto &nodes = pipeEle->getExternalNodes();
    output.attr("node1", nodes(0));
    output.attr("node2", nodes(1));

    if (strcmp(argv[0], "sectionX") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        if (argc > 1) {
            float xL = atof(argv[1]);
            if (xL < -1.0) xL = -.0;
            if (xL > 1.0) xL = 1.0;
            theResponse = new ElementResponse(this, 1, Vector(6));
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
        theResponse = new ElementResponse(this, 2, Vector(6));

    } else if (strcmp(argv[0], "sectionJ") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        theResponse = new ElementResponse(this, 3, Vector(6));

    } else if (strcmp(argv[0], "sectionC") == 0) {
        output.tag("ResponseType", "N");
        output.tag("ResponseType", "Vy");
        output.tag("ResponseType", "Vz");
        output.tag("ResponseType", "T");
        output.tag("ResponseType", "My");
        output.tag("ResponseType", "Mz");
        theResponse = new ElementResponse(this, 4, Vector(6));
    } else if (strncmp(argv[0], "center", 80) == 0) {
      output.tag("ResponseType", "XC");
      output.tag("ResponseType", "YC");
      output.tag("ResponseType", "ZC");      
      theResponse = new ElementResponse(this, 5, Vector(3));
    }
    return theResponse;
}

int CurvedPipe::getResponse(int responseID, Information &info) {
    if (responseID == 1) {
        // section forces
        double xL = info.theDouble;
        double theta = xL * theta0;
        Vector s(6);

        if (getSectionForce(theta, s) < 0) {
            opserr << "WARNING: failed to get section force\n";
            return -1;
        }
        return info.setVector(s);
    } else if (responseID == 2) {
        // section forces at I
        Vector s(6);
        if (getSectionForce(-theta0, s) < 0) {
            opserr << "WARNING: failed to get section force at I\n";
            return -1;
        }
        return info.setVector(s);

    } else if (responseID == 3) {
        // section forces at J
        Vector s(6);
        if (getSectionForce(theta0, s) < 0) {
            opserr << "WARNING: failed to get section force at J\n";
            return -1;
        }
        return info.setVector(s);

    } else if (responseID == 4) {
        // section forces at C
        Vector s(6);
        if (getSectionForce(0.0, s) < 0) {
            opserr
                << "WARNING: failed to get section force at center\n";
            return -1;
        }
        return info.setVector(s);
    } else if (responseID == 5) {
      return info.setVector(center);
    }
    return 0;
}

int CurvedPipe::setParameter(const char **argv, int argc,
                             Parameter &param) {
    if (argc < 1) return -1;

    return 0;
}

int CurvedPipe::updateParameter(int parameterID, Information &info) {
    return 0;
}

int CurvedPipe::getTheta0() {
    // nodes
    auto theNodes = pipeEle->getNodePtrs();
    const auto &crdsI = theNodes[0]->getCrds();
    const auto &crdsJ = theNodes[1]->getCrds();
    double thk = pipeEle->getSection()->WALL();

    // check intersection point
    if (intersection) {
        // get T
        Vector crdsT = center;

        // I->T
        Vector IT = crdsT;

        // J->T
        Vector JT = crdsT;

        // out of plane direction
        Vector N(crdsT.Size());

        IT -= crdsI;
        JT -= crdsJ;
        if (pipeEle->crossProduct(IT, JT, N) < 0) {
            return -1;
        }
        N.Normalize();

        if (N.Norm() < tolWall * thk) {
            opserr << "The intersection point is too close to the "
                      "chord\n";
            return -1;
        }

        // system to solve for center
        Matrix tempC(center.Size(), center.Size());
        for (int j = 0; j < center.Size(); ++j) {
            tempC(0, j) = IT(j);
            tempC(1, j) = JT(j);
            tempC(2, j) = N(j);
        }

        Vector rhs(center.Size());
        rhs(0) = IT ^ crdsI;
        rhs(1) = JT ^ crdsJ;
        rhs(2) = N ^ crdsT;

        if (tempC.Solve(rhs, center) < 0) {
            opserr << "WANRING: failed to compute center from "
                      "tangent intersection\n";
            return -1;
        }

        intersection = false;
    }

    // transformation matrix alg
    Vector yAxis = crdsI;
    yAxis += crdsJ;
    yAxis *= 0.5;
    yAxis = center - yAxis;
    yAxis.Normalize();

    Vector xAxis = crdsJ;
    xAxis -= crdsI;
    Length = xAxis.Norm();
    if (Length < 1e-16) {
        opserr << "WARNING: element length < 1e-16\n";
        return -1;
    }
    xAxis.Normalize();

    Vector zAxis;
    if (pipeEle->crossProduct(xAxis, yAxis, zAxis) < 0) {
        return -1;
    }

    alg.resize(12, 12);
    alg.Zero();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            alg(3 * i, 3 * i + j) = xAxis(j);
            alg(3 * i + 1, 3 * i + j) = yAxis(j);
            alg(3 * i + 2, 3 * i + j) = zAxis(j);
        }
    }

    // abl
    double invL = 1.0 / Length;
    abl.resize(6, 12);
    abl.Zero();
    abl(0, 0) = -1;
    abl(0, 6) = 1;
    abl(1, 1) = invL;
    abl(1, 5) = 1;
    abl(1, 7) = -invL;
    abl(2, 1) = invL;
    abl(2, 7) = -invL;
    abl(2, 11) = 1;
    abl(3, 2) = -invL;
    abl(3, 4) = 1;
    abl(3, 8) = invL;
    abl(4, 2) = -invL;
    abl(4, 8) = invL;
    abl(4, 10) = 1;
    abl(5, 3) = -1;
    abl(5, 9) = 1;

    // compute radius
    double R1 = (center - crdsI).Norm();
    double R2 = (center - crdsJ).Norm();
    radius = (R1 + R2) / 2.0;
    if (radius <= 0) {
        opserr << "WARNING: radius <= 0\n";
        return -1;
    }

    if (fabs(R1 - R2) > tolWall * thk) {
        opserr << "WARNING: the computed radius from node I is "
                  "different to "
                  "the one computed from node J more than "
               << tolWall << " * wall thickness\n";
        return -1;
    }

    // compute theta0
    double Lhalf = Length / 2.0;
    if (Lhalf > 0.9848 * radius) {
        opserr << "WARNING: the angle of the curve should be < 160 "
                  "degree\n";
        return -1;
    }
    theta0 = asin(Lhalf / radius);

    // update material data
    double E, nu, G, alp;
    if (pipeEle->updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return -1;
    }

    // compute kp
    double dout = pipeEle->getSection()->DOUT();
    double RM = (dout - thk) * 0.5;
    double h = thk * radius / (RM * RM);
    double DUM2 = pow(radius / thk, 4.0 / 3);
    double DUM = 6 * pipeEle->getPressure() / (E * h);
    DUM = 1.0 + DUM * DUM2;
    this->kp = 1.65 / h / DUM;
    if (this->kp < 1) {
        this->kp = 1.0;
    }

    return 0;
}

void CurvedPipe::bx(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();
    double c = cos(theta);
    double s = sin(theta);
    double L = Length;
    double R = radius;
    double H = R * (c - cos(theta0));
    double H0 = R * cos(theta0);
    double invL = 1.0 / L;
    mat(0, 0) = c;
    mat(0, 1) = -s * invL;
    mat(0, 2) = -s * invL;
    mat(1, 0) = -H;
    mat(1, 1) = R * s * invL - 0.5;
    mat(1, 2) = R * s * invL + 0.5;
    mat(2, 3) = s * H0 * invL - 0.5 * c;
    mat(2, 4) = s * H0 * invL + 0.5 * c;
    mat(2, 5) = -s;
    mat(3, 3) = R * invL * (1 - cos(theta - theta0));
    mat(3, 4) = R * invL * (1 - cos(theta + theta0));
    mat(3, 5) = c;
    mat(4, 0) = s;
    mat(4, 1) = c * invL;
    mat(4, 2) = c * invL;
    mat(5, 3) = invL;
    mat(5, 4) = invL;
}

void CurvedPipe::Spx(double theta, Vector &vec) {
    vec.resize(6);
    vec.Zero();
    double wx = pipeEle->refwx();
    double wy = pipeEle->refwy();
    double wz = pipeEle->refwz();
    double c = cos(theta);
    double s = sin(theta);
    double L = Length;
    double R = radius;
    double R2 = R * R;
    double H0 = R * cos(theta0);
    double invL = 1.0 / L;
    double stt0 = sin(theta - theta0);
    double ctt0 = cos(theta - theta0);
    vec(0) = wx * R *
                 (c * theta0 - s - c * theta +
                  2 * s * theta0 * H0 * invL) -
             wy * R * s * theta;
    vec(1) = wx * R *
                 (c * R * theta - 2 * s * R * theta0 * H0 * invL -
                  c * R * theta0 + theta0 * H0) +
             wy * R * (s * R * theta - theta0 * L * 0.5 + c * R - H0);
    vec(2) = wz * R2 * (-theta0 * stt0 + ctt0 - 1);
    vec(3) = wz * R2 * (-theta + theta0 * ctt0 + stt0);
    vec(4) = wx * R *
                 (c + s * theta0 - s * theta -
                  2 * c * theta0 * H0 * invL) +
             wy * R * c * theta;
    vec(5) = -wz * R * theta;
}

void CurvedPipe::fs(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();

    if (kp < 1) {
        kp = 1;
    }

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (pipeEle->updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return;
    }

    // update material data
    double E, nu, G, alp;
    if (pipeEle->updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return;
    }

    mat(0, 0) = 1.0 / (E * A);
    mat(1, 1) = kp / (E * Iz);
    mat(2, 2) = kp / (E * Iy);
    mat(3, 3) = 1.0 / (G * Jx);

    if (alphaV > 99) {
        alphaV = 0.0;
    }
    if (alphaV > 0) {
        double Ay = A / alphaV;
        double Az = A / alphaV;
        mat(4, 4) = 1.0 / (G * Ay);
        mat(5, 5) = 1.0 / (G * Az);
    }
}

void CurvedPipe::fb(double theta, Matrix &mat) {
    mat.resize(6, 6);
    mat.Zero();

    Matrix bxmat, fsmat;
    bx(theta, bxmat);
    fs(theta, fsmat);
    mat.addMatrixTripleProduct(0.0, bxmat, fsmat, 1.0);
}

void CurvedPipe::ubno(double theta, Vector &vec) {
    vec.resize(6);
    vec.Zero();

    Matrix fsmat, bxmat;
    Vector spxvec;
    fs(theta, fsmat);
    bx(theta, bxmat);
    Spx(theta, spxvec);

    Vector temp(6);
    temp.addMatrixVector(0.0, fsmat, spxvec, 1.0);
    vec.addMatrixTransposeVector(0.0, bxmat, temp, 1.0);
}

int CurvedPipe::kb(Matrix &mat, Vector &vec) {
    Matrix fbmat;
    Vector ubnovec;

    // update section data
    double A, Jx, Iy, Iz, rho, alphaV;
    if (pipeEle->updateSectionData(A, Jx, Iy, Iz, rho, alphaV) < 0) {
        opserr << "Pipe::setDomain failed to update section data\n";
        return -1;
    }

    // update material data
    double E, nu, G, alp;
    if (pipeEle->updateMaterialData(E, nu, G, alp) < 0) {
        opserr << "Pipe::setDomain failed to update material data\n";
        return -1;
    }

    // radious
    double R = radius;

    // integrate fb and ub0
    integrateGauss(-theta0, theta0, fbmat, ubnovec);

    // ub0 due to thermal
    double dT = pipeEle->aveTemp();
    if (dT > 0 && pipeEle->inclThermalLoad()) {
        ubnovec(0) += 2 * R * alp * dT * sin(theta0);
    }

    // ub0 due to internal pressure
    double pressure = pipeEle->getPressure();
    if (pressure != 0 && pipeEle->inclPressureLoad()) {
        double dout = pipeEle->getSection()->DOUT();
        double thk = pipeEle->getSection()->WALL();

        // strain due to pressure converted from SAP 5
        double RM = (dout - thk) * 0.5;

        int MEL = 1;
        double BTA = 0.0;

        if (MEL) {
            // MEL REPORT 10-66, EQUATION (3-29).
            double DUM = 3.14159265 * pow(RM, 4) * pressure;
            DUM *= 0.5 / (E * Iz);
            double DU2 = RM / R;
            BTA = DUM * (2 - 2 * nu + (3 + 1.5 * nu) * DU2 * DU2);
            BTA = -BTA / R;

        } else {
            // C. S. PARKER, EQUATION (10), 2-28-69.
            double DU2 = R / RM;
            double DUM = pressure * RM * 0.5 / (E * thk);
            double DU3 =
                1.0 + DUM * (1 - nu * (2 * DU2 - 1) / (DU2 - 1));
            BTA = DU3 / (1.0 + DUM * (2 - nu));
            BTA = -(1.0 - BTA) / R;
        }

        // axial strain
        ubnovec(0) += 0.5 * pressure * R * (dout - thk) *
                      (1. - 2 * nu) * sin(theta0) / (E * thk);

        // curvature
        ubnovec(0) +=
            2 * R * R * BTA * (theta0 * cos(theta0) - sin(theta0));
        ubnovec(1) += -R * BTA * theta0;
        ubnovec(2) += R * BTA * theta0;
    }

    // kb
    if (fbmat.Invert(mat) < 0) {
        return -1;
    }

    // pb0
    vec.resize(6);
    vec.addMatrixVector(0.0, mat, ubnovec, -1.0);

    return 0;
}

void CurvedPipe::integrateGauss(double a, double b, Matrix &resm,
                                Vector &resv) {
    double p = (b - a) / 2.0;
    double q = (b + a) / 2.0;
    resm.resize(6, 6);
    resm.Zero();
    resv.resize(6);
    resv.Zero();
    for (int i = 0; i < (int)gaussPts.size() / 2; ++i) {
        double xi = gaussPts[2 * i + 1] * p + q;
        double wi = gaussPts[2 * i] * p;
        Matrix mat;
        Vector vec;
        fb(xi, mat);
        resm.addMatrix(1.0, mat, wi);
        ubno(xi, vec);
        resv.addVector(1.0, vec, wi);
    }
    resm *= radius;
    resv *= radius;
}

void CurvedPipe::plw(Vector &vec) {
    vec.resize(12);
    vec.Zero();

    // member load
    double wx = pipeEle->refwx();
    double wy = pipeEle->refwy();
    double wz = pipeEle->refwz();
    double R = radius;
    double H0 = R * cos(theta0);
    double L = Length;

    // plw1
    vec(0) = -2 * wx * R * theta0;

    // plw2
    vec(1) = wx * R * (1.0 - 2 * theta0 * H0 / L) - wy * R * theta0;

    // plw8
    vec(7) = -wx * R * (1.0 - 2 * theta0 * H0 / L) - wy * R * theta0;

    // plw3
    vec(2) = -wz * R * theta0;

    // plw9
    vec(8) = -wz * R * theta0;

    // plw4
    vec(3) = wz * R * (L - 2 * theta0 * H0);
}

int CurvedPipe::getSectionForce(double theta, Vector &s) {
    // b
    Matrix b;
    this->bx(theta, b);

    // sp
    this->Spx(theta, s);

    // kb and pb0
    Matrix kbm;
    Vector pb0;
    if (kb(kbm, pb0) < 0) {
        opserr << "WARNING: failed to compute kb -- "
                  "getSectionForce\n ";
        return -1;
    }

    // ug
    Vector ug(12);
    auto theNodes = pipeEle->getNodePtrs();
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();
    for (int i = 0; i < 6; i++) {
        ug[i] = disp1(i);
        ug[i + 6] = disp2(i);
    }

    // ul
    Vector ul(12);
    ul.addMatrixVector(0.0, alg, ug, 1.0);

    // ub
    Vector ub(6);
    ub.addMatrixVector(0.0, abl, ul, 1.0);

    // vector pb
    Vector pb(6);
    pb.addMatrixVector(0.0, kbm, ub, 1.0);
    pb += pb0;

    // s = bpb +sp
    s.addMatrixVector(1.0, b, pb, 1.0);

    // from: s = [N, Mz, My, T, Vy, Vz]
    double Vy = s[4];
    double Vz = s[5];
    double My = s[2];
    double Mz = s[1];

    // to: s = [N, Vy, Vz, T, My, Mz]
    s[1] = Vy;
    s[2] = Vz;
    s[4] = My;
    s[5] = Mz;

    return 0;
}
