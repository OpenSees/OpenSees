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
#include <Node.h>
#include <Pipe.h>

#include <cmath>

void *OPS_PipeElement() {
    // line: 9343
    //               TAG  TP  I   J   MAT SEC TO   P  VECY     INC
    // READ (5,1040) INEL,R1,INI,INJ,IMAT,ISP,TRI,PRI,X2,X3,X4,INC
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "Invalid #args,  want: element pipe "
                  "tag? nd1? nd2? matTag? secTag? <vecyx? vecyy? "
                  "vecyz?> "
                  "<To? p? gDir?>\n";
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
    double data[] = {0, 0, 0, 0, 0};
    numData = OPS_GetNumRemainingInputArgs();
    if (numData > 5) {
        numData = 5;
    }
    if (numData > 0) {
        if (OPS_GetDoubleInput(&numData, data) < 0) {
            opserr << "WARNING: invalid data for pipe element\n";
            return 0;
        }
    }

    // get gravity direction
    int gDir = -3;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &gDir) < 0) {
            opserr << "WARNING invalid gDir for pipe element\n";
            return 0;
        }
    }

    auto *theSect = dynamic_cast<PipeSection *>(
        OPS_getSectionForceDeformation(iData[4]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[4]
               << " is not found\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[3]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[3]
               << " is not found\n";
        return 0;
    }

    auto *ele =
        new Pipe(iData[0], iData[1], iData[2], *theMat, *theSect,
                 data[0], data[1], data[2], data[3], data[4], gDir);

    return ele;
}

Pipe::Pipe()
    : Element(0, ELE_TAG_Pipe),
      theNodes(2, 0),
      connectedExternalNodes(2),
      theMat(0),
      theSect(0),
      M(),
      K(),
      P(),
      localY(3),
      T0(0.0),
      pressure(0.0),
      gDir(-3),
      currLN(0.0),
      currTavg(0.0),
      currT() {}

Pipe::Pipe(int tag, int nd1, int nd2, PipeMaterial &mat,
           PipeSection &sect, double vecyx, double vecyy,
           double vecyz, double t0, double pre, int g)
    : Element(tag, ELE_TAG_Pipe),
      theNodes(2),
      connectedExternalNodes(2),
      theMat(0),
      theSect(0),
      M(),
      K(),
      P(),
      localY(3),
      T0(t0),
      pressure(pre),
      gDir(g),
      currLN(0.0),
      currTavg(0.0),
      currT() {
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

    // local y
    localY(0) = vecyx;
    localY(1) = vecyy;
    localY(2) = vecyz;
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

int Pipe::getNumExternalNodes(void) const { return 2; }

const ID &Pipe::getExternalNodes(void) {
    return connectedExternalNodes;
}

Node **Pipe::getNodePtrs(void) { return &theNodes[0]; }

int Pipe::getNumDOF(void) { return 12; }

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

    // set data
    M.resize(getNumDOF(), getNumDOF());
    M.Zero();
    K.resize(getNumDOF(), getNumDOF());
    K.Zero();
    P.resize(getNumDOF());
    P.Zero();
    currT.resize(ndm, ndm);
    currT.Zero();
    currLN = 0;
    currTavg = 0;

    // set domain
    this->DomainComponent::setDomain(theDomain);
}

int Pipe::commitState(void) {
    int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "Pipe::commitState () - failed in "
                  "base class";
    }
    retVal += theMat->commitState();
    retVal += theSect->commitState();

    M.resize(getNumDOF(), getNumDOF());
    M.Zero();

    K.resize(getNumDOF(), getNumDOF());
    K.Zero();

    P.resize(getNumDOF());
    P.Zero();

    int ndm = OPS_GetNDM();
    currT.resize(ndm, ndm);
    currT.Zero();
    currLN = 0;
    currTavg = 0;

    return retVal;
}

int Pipe::revertToLastCommit(void) {
    int retVal = 0;
    retVal += theMat->revertToLastCommit();
    retVal += theSect->revertToLastCommit();
    return retVal;
}

int Pipe::revertToStart(void) {
    int retVal = 0;
    retVal += theMat->revertToStart();
    retVal += theSect->revertToStart();
    return retVal;
}

int Pipe::update(void) {
    int retVal = 0;

    // line 9460
    // transformation matrix
    Matrix T;
    retVal = tangdc(T);
    if (retVal < 0) {
        return retVal;
    }

    // line 9462-9466
    // get average element temperature
    double Ti = theNodes[0]->getTemp();
    double Tj = theNodes[1]->getTemp();
    double Tavg = 0.5 * (Ti + Tj);
    auto Tpt = theMat->selectPoint(Tavg, retVal);
    if (retVal < 0) {
        return retVal;
    }

    // line 9474-9484
    // test if new K is needed
    const auto &crds1 = theNodes[0]->getCrds();
    const auto &crds2 = theNodes[1]->getCrds();
    auto localX = crds2 - crds1;
    double xln = localX.Norm();
    if (fabs(currLN - xln) + fabs(currTavg - Tavg) < 1e-6) {
        // check transformation matrix
        double du2 = 0.0;
        double du3 = 0.0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                du2 += fabs(currT(i, j) - T(i, j));
                du3 += fabs(T(i, j) * 1e-6);
            }
        }

        if (du2 <= du3) {
            // no change of the stiffness matrix
            return retVal;
        }
    }

    // update curr
    currLN = xln;
    currTavg = Tavg;
    currT = T;

    // compute matrices
    Matrix Flex, invF, B, H, FEF, FEFC, SA;
    retVal = tangks(Tpt, Flex, invF, B, H, FEF, FEFC, SA);
    if (retVal < 0) {
        return retVal;
    }

    // 9521
    double NS = 12;
    double delt = currTavg - T0;
    double wgt = theSect->WGT();

    // local to global for K
    retVal = localGlobal();
    if (retVal < 0) {
        return retVal;
    }

    // TODO: 9682

    return retVal;
}

const Matrix &Pipe::getTangentStiff(void) { return K; }
const Matrix &Pipe::getInitialStiff(void) { return K; }
const Matrix &Pipe::getMass(void) { return M; }

void Pipe::zeroLoad(void) { P.Zero(); }

int Pipe::addLoad(ElementalLoad *theLoad, double loadFactor) {
    return 0;
}

int Pipe::addInertiaLoadToUnbalance(const Vector &accel) { return 0; }

const Vector &Pipe::getResistingForce(void) { return P; }
const Vector &Pipe::getResistingForceIncInertia(void) { return P; }

int Pipe::sendSelf(int commitTag, Channel &theChannel) { return 0; }
int Pipe::recvSelf(int commitTag, Channel &theChannel,
                   FEM_ObjectBroker &theBroker) {
    return 0;
}

void Pipe::Print(OPS_Stream &s, int flag) {}

Response *Pipe::setResponse(const char **argv, int argc,
                            OPS_Stream &s) {
    return 0;
}

int Pipe::getResponse(int responseID, Information &info) { return 0; }

int Pipe::setParameter(const char **argv, int argc,
                       Parameter &param) {
    return -1;
}
int Pipe::updateParameter(int parameterID, Information &info) {
    return -1;
}

// utility
// default convention for local y-axis
void Pipe::defaultLocalY(const Vector &localX, Vector &deLocalY) {
    // line 17356-17378
    // default convention
    double du2 = localX(0) * localX(0) + localX(2) * localX(2);
    if (du2 <= 1e-12) {
        // vertical member
        deLocalY(0) = 0.0;
        deLocalY(1) = 0.0;
        deLocalY(2) = 1.0;
    } else {
        // non-vertical member
        du2 = sqrt(du2);
        deLocalY(0) = -localX(0) * localX(1) / du2;
        deLocalY(1) = du2;
        deLocalY(2) = -localX(2) * localX(1) / du2;
    }
}

// computation of direction cosine array for the local axes of
// pipe tangent element
int Pipe::tangdc(Matrix &dircos) {
    // line 17296-17387
    int ndm = OPS_GetNDM();
    dircos.resize(ndm, ndm);
    dircos.Zero();

    // local x-axis from node I to node J
    const auto &crds1 = theNodes[0]->getCrds();
    const auto &crds2 = theNodes[1]->getCrds();
    auto localX = crds2 - crds1;
    double xln = localX.Norm();
    if (xln <= 1e-8) {
        opserr << "WARNING: element length < 1e-8\n";
        return -1;
    }
    localX /= xln;

    // local y-axis
    if (localY.Norm() < 1e-8) {
        // get default local y
        defaultLocalY(localX, localY);
    } else {
        // direct user input of the local y-axis
        localY.Normalize();

        // test if input local y-axis is vertical to local x-axis
        double dotxy = 0.0;
        for (int i = 0; i < ndm; ++i) {
            dotxy += localX(i) * localY(i);
        }
        if (fabs(dotxy) >= 1e-6) {
            // get default local y
            defaultLocalY(localX, localY);
        }
    }

    // local z-axis
    Vector localZ(ndm);
    localZ(0) = localX(1) * localY(2) - localX(2) * localY(1);
    localZ(1) = localX(2) * localY(0) - localX(0) * localY(2);
    localZ(2) = localX(0) * localY(1) - localX(1) * localY(0);

    // copy to the matrix
    dircos.resize(ndm, ndm);
    for (int i = 0; i < ndm; ++i) {
        dircos(i, 0) = localX(i);
        dircos(i, 1) = localY(i);
        dircos(i, 2) = localZ(i);
    }

    return 0;
}

// computation of the element stiffness and load matrices
// for a tangent (straight) pipe element
int Pipe::tangks(const PipeMaterialTemperaturePoint &Tpt,
                 Matrix &Flex, Matrix &invF, Matrix &B, Matrix &H,
                 Matrix &FEF, Matrix &FEFC, Matrix &SA) {
    // 17395-17742

    // set the factor for axial deformations
    double axial = 1.0;

    // set the factor for shear deformations (=0, neglect)
    double xkap = theSect->ALFAV();
    if (xkap > 99.99) {
        xkap = 0.0;
    }

    // compute the material factors
    double re = Tpt.E;
    double xnu1 = Tpt.xnu + 1.0;

    // compute section property constants
    double ra = axial * currLN * re / theSect->AREA();
    double rv = 2 * xkap * xnu1 * currLN * re / theSect->AREA();
    double rt = xnu1 * currLN * re / theSect->XMI();
    double rb = currLN * re / theSect->XMI();
    double xl2 = currLN * currLN;

    // form the node flexibility matrix at node J
    // referenced to the local (x,y,z) coordinate system at node I
    // x - direction: axial from node I to node J
    // y - direction: transverse bending axis
    // z - direction: transverse bending axis orthogonal to x and y
    Flex.resize(6, 6);
    Flex.Zero();

    // axial for Temp matrix
    Flex(0, 0) += ra;

    // shear for flex matrix
    Flex(1, 1) += rv;
    Flex(2, 2) += rv;

    // torsion for flex matrix
    Flex(3, 3) += rt;

    // bending for flex matrix
    Flex(1, 1) += rb * xl2 / 3.0;
    Flex(2, 2) += rb * xl2 / 3.0;
    Flex(4, 4) += rb;
    Flex(5, 5) += rb;
    Flex(1, 5) += rb * currLN * 0.5;
    Flex(2, 4) -= rb * currLN * 0.5;

    // fill flex matrix
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            Flex(j, i) = Flex(i, j);
        }
    }

    // form the node stiffness matrix by inverting Flex
    invF.resize(6, 6);
    invF.Zero();
    Flex.Invert(invF);

    // compute the deflections/rotations (measured in
    // the x,y,z system at node I) at node J to uniform loads
    // in each of the x,y,z directions (at I).
    // the uniform loads are direction invariant
    // with position along the length, and node I
    // is completely fixed while node J is free.

    // free end deflections at node J due to
    // 1. uniform load in x direction
    // 2. uniform load in y direction
    // 3. uniform load in z direction
    // 4. uniform thermal expansion
    // 5. internal pressure
    B.resize(6, 6);
    B.Zero();

    // axial for B
    ra *= 0.5 * currLN;
    B(0, 0) += ra;

    // shear for B
    rv *= 0.5 * currLN;
    B(1, 1) += rv;
    B(2, 2) += rv;

    // bending for B
    rb *= xl2 / 6.0;
    B(1, 1) += rb * currLN * 0.75;
    B(2, 2) += rb * currLN * 0.75;
    B(4, 2) -= rb;
    B(5, 1) += rb;

    // compute the free node deflections at end J due to a uniform
    // thermal expansion
    B(0, 3) = currLN * Tpt.alp;

    // 17577: compute the free node deflections at end J due to
    // pressure
    double rm = (theSect->DOUT() - theSect->WALL()) * 0.5;
    B(0, 4) = 0.5 * pressure * rm * re / theSect->WALL() *
              (1.0 - 2.0 * Tpt.xnu) * currLN;

    // 17585: set up the force transformation relating reactions at
    // node I acting on the member end due to unit loads applied to
    // the member end at node J

    // force transformation relating reactions at node I
    // due to unit loads at node J
    H.resize(6, 6);
    H.Zero();
    for (int i = 0; i < 6; ++i) {
        H(i, i) = -1.0;
    }
    H(5, 1) = -currLN;
    H(4, 2) = currLN;

    // 17601: form the upper triangular portion of the local element
    // stiffness matrix for the tangent

    // local stiffness matrix
    K.resize(getNumDOF(), getNumDOF());
    K.Zero();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            // 17606
            K(i + 6, j + 6) = invF(i, j);
        }
    }

    for (int IR = 0; IR < 6; ++IR) {
        for (int IC = 0; IC < 6; ++IC) {
            // 17611
            K(IR, IC + 6) = 0.0;

            // 17612
            for (int IN = 0; IN < 6; ++IN) {
                // 17613
                K(IR, IC + 6) += H(IR, IN) * invF(IN, IC);
            }
        }
    }

    for (int IR = 0; IR < 6; ++IR) {
        for (int IC = IR; IC < 6; ++IC) {
            // 17619
            K(IR, IC) = 0.0;

            // 17620
            for (int IN = 0; IN < 6; ++IN) {
                // 17621
                K(IR, IC) += K(IR, IN + 6) * H(IC, IN);
            }
        }
    }

    // 17625: reflect for symmetry
    for (int i = 0; i < getNumDOF(); ++i) {
        for (int k = 0; k < getNumDOF(); ++k) {
            K(k, i) = K(i, k);
        }
    }

    // 17639: compute the restrained node forces
    // acting on the nodes of the tangent due
    // to the member loadings

    // fixed end forces (acting on the nodes) due to
    // 1. uniform load in x direction
    // 2. uniform load in y direction
    // 3. uniform load in z direction
    // 4. uniform thermal expansion
    // 5. internal pressure
    FEF.resize(12, 5);
    FEF.Zero();
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 12; ++j) {
            for (int k = 0; k < 6; ++k) {
                // 17646
                FEF(j, i) -= K(j, k + 6) * B(k, i);
            }
        }
    }

    // 17650: for the distributed loads superimpose
    // the cantilever reactions acting on the element at node I
    double dum = 0.5 * xl2;
    FEF(0, 0) -= currLN;
    FEF(1, 1) -= currLN;
    FEF(5, 1) -= dum;
    FEF(2, 2) -= currLN;
    FEF(4, 2) += dum;

    // 17668: form the lumped mass matrix
    dum = 0.5 * currLN * theSect->RHO();
    M.resize(getNumDOF(), getNumDOF());
    M.Zero();
    for (int k = 0; k < 3; ++k) {
        M(k, k) = dum;
        M(k + 6, k + 6) = dum;
    }

    // 17678: compute the fixe-node corrections
    // to the section stresses due to element loads.
    // forces act on the segment between the point
    // of evaluation and node I

    // fixed-end force corrections to the member loads due to
    // the five types of element loads
    FEFC.resize(18, 5);
    FEFC.Zero();

    // 1. at node I
    for (int i = 0; i < 5; ++i) {
        for (int k = 0; k < 6; ++k) {
            // 17686
            FEFC(k, i) = -FEF(k, i);

            // 2. at node J
            // 17690
            FEFC(k + 6, i) = FEF(k + 6, i);
        }
    }

    // 17699: form the transformation relating global
    // displacements and member stress resultants at node I and J

    // stress-displacement transformation relating
    // the 12 global components of displacement to
    // the 6 local components of member loads
    // located at node I and at node J.
    SA.resize(18, 12);
    SA.Zero();
    for (int k1 = 1; k1 <= 10; k1 += 3) {
        double fac = -1.0;
        if (k1 > 4) {
            fac = 1.0;
        }
        int nrs = k1 - 1;
        for (int k2 = 1; k2 <= 10; k2 += 3) {
            int ncs = k2 - 1;
            for (int ir = 1; ir <= 3; ++ir) {
                int nr = nrs + ir;
                for (int ic = 1; ic <= 3; ++ic) {
                    int nc = ncs + ic;
                    SA(nr - 1, nc - 1) = 0.0;
                    for (int IN = 1; IN <= 3; ++IN) {
                        int n = ncs + IN;
                        SA(nr - 1, nc - 1) += fac * K(nr - 1, n - 1) *
                                              currT(ic - 1, IN - 1);
                    }
                }
            }
        }
    }

    return 0;
}

int Pipe::localGlobal() {
    Matrix Q(3, 3), QQ(3, 3);

    // 9643
    for (int IR = 1; IR <= 10; IR += 3) {
        int IRS = IR - 1;
        for (int IC = IR; IC <= 10; IC += 3) {
            int ICS = IC - 1;

            // 9648
            for (int I = 1; I <= 3; ++I) {
                int II = IRS + I;
                for (int J = 1; J <= 3; ++J) {
                    int JJ = ICS + J;

                    // 9652
                    Q(I - 1, J - 1) = K(II - 1, JJ - 1);
                }
            }

            // 9655
            for (int I = 1; I <= 3; ++I) {
                for (int J = 1; J <= 3; ++J) {
                    QQ(I - 1, J - 1) = 0.0;
                    for (int KN = 1; KN <= 3; ++KN) {
                        // 9659
                        QQ(I - 1, J - 1) +=
                            Q(I - 1, KN - 1) * currT(J - 1, KN - 1);
                    }
                }
            }

            // 9663
            for (int I = 1; I <= 3; ++I) {
                int II = IRS + I;
                for (int J = 1; J <= 3; ++J) {
                    int JJ = ICS + J;
                    K(II - 1, JJ - 1) = 0.0;
                    for (int KN = 1; KN <= 3; ++KN) {
                        // 9659
                        K(II - 1, JJ - 1) +=
                            currT(I - 1, KN - 1) * QQ(KN - 1, J - 1);
                    }
                }
            }
        }
    }

    // 9676
    for (int I = 0; I < 12; ++I) {
        for (int J = 0; J < 12; ++J) {
            K(J, I) = K(I, J);
        }
    }
    return 0;
}

int Pipe::loading(Vector &FAC) {
    // 9682
    // 1 - x gravity
    // 2 - y gravity
    // 3 - z gravity
    // 4 - thermal
    // 5 - internal pressure

    // form the participation factors from the five types of loading
    

    return 0;
}