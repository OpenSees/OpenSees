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
                  "<To? p?>\n";
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

    auto *ele = new Pipe(iData[0], iData[1], iData[2], *theMat,
                         *theSect, data[0], data[1]);

    return ele;
}

Pipe::Pipe()
    : Element(0, ELE_TAG_Pipe),
      theNodes(2, 0),
      connectedExternalNodes(2),
      theMat(0),
      theSect(0),
      K(),
      P(),
      localY(3),
      T0(0.0),
      pressure(0.0),
      currLN(0.0),
      currTavg(0.0),
      currT() {}

Pipe::Pipe(int tag, int nd1, int nd2, PipeMaterial &mat,
           PipeSection &sect, double vecyx, double vecyy,
           double vecyz, double t0, double pre)
    : Element(tag, ELE_TAG_Pipe),
      theNodes(2),
      connectedExternalNodes(2),
      theMat(0),
      theSect(0),
      K(),
      P(),
      localY(3),
      T0(t0),
      pressure(pre),
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

    // compute stiffness

    return retVal;
}

const Matrix &Pipe::getTangentStiff(void) { return K; }
const Matrix &Pipe::getInitialStiff(void) { return K; }
const Matrix &Pipe::getMass(void) { return K; }

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
int Pipe::tangks(const PipeMaterialTemperaturePoint &Tpt) {
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
}
