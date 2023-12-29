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
#ifndef Pipe_h
#define Pipe_h

#include <CrdTransf.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <elementAPI.h>

#include <vector>

// straight pipe element
class Pipe : public Element {
    // Line: 9369
   private:
    std::vector<Node *> theNodes;
    ID connectedExternalNodes;
    CrdTransf *theCoordTransf;
    UniaxialMaterial *theMat;
    SectionForceDeformation *theSect;

    Matrix K;
    Vector P;

    double T0;        // stress free temperature
    double pressure;  // internal pressure

   public:
    Pipe()
        : Element(0, ELE_TAG_Pipe),
          theNodes(2, 0),
          connectedExternalNodes(2),
          theCoordTransf(0),
          theMat(0),
          theSect(0),
          K(),
          P(),
          T0(0.0),
          pressure(0.0) {}
    Pipe(int tag, int nd1, int nd2, CrdTransf &coordTransf,
         UniaxialMaterial &mat, SectionForceDeformation &sect,
         double t0 = 0, double pre = 0)
        : Element(tag, ELE_TAG_Pipe),
          theNodes(2),
          connectedExternalNodes(2),
          theCoordTransf(0),
          theMat(0),
          theSect(0),
          K(),
          P(),
          T0(t0),
          pressure(pre) {
        // nodes
        connectedExternalNodes(0) = nd1;
        connectedExternalNodes(1) = nd2;

        // transf tag
        int ndm = OPS_GetNDM();
        if (ndm == 2) {
            theCoordTransf = coordTransf.getCopy2d();
        } else if (ndm == 3) {
            theCoordTransf = coordTransf.getCopy3d();
        }

        if (!theCoordTransf) {
            opserr << "Pipe element -- failed to get "
                      "copy of coordinate transformation"
                   << coordTransf.getTag() << "\n";
            exit(-1);
        }

        // section
        theSect = sect.getCopy();
        if (theSect == 0) {
            opserr << "Pipe element - failed to "
                      "get a copy of section "
                   << sect.getTag() << "\n";
            exit(-1);
        }

        // material
        theMat = mat.getCopy();
        if (theMat == 0) {
            opserr << "Pipe element - failed to get a copy of "
                      "material with tag "
                   << mat.getTag() << "\n";
            exit(-1);
        }
    }

    ~Pipe() {
        if (theCoordTransf != 0) {
            delete theCoordTransf;
        }
        if (theMat != 0) {
            delete theMat;
        }
        if (theSect != 0) {
            delete theSect;
        }
    }

    const char *getClassType(void) const { return "Pipe"; };

    int getNumExternalNodes(void) const { return 2; }
    const ID &getExternalNodes(void) {
        return connectedExternalNodes;
    }
    Node **getNodePtrs(void) { return &theNodes[0]; }

    int getNumDOF(void) {
        int ndm = OPS_GetNDM();
        if (ndm == 2) {
            return 6;
        } else if (ndm == 3) {
            return 12;
        }
        return 0;
    }

    void setDomain(Domain *theDomain) {
        // check domain
        if (theDomain == 0) {
            opserr << "Pipe::setDomain -- Domain is null\n";
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

        // set domain
        this->DomainComponent::setDomain(theDomain);

        // set crdTransf
        if (theCoordTransf->initialize(theNodes[0], theNodes[1]) !=
            0) {
            opserr << "Pipe::setDomain  tag: " << this->getTag()
                   << " -- Error initializing coordinate "
                      "transformation\n";
            exit(-1);
        }
        double L = theCoordTransf->getInitialLength();
        if (L == 0.0) {
            opserr << "Pipe::setDomain  tag: " << this->getTag()
                   << " -- Element has zero length\n";
            exit(-1);
        }
    }

    int commitState(void) {
        int retVal = 0;
        // call element commitState to do any base class stuff
        if ((retVal = this->Element::commitState()) != 0) {
            opserr << "Pipe::commitState () - failed in "
                      "base class";
        }
        retVal += theCoordTransf->commitState();
        retVal += theMat->commitState();
        retVal += theSect->commitState();

        K.resize(getNumDOF(), getNumDOF());
        K.Zero();

        P.resize(getNumDOF());
        P.Zero();

        return retVal;
    }

    int revertToLastCommit(void) {
        int retVal = 0;
        retVal += theCoordTransf->revertToLastCommit();
        retVal += theMat->revertToLastCommit();
        retVal += theSect->revertToLastCommit();
        return retVal;
    }

    int revertToStart(void) {
        int retVal = 0;
        retVal += theCoordTransf->revertToStart();
        retVal += theMat->revertToStart();
        retVal += theSect->revertToStart();
        return retVal;
    }

    int update(void) {
        int retVal = 0;
        retVal += theCoordTransf->update();

        // maybe: material setTrialStrain
        // maybe: section setTrialSectionDeformation
        return retVal;
    }

    const Matrix &getTangentStiff(void) { return K; }
    const Matrix &getInitialStiff(void) { return K; }
    const Matrix &getMass(void) { return K; }

    void zeroLoad(void) { P.Zero(); }

    int addLoad(ElementalLoad *theLoad, double loadFactor) {
        return 0;
    }

    int addInertiaLoadToUnbalance(const Vector &accel) { return 0; }

    const Vector &getResistingForce(void) { return P; }
    const Vector &getDampingForce(void) { return P; }
    const Vector &getResistingForceIncInertia(void) { return P; }

    int sendSelf(int commitTag, Channel &theChannel) { return 0; }
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) {
        return 0;
    }

    void Print(OPS_Stream &s, int flag = 0) {}
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes = 0, int numModes = 0) {
        return 0;
    }

    Response *setResponse(const char **argv, int argc,
                          OPS_Stream &s) {
        return 0;
    }

    int getResponse(int responseID, Information &info) { return 0; }

    int setParameter(const char **argv, int argc, Parameter &param) {
        return -1;
    }
    int updateParameter(int parameterID, Information &info) {
        return -1;
    }
};

void *OPS_PipeElement() {
    // line: 9343
    //               TAG  TP  I   J   MAT SEC TO   P  VECXZ     INC
    // READ (5,1040) INEL,R1,INI,INJ,IMAT,ISP,TRI,PRI,X2,X3,X4,INC 
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "Invalid #args,  want: element pipe "
                  "tag? nd1? nd2? matTag? secTag? transfTag? "
                  "<To? p?>\n";
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
    double data[2] = {0, 0};
    numData = OPS_GetNumRemainingInputArgs();
    if (numData > 2) {
        numData = 2;
    }
    if (numData > 0) {
        if (OPS_GetDoubleInput(&numData, data) < 0) {
            opserr << "WARNING: invalid data for pipe element\n";
            return 0;
        }
    }

    // check objects
    auto *theTrans = OPS_getCrdTransf(iData[5]);
    if (theTrans == 0) {
        opserr << "WARNING: CrdTransf " << iData[5]
               << " is not found\n";
        return 0;
    }

    auto *theSect = OPS_getSectionForceDeformation(iData[4]);
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[4]
               << " is not found\n";
        return 0;
    }

    auto *theMat = OPS_getUniaxialMaterial(iData[3]);
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[3]
               << " is not found\n";
        return 0;
    }

    auto *ele = new Pipe(iData[0], iData[1], iData[2], *theTrans,
                         *theMat, *theSect, data[0], data[1]);

    return ele;
}

#endif
