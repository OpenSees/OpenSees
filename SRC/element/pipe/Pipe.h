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
   private:
    std::vector<Node *> theNodes;
    ID connectedExternalNodes;
    CrdTransf *theCoordTransf;
    UniaxialMaterial *theMat;
    SectionForceDeformation *theSect;

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

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    int setDamping(Domain *theDomain, Damping *theDamping);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getDampingForce(void);
    const Vector &getResistingForceIncInertia(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes = 0, int numModes = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

   private:
    double A, E, G, Jx, Iy, Iz;

    double rho;
    int cMass;

    int releasez;  // moment release for bending about z-axis 0=none,
                   // 1=I, 2=J, 3=I,J
    int releasey;  // same for y-axis

    static Matrix K;
    static Vector P;
    Vector Q;

    static Matrix kb;
    Vector q;
    double q0[5];  // Fixed end forces in basic system (no torsion)
    double p0[5];  // Reactions in basic system (no torsion)

    double wx;
    double wy;
    double wz;

    Node *theNodes[2];

    ID connectedExternalNodes;

    CrdTransf *theCoordTransf;

    Damping *theDamping;
};

void *OPS_PipeElement() {
    // 7722 - 7747
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 7) {
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
