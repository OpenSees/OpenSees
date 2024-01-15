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

#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <PipeMaterial.h>
#include <PipeSection.h>
#include <Vector.h>
#include <elementAPI.h>

#include <vector>

class Node;

// straight pipe element
class Pipe : public Element {
    // Line: 9369
   private:
    std::vector<Node *> theNodes;
    ID connectedExternalNodes;
    PipeMaterial *theMat;
    PipeSection *theSect;

    Matrix M;
    Matrix K;
    Vector P;
    Vector localY;

    double T0;        // stress free temperature
    double pressure;  // internal pressure

    double currLN;    // line length
    double currTavg;  // average temperature
    Matrix currT;     // transformation matrix

   public:
    Pipe();
    Pipe(int tag, int nd1, int nd2, PipeMaterial &mat,
         PipeSection &sect, double vecyx = 0, double vecyy = 0,
         double vecyz = 0, double t0 = 0, double pre = 0);

    ~Pipe();

    const char *getClassType(void) const;

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);

    void setDomain(Domain *theDomain);

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
    const Vector &getResistingForceIncInertia(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &info);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // utility
   private:
    // default convention for local y-axis
    void defaultLocalY(const Vector &localX, Vector &deLocalY);

    // computation of direction cosine array for the local axes of
    // pipe tangent element
    int tangdc(Matrix &dircos);

    // computation of the element stiffness and load matrices
    // for a tangent (straight) pipe element
    int tangks(const PipeMaterialTemperaturePoint &Tpt);
};

#endif
