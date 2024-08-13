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
class Channel;
class Information;
class CrdTransf;
class Response;
class Renderer;
class SectionForceDeformation;

// straight pipe element
class Pipe : public Element {
   private:
    PipeMaterial *theMat;
    PipeSection *theSect;

    double T0;          // stress free temperature
    double pressure;    // internal pressure
    bool thermalLoad;   // include thermal load
    bool pressureLoad;  // include internal pressure load

    Matrix K;
    Vector P;
    Vector Q;

    Matrix kb;
    Vector q;
    double q0[5];  // Fixed end forces in basic system (no torsion)
    double p0[5];  // Reactions in basic system (no torsion)

    double wx;
    double wy;
    double wz;

    Node *theNodes[2];
    ID connectedExternalNodes;
    CrdTransf *theCoordTransf;

   public:
    Pipe();
    Pipe(int tag, int Nd1, int Nd2, PipeMaterial &mat,
         PipeSection &sect, double to = 0.0, double pre = 0.0,
         bool tl = true, bool pl = true);

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
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes = 0, int numModes = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

   public:
    double aveTemp();
    int updateSectionData(double &A, double &Jx, double &Iy,
                          double &Iz, double &rho, double &alphaV);
    int updateMaterialData(double &E, double &nu, double &G,
                           double &alp);
    static int crossProduct(const Vector &A, const Vector &B,
                            Vector &res);
    static int nextTransfTag();
    void shearCoefficients(double &B1, double &B2, double &C1,
                           double &C2);
    Matrix &refK() { return K; }
    Vector &refQ() { return Q; }
    Vector &refP() { return P; }
    double &refwx() { return wx; }
    double &refwy() { return wy; }
    double &refwz() { return wz; }
    PipeMaterial *getMaterial() { return theMat; }
    PipeSection *getSection() { return theSect; }
    double getT0() { return T0; }
    double getPressure() { return pressure; }
    bool inclThermalLoad() { return thermalLoad; }
    bool inclPressureLoad() { return pressureLoad; }
    void initLoad();
    void getSectionForce(double xL, Vector &s);
};

#endif
