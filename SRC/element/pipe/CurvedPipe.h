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
#ifndef CurvedPipe_h
#define CurvedPipe_h

#include <Element.h>
#include <Matrix.h>
#include <Pipe.h>
#include <PipeMaterial.h>
#include <PipeSection.h>
#include <Vector.h>

#include <vector>

class Node;

// CurvedPipe element
class CurvedPipe : public Element {
   private:
    Pipe *pipeEle;

    Vector center;
    double radius;
    double theta0;
    double tolWall;
    double kp;          // flexibility factor
    double Length;      // length of chord
    bool intersection;  // if center is an intersection point
    Matrix alg;         // local-global transformation
    Matrix abl;         // basic-local transformation

    static std::vector<double> gaussPts;

   public:
    CurvedPipe();
    CurvedPipe(int tag, int Nd1, int Nd2, PipeMaterial &mat,
               PipeSection &sect, const Vector &c, double to = 0.0,
               double pre = 0.0, double tol = 0.1,
               bool inter = false,
               bool tl = true, bool pl = true);

    ~CurvedPipe();

    const char *getClassType() const;

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    int update();
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();

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
    int getTheta0();

    void plw(Vector &vec);
    void bx(double theta, Matrix &mat);
    void Spx(double theta, Vector &vec);
    void fs(double theta, Matrix &mat);
    void fb(double theta, Matrix &mat);
    void ubno(double theta, Vector &);
    int kb(Matrix &mat, Vector &vec);
    void integrateGauss(double a, double b, Matrix &res,
                        Vector &resv);
    int getSectionForce(double theta, Vector &s);
};

#endif
