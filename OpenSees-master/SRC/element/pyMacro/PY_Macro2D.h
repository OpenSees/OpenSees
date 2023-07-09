#ifndef PY_Macro2D_h
#define PY_Macro2D_h

// Written:
//
// Description:
//

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class PY_Macro2D : public Element
{
  public:
    PY_Macro2D(int tag,
	       int Nd1,
	       int Nd2,
	       double K,
	       double py,
	       double a,
	       double b,
	       double g,
	       double m1,
	       double m2,
	       double w1,
	       double p1,
	       double S1,
	       double beta,
	       double s1,
	       double tolerance,
	       int maxNumIter);

    PY_Macro2D();
    ~PY_Macro2D();

    const char *getClassType(void) const {return "PY_Macro2D";};

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information
 //   const Matrix &getKi(void);		//!!!!!!!!!!!
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);
    double signum(double);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:

  private:
    double K;
    double py;
    double a;
    double b;
    double g;
    double m1;
    double m2;
    double w1;
    double p1;
    double S1;
    double beta;
    double s1;
    double tolerance;
    int maxNumIter;
    

    double Ttangent;
    double Tforce;
    double Tz; 
    double TU; 
    double TW;
    double TS;
    double TS0;
    double CW;
    double CS;
    double CS0;

    double Ctangent;
    double Cforce;
    double Cz; 
    double CU; 

    double Tt;
    double Ct;
    
    Matrix trans; //!!!!!!!!!!!!

    ID  connectedExternalNodes;     // contains the tags of the end nodes

    Node *theNodes[2];

    // static data - single copy for all objects of the class
    static Matrix theMatrix;   // class wide matrix for 5*5
    static Vector theVector;   // class wide Vector for size 5
};

#endif




