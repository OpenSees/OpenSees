#include "PY_Macro2D.h"

#include <elementAPI.h>
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Renderer.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <ElementResponse.h>
#include <math.h>

static int nDOF = 4; // 2 dof at pile, 2 at free field	!!!!!!!!

// initialise the class wide variables
Matrix PY_Macro2D::theMatrix(nDOF,nDOF);
Vector PY_Macro2D::theVector(nDOF);

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numPY_Macro2D = 0;

OPS_Export void *
OPS_PY_Macro2D(void)
{
  if (numPY_Macro2D == 0) {
    opserr << "PY_Macro2D element - Written by V.Varun and A.Shafiee, Georgia Tech Copyright 2009\n";
    numPY_Macro2D++;
  }
  // get the id and end nodes
  int iData[4];		//!!!!!!!!!!
  double dData[14];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data for PY_Macro2D\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 14;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element data for PY_Macro2D element with tag: " << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  Element *theEle = new PY_Macro2D(eleTag, iData[1], iData[2], dData[0], dData[1], dData[2], dData[3], dData[4],
	  dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], iData[13]);


  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating PY_Macro2D element with tag " << eleTag << endln;
    return 0;
  }

  return theEle;
}


// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.
PY_Macro2D::PY_Macro2D(int tag,
		   int node1,
		   int node2,
		   double _K,
		   double _py,
		   double _a,
		   double _b,
		   double _g,
		   double _m1,
		   double _m2,
		   double _w1,
		   double _p1,
		   double _S1,
		   double _beta,
		   double _s1,
		   double _tolerance,
		   int _maxNumIter)
 :Element(tag, ELE_TAG_PY_MACRO2D),
  K(_K), py(_py), a(_a), b(_b), g(_g), m1(_m1), m2(_m2),
  w1(_w1), p1(_p1), S1(_S1), beta(_beta), s1(_s1), 
  tolerance(_tolerance), maxNumIter(_maxNumIter),
  Ttangent(0.0), Tforce(0.0), Tz(0.0), TU(0.0), TW(0.0), TS(1.0), TS0(1.0), CW(0.0), CS(1.0), CS0(1.0), 
  Ctangent(0.0), Cforce(0.0), Cz(0.0), CU(0.0), Tt(0.0), Ct(0.0), trans(1,4),
  connectedExternalNodes(2)
{
  connectedExternalNodes(0) = node1;
  connectedExternalNodes(1) = node2;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
PY_Macro2D::PY_Macro2D()
 :Element(0, ELE_TAG_PY_MACRO2D),
  K(0), py(0), a(0), b(0), g(0), m1(0), m2(0),
  w1(0), p1(0), S1(0), beta(0), s1(0),tolerance(0), maxNumIter(0),
  Ttangent(0.0), Tforce(0.0), Tz(0.0), TU(0.0), TW(0.0), TS(1.0), TS0(1.0), CW(0.0), CS(1.0), CS0(1.0),
  Ctangent(0.0), Cforce(0.0), Cz(0.0), CU(0.0), Tt(0.0), Ct(0.0), trans(1,4),
  connectedExternalNodes(2)
{
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
PY_Macro2D::~PY_Macro2D()
{
  // Does nothing
}


int
PY_Macro2D::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
PY_Macro2D::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
PY_Macro2D::getNodePtrs(void)
{
  return theNodes;
}

int
PY_Macro2D::getNumDOF(void)
{
  return nDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the truss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
PY_Macro2D::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }

  // first set the node pointers
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  // if can't find both - send a warning message
  if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
    if (theNodes[0] == 0)
      opserr <<"PY_Macro2D::setDomain() - truss" << this->getTag() << " node " << Nd1 <<
	"does not exist in the model\n";
    else
      opserr <<"PY_Macro2D::setDomain() - truss" << this->getTag() << " node " << Nd2 <<
	"does not exist in the model\n";

    return;
  }

  // now determine the number of dof and the dimesnion
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();

  if (dofNd1 != 2) {		//!!!!!!!!!!!!!
    opserr <<"WARNING PY_Macro2D::setDomain(): node 1: " << Nd1 << " needs 3 dof\n ";
    return;
  }

  if (dofNd2 != 2) {
    opserr <<"WARNING PY_Macro2D::setDomain(): node 2: " << Nd2 << " needs 2 dof\n ";
    return;
  }

  // call the base class method
  this->DomainComponent::setDomain(theDomain);

   // now determine the length & transformation matrix
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);

    double L = sqrt(dx*dx + dy*dy);

    if (L == 0.0) {
      opserr << "WARNING PY_Macro2D::setDomain() - PY_Macro2D " << this->getTag() <<
	" has zero length\n";
      return;  // don't go any further - otherwise divide by 0 error
    }

    double cs = dx/L;
    double sn = dy/L;

    trans(0,0) = -cs;
    trans(0,1) = -sn;
    trans(0,2) = cs;
    trans(0,3) = sn;
}


int
PY_Macro2D::commitState()
{
  CU = TU;
  Cz = Tz;
  Cforce = Tforce;
  Ctangent = Ttangent;
  CW = TW;
  CS = TS;
  CS0 = TS0;
  Ct = Tt;	

  return 0;
}

int
PY_Macro2D::revertToLastCommit()
{
  TU = CU;
  Tz = Cz;
  Tforce = Cforce;
  Ttangent = Ctangent;

  TW= CW;
  TS = CS;
  TS0 = CS0;
//  Tt = Ct;

  return 0;
}

int
PY_Macro2D::revertToStart()
{

  TU = 0.0;
  Tz = 0.0;
  Tforce = 0.0;
  Ttangent = K;

  TW = 0.0;
  TS = 1.0;
  TS0 = 1.0;
  Tt = 0.0;

  return 0;
}

double
PY_Macro2D::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
}

int
PY_Macro2D::update(void)
{
  Domain *theDomain=this->getDomain();
  Tt = theDomain->getCurrentTime();
  double dt = Tt - Ct;

  // determine the strain
  const Vector &disp1 = theNodes[0]->getTrialDisp();
  const Vector &disp2 = theNodes[1]->getTrialDisp();

  double Ru = disp1(1);          // pore water pressure

  // Use displacements to find the total strain in the element
  TU = 0.0;
  for (int i=0; i<2; i++)
  {
    TU -= (disp2(i)-disp1(i)) * trans(0,i);
  }
  // Find the change in strain
  double dU = TU - CU;

  // Declare the other variables required
  double dz, f, f_, Tzold, Tznew;

  dz = K/py*(1-(tanh(a*abs(Cz))/tanh(a))*(b+g*signum(dU*Cz)))*dU;
  Tz = Cz+dz;


 // Pore Pressure Generation Model
	Tforce = py*Tz*CS;
	Ttangent = K*(1-(tanh(a*fabs(Tz))/tanh(a))*(b+g*signum(dU*Tz)))*TS;

	TW = CW;
	double dSb = 0.0;
	if (fabs(Tz) <= 0.67*m2/m1)
	{
		TW = CW+fabs(Tforce*dU)/py/(py/K);
		dSb = exp(-1*pow(TW/w1,1.4))*1.4*pow(TW/w1,0.4)*fabs(Tforce*dU)/py/(py/K)/w1;
	}

	double Sff = 1-Ru;
	double dSd = beta/(0.01+0.99*fabs(Sff-CS0))*pow(CS,p1) *dt/(1+beta/(0.01+0.99*fabs(Sff-CS0))*pow(CS,p1) *dt)*(Sff-CS);
	TS0 = CS0 - dSb + dSd;

	if (fabs(Tz) <= 0.67*m2/m1)
	{
		TS = TS0;
	} else
	{
		double alp = 0.67*m2/m1;
		TS = TS0*(1+alp*alp)/(fabs(Tz)*alp+pow((Tz*alp)*(Tz*alp)+(1-Tz*Tz)*(1+alp*alp),0.5));
	}

 // Compute force and tangent
//	Tforce = py*Tz*TS;
//	Ttangent = K*(1-(tanh(a*fabs(Tz))/tanh(a))*(b+g*signum(dU*Tz)))*TS;



  return 0;
}


const Matrix &
PY_Macro2D::getTangentStiff(void)
{
  theMatrix.Zero();
  theMatrix = trans^trans;
  theMatrix *= Ttangent;

  return theMatrix;
}


const Matrix &
PY_Macro2D::getInitialStiff(void)
{
  theMatrix.Zero();

  theMatrix = trans^trans;
  theMatrix *= K;

  return theMatrix;
}

const Matrix &
PY_Macro2D::getDamp(void)
{
  theMatrix.Zero();
  return theMatrix;
}


const Matrix &
PY_Macro2D::getMass(void)
{
  theMatrix.Zero();
  return theMatrix;
}

const Vector &
PY_Macro2D::getResistingForce()
{
  theVector.Zero();
  for (int i=0; i<4 ; i++)
	  theVector(i) = trans(0,i)*Tforce;

  return theVector;
}


const Vector &
PY_Macro2D::getResistingForceIncInertia()
{
  return this->getResistingForce();
}

int
PY_Macro2D::sendSelf(int commitTag, Channel &theChannel)
{
  return 1;
}

int
PY_Macro2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
PY_Macro2D::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag();
  s << " type: PY_Macro2D  iNode: " << connectedExternalNodes(0);
  s << " jNode: " << connectedExternalNodes(1) << endln;
}


Response*
PY_Macro2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","PY_Macro2D");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types for the PY_Macro2D
  //


  if (strcmp(argv[0],"S")) {
    output.tag("ResponseType", "S1");
    theResponse = new ElementResponse(this, 3, 0.0);
    output.endTag();
  }

  return theResponse;
}

 int
   PY_Macro2D::getResponse(int responseID, Information &eleInfo)
 {
   double strain;

   switch (responseID) {
   case 1:
     return eleInfo.setDouble(S1);

   default:
     return 0;
   }
 }

int 
PY_Macro2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}
