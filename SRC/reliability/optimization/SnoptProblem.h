/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **                                 **
** ****************************************************************** */
                                                                        

//
// Written by Quan Gu (qgu@ucsd.edu)
//
#ifndef SNOPTPROBLEM_H
#define SNOPTPROBLEM_H

#include "snopt.h"
#include <ReliabilityDomain.h>
#include <ProbabilityTransformation.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <Matrix.h>
#include <Vector.h>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <fstream>
using std::ofstream;

// class SnoptProblemoblem performs problem set-up, initialization,
// and problem-specific variable storage.
class SnoptProblem: public FindDesignPointAlgorithm {

protected:
  //************************************************************
  // The variables below must be set by the user with initialize
  char       Prob[200];

  integer    n, neF;
  integer    ObjRow;
  doublereal ObjAdd;

  integer    lenA, neA;
  integer    *iAfun, *jAvar;
  doublereal *A;
  integer    lenG, neG;
  integer    *iGfun, *jGvar;

  doublereal *x, *xlow, *xupp, *xmul;
  doublereal *F, *Flow, *Fupp, *Fmul;

  integer *xstate, *Fstate;

  char *xnames, *Fnames;
  integer nxnames, nFnames;

  My_fp usrfun;
  //***********************************************************
  //***********************************************************
  void userDataSet();
  void errMsgExit( char *var );

  integer    inform, fortranStyleObj, fortranStyleAG;
  integer    initCalled;
  integer    minrw, miniw, mincw;
  integer    lenrw, leniw, lencw;
  doublereal *rw;
  integer    *iw;
  char       *cw;

  integer    iSpecs, iSumm, iPrint;
  char       specname[200], printname[200];
  integer    spec_len, prnt_len;

  void init2zero();
  void setMemory();
  void alloc    ( integer lencw, integer leniw, integer lenrw );
  void realloc  ( integer lencw, integer leniw, integer lenrw );
  void memcpyIn ( char *tcw, integer *tiw, doublereal *trw,
		  integer tlencw, integer tleniw, integer tlenrw);
  void memcpyOut( char *tcw, integer *tiw, doublereal *trw,
		  integer tlencw, integer tleniw, integer tlenrw);
public:
   SnoptProblem(	int passedMaxNumberOfIterations, 
					GFunEvaluator *passedGFunEvaluator,
					GradGEvaluator *passedGradGEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					int printFlag,
					char *fileNamePrint,
					Vector *startPoint, char * probType,
					ReliabilityDomain * passedReliabilityDomain);

   SnoptProblem(){this->outputFile =0;};


  ~SnoptProblem();
  void increment();
  void decrement();
  void computeJac();

  int snmema( integer &mincw, integer &miniw, integer &minrw);

  void init();
  void setParameter    ( char *stroptin );
  void getParameter    ( char *stroptin, char *stroptout );
  void setIntParameter ( char *stropt,   integer     opt );
  void getIntParameter ( char *stropt,   integer    &opt );
  void setRealParameter( char *stropt,   doublereal  opt );
  void getRealParameter( char *stropt,   doublereal &opt );
  void solve           ( integer starttype );
  void setPrintFile    ( char printname[] );
  void setSpecFile     ( char specname[] );

  // Functions that set up the problem data:
  void setProblemSize( integer n, integer neF );
  void setObjective  ( integer ObjRow, doublereal ObjAdd );
  void setA          ( integer lenA, integer *iAfun, integer *jAvar,
		       doublereal *A );
  void setNeA        ( integer neA );
  void setG          ( integer lenG, integer *iGfun, integer *jGvar );
  void setNeG        ( integer neG );
  void setX          ( doublereal *x, doublereal *xlow, doublereal *xupp,
		       doublereal *xmul, integer *xstate );
  void setF          ( doublereal *F, doublereal *Flow, doublereal *Fupp,
		       doublereal *Fmul, integer *Fstate );
  void setXNames     ( char *xnames, integer nxnames );
  void setFNames     ( char *Fnames, integer nFnames );
  void setProbName   ( char *Prob );
  void setUserFun    ( My_fp usrfun );

// ---------------- for reliability part only --------------------------
// ----------------------------------------------------------------  
public:
	int findDesignPoint(ReliabilityDomain *theReliabilityDomain);

	const Vector & get_x();
	const Vector & get_u();
	const Vector & get_alpha();
	const Vector & get_gamma();
	int getNumberOfSteps();
	const Vector & getSecondLast_u();
	const Vector & getSecondLast_alpha();
	const Vector & getLastSearchDirection();
	double getFirstGFunValue();
	double getLastGFunValue();
	const Vector & getGradientInStandardNormalSpace();
	int    getNumberOfEvaluations();
    double * getHessian();
	// Quan and Michele Jan 2006  
   int setStartPt(Vector *);
protected:
public:  // not good, but avoid large amount of functions.
	int  setUSecondLast(const Vector &);
	int  setU(const Vector &);
	int  setAlphaSecondLast(const Vector &);

//private:	
	
	Vector reliability_x;
	Vector u;
	

	// The reliability domain and tools for the analysis
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	GradGEvaluator *theGradGEvaluator;
	ProbabilityTransformation *theProbabilityTransformation;

	// Data members set when the object is created
	int maxNumberOfIterations;

	// Data members where the results are to be stored

	Vector alpha;
	Vector gradientInStandardNormalSpace;
	Vector gamma;
	Vector uSecondLast;
	Vector alphaSecondLast;
	int ii;
//	Vector searchDirection;
	double Gfirst;
	double Glast;

	// Data members set through the call when a job is to be done
	Vector *startPoint;
	Vector *designPoint_uStar;

	int printFlag;
	char *fileNamePrint;
	int numberOfEvaluations;

};

#endif
