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

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "snopt.h"

#include <SnoptProblem.h>
#include <snoptfilewrapper.h>
#include <ToyFunction.h>
using namespace std;

//------------------------------------
#include <FindDesignPointAlgorithm.h>
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <MatrixOperations.h>
#include <Matrix.h>
#include <Vector.h>
#include <GammaRV.h>

#include <fstream>
#include <iomanip>
#include <iostream>

using std::ifstream;
using std::ios;

using std::setw;
using std::setprecision;


SnoptProblem::SnoptProblem(int passedMaxNumberOfIterations, 
					GFunEvaluator *passedGFunEvaluator,
					GradGEvaluator *passedGradGEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
			   bool passedStartAtOrigin,
					int pPrintFlag,
					char *pFileNamePrint,
					char * probType, 
					ReliabilityDomain * passedReliabilityDomain)
:FindDesignPointAlgorithm(passedReliabilityDomain), 
  iSpecs(0), iSumm(6), iPrint(0), initCalled(0)
{



//-------------------- step 1 SNOPT init() -------------------------------------
  init2zero();

  // Nothing incremented yet
  fortranStyleObj = 0;
  fortranStyleAG  = 0;

  //  iSpecs =  4;
  //  iSumm  = 6;
  //  iPrint = 15;

  // Create temporary memory for the call to sninit_.
  // Lengths must all be >= 500.
  lencw  = 500;
  leniw  = 500;
  lenrw  = 500;
  this->alloc( 500, 500, 500 );

  // sninit_ "undefines" the optional parameters

  this->init();
// ---------------------- step 2 save reliabiity pointer ------------------------------------
// if probType =="reliability" then
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	theGFunEvaluator				= passedGFunEvaluator;
	theGradGEvaluator				= passedGradGEvaluator;
	theProbabilityTransformation	= passedProbabilityTransformation;
	startAtOrigin = passedStartAtOrigin;
	printFlag						= pPrintFlag;
	numberOfEvaluations =0;
	fileNamePrint = new char[256];
	if (printFlag != 0) {
		strcpy(fileNamePrint,pFileNamePrint);
	}
	else {
		strcpy(fileNamePrint,"SNOPSearchpoints.out");
	}




// ---------------------- step 3 snopt construction for reliability -------------------------------
//	this->theReliabilityDomain = passedReliabilityDomain;

 
//--------------------- step 3.1 for snopt memory alloc -------------------------------------
  // Allocate and initialize;
  integer numOfMaxIter = this->maxNumberOfIterations;  // this number need to be input
  this->n     =  theReliabilityDomain->getNumberOfRandomVariables();   // n need to be got from numOfRandomVariable x1 x2 x3 ..
  this->neF   =  2;   // in reliability proble, is 2   Function number is 2
  this->lenA  = 1;   // not needed

  this->iAfun = new integer[lenA];
  this->jAvar = new integer[lenA];
  this->A  = new doublereal[lenA];

  this->lenG   = n*2; //// for reliability only
  this->iGfun = new integer[lenG];
  this->jGvar = new integer[lenG];

  this->x      = new doublereal[n];
  this->xlow   = new doublereal[n];
  this->xupp   = new doublereal[n];
  this->xmul   = new doublereal[n];
  this->xstate = new    integer[n];

  this->F      = new doublereal[neF];
  this->Flow   = new doublereal[neF];
  this->Fupp   = new doublereal[neF];
  this->Fmul   = new doublereal[neF];
  this->Fstate = new integer[neF];

  this-> nxnames = 1;
  this-> nFnames = 1;
  this->xnames = new char[nxnames*8];
  this->Fnames = new char[nFnames*8];

//--------------------- step 3.2 set problem bound for reliability --------------------------------------

  this-> ObjRow = 0;
  this-> ObjAdd = 0;

  // Set the upper and lower bounds.
  int i;
  for (i=0;i<n;i++){
	xlow[i]   = -1e20;
	xupp[i]   = 1e20;
	xstate[i] =    0;
  }
 
  
  Flow[0] = 0.0;   Flow[1] = 0.0;
  Fupp[0] =  1e20; Fupp[1] = 0.0;


//--------------------- step 3.3 init start point to 0 note reliability_x has no memory --------------------
  
  //set x, x[i] is in fact u in stand normal space, while Vector reliability_x is true x in x-space !!
  // this->u should copy value from x from  time to time to make sure u=x


  for (int ii=0;ii<n;ii++) { x[ii]=0.0;}
  this->u.Zero();



//--------------------- step 3.4 set G sequence -----------------------------------------
  for(i=0;i<neF;i++) Fmul[i]   = 0;
   
// should be: (0,0) (0,1) (0,2) (0,3)..  (1,0) (1,1) (1,2) (1,3)..
  this->neG = 0;
  for (i=0; i<neF;i++){
	  for(int j =0; j<n; j++){
		 iGfun[neG] = i;
		 jGvar[neG] = j;
		 neG++;
		}
  }








  this->neA = 0;

//--------------------- step 3.5 set others  -----------------------------------------

  this->usrfun= toyusrfg_ ;  // Sets the usrfun that supplies G and F.

    this->setPrintFile   ( "Toy1.out" );
    this->setSpecFile    ( "sntoya.spc" );
    this->setIntParameter( "Derivative option", 1 );
    this->setIntParameter( "Major Iteration limit", maxNumberOfIterations );
    this->setProbName   ( probType );
	this->ii=0;

	// ---

	this->outputFile =0;

}

SnoptProblem::~SnoptProblem()
{
  //Close print and spec files if necessary.
  if (iPrint != 0 ) {
    snoptclose_( &iPrint );
  }
  if (iSpecs != 0 ) {
    snoptclose_( &iSpecs );
  }

  //Delete work arrays.
  delete [] rw;
  delete [] iw;
  delete [] cw;

  // delete []iAfun;  delete []jAvar;  delete []A;
  delete []iGfun;  delete []jGvar;

  delete []x;      delete []xlow;   delete []xupp;
  delete []xmul;   delete []xstate;

  delete []F;      delete []Flow;   delete []Fupp;
  delete []Fmul;   delete []Fstate;

  delete []xnames; delete []Fnames;

  if (this->outputFile !=0) delete outputFile; 
}

void SnoptProblem::init2zero()
{
  // Data that must be set by user.

  n = 0; neF = 0;

  ObjRow =  0;
  ObjAdd =  0;

  iAfun = 0; jAvar = 0;  A = 0;
  iGfun = 0; jGvar = 0;

  x = 0; xlow = 0; xupp = 0; xmul = 0;
  F = 0; Flow = 0; Fupp = 0; Fmul = 0;

  xstate  = 0; Fstate  = 0;
  nxnames = 0; nFnames = 0;

  usrfun = 0;

  neA  = -1;  // Indicate that neA is yet to be assigned
  neG  = -1;
  lenA = -1;
  lenG = -1;
}

void SnoptProblem::userDataSet()
{
  if ( n    == 0)  errMsgExit( "n"  );
  if ( neF  == 0)  errMsgExit( "neF");

  if ( x    == 0 ) errMsgExit( "x"    );
  if ( xlow == 0 ) errMsgExit( "xlow" );
  if ( xupp == 0 ) errMsgExit( "xupp" );
  if ( xmul == 0 ) errMsgExit( "xmul" );

  if ( F    == 0 ) errMsgExit( "F"    );
  if ( Flow == 0 ) errMsgExit( "Flow" );
  if ( Fupp == 0 ) errMsgExit( "Fupp" );
  if ( Fmul == 0 ) errMsgExit( "Fmul" );

  if ( xnames  == 0 ) errMsgExit( "xnames"  );
  if ( Fnames  == 0 ) errMsgExit( "Fnames"  );
  if ( nxnames == 0 ) errMsgExit( "nxnames" );
  if ( nFnames == 0 ) errMsgExit( "nFnames" );

  if ( usrfun ==  0 ) errMsgExit( "usrfun" );
  if ( lenA   == -1 ) errMsgExit( "lenA" );
  if ( lenG   == -1 ) errMsgExit( "lenG" );

  if (( neA > 0) & (iAfun == 0 )) errMsgExit( "iAfun" );
  if (( neA > 0) & (jAvar == 0 )) errMsgExit( "jAvar" );
  if (( neA > 0) & (A     == 0 )) errMsgExit( "A"     );

  if (( neG > 0) & (iGfun == 0 )) errMsgExit( "iGfun" );
  if (( neG > 0) & (jGvar == 0 )) errMsgExit( "jGvar" );
}

void SnoptProblem::errMsgExit( char *var )
{
  cerr << "****************************************************\n";
  cerr << "Error: " << var << " must be set prior to call to " << endl
       << "SnoptProblem::solve() or SnoptProblem::computeJac()!\n";
  cerr << "****************************************************\n";
  exit(1);
}

void SnoptProblem::setMemory()
{
  int memoryGuess;
  memoryGuess = this->snmema(mincw, miniw, minrw);
  if ( mincw > lencw | miniw > leniw | minrw > lenrw ) {
    // Reallocate memory while retaining the values set in sninit_
    this->realloc( mincw, miniw, minrw );
    // Save the lengths of the new work arrays.
    this->setIntParameter("Total real workspace   ", lenrw );
    this->setIntParameter("Total integer workspace", leniw );

    // Did we have to guess values of neA and neG for snmema()
    if ( memoryGuess == 1 ) {
      this->computeJac();
      memoryGuess = this->snmema(mincw, miniw, minrw);
      assert( memoryGuess == 0 );
      this->realloc( mincw, miniw, minrw );
      this->setIntParameter("Total real workspace   ", lenrw );
      this->setIntParameter("Total integer workspace", leniw );
    }
  }
}

void SnoptProblem::alloc( integer alencw, integer aleniw, integer alenrw )
{
  // Reset work array lengths.
  lencw = alencw;
  leniw = aleniw;
  lenrw = alenrw;

  // Allocate new memory for work arrays.
  cw = new char[8*lencw];
  iw = new integer[leniw];
  rw = new doublereal[lenrw];
}

void SnoptProblem::realloc(  integer alencw, integer aleniw, integer alenrw )
{
  // Call to this->alloc will overwrite these values => must save.
  integer tlencw = lencw;
  integer tleniw = leniw;
  integer tlenrw = lenrw;

  // Call to this->alloc will create new values for cw, iw, rw => must save.
  char       *tcw = cw;
  integer    *tiw = iw;
  doublereal *trw = rw;

  // Allocate new memory
  this->alloc   ( alencw, aleniw, alenrw );
  // Copy in old values, previously set
  this->memcpyIn( tcw, tiw, trw, tlencw, tleniw, tlenrw );

  // Delete temporary work arrays
  delete [] tcw;
  delete [] tiw;
  delete [] trw;
}

void SnoptProblem::memcpyIn( char *tcw, integer *tiw, doublereal *trw,
			     integer tlencw, integer tleniw, integer tlenrw )
{
  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy( cw, tcw, 8*mlencw*sizeof( char ) );
  memcpy( iw, tiw, mleniw*sizeof(integer));
  memcpy( rw, trw, mlenrw*sizeof(doublereal));
}

void SnoptProblem::memcpyOut( char *tcw, integer *tiw, doublereal *trw,
			     integer tlencw, integer tleniw, integer tlenrw )
{
  integer mlencw = lencw < tlencw ? lencw : tlencw;
  integer mleniw = leniw < tleniw ? leniw : tleniw;
  integer mlenrw = lenrw < tlenrw ? lenrw : tlenrw;

  memcpy( tcw, cw, 8*mlencw*sizeof( char ) );
  memcpy( tiw, iw, mleniw*sizeof(integer));
  memcpy( trw, rw, mlenrw*sizeof(doublereal));
}

void SnoptProblem::increment()
{
  if( !fortranStyleObj ) {
    //Increment row indicator.
    ObjRow++;
    fortranStyleObj = 1;
  }

  if( !fortranStyleAG ) {
    //Increment A indices.
    for( int k = 0; k < neA; k++ ) {
      iAfun[k]++; jAvar[k]++;
    }
    //Increment G indices.
    for(int  kk = 0; kk < neG; kk++ ) {
      iGfun[kk]++; jGvar[kk]++;
    }
    fortranStyleAG = 1;
  }
}

void SnoptProblem::decrement()
{
  if( fortranStyleObj ) {
    //Decrement row indicator.
    ObjRow--;
    fortranStyleObj = 0;
  }

  if (fortranStyleAG) {
    //Decrement A indices.
    for( int k = 0; k < neA; k++ ) {
      iAfun[k]--; jAvar[k]--;
    }
    //Decrement G indices.
    for(int  kk = 0; kk < neG; kk++ ) {
      iGfun[kk]--; jGvar[kk]--;
    }
    fortranStyleAG = 0;
  }
}

void SnoptProblem::computeJac()
{
  //Ensures all user data has been initialized.
  userDataSet();
  this->snmema( mincw, miniw, minrw );
  if ( mincw > lencw | miniw > leniw | minrw > lenrw ) {
    // Reallocate memory while retaining the values set in sninit_
    this->realloc( mincw, miniw, minrw );
    // Save the lengths of the new work arrays.
    this->setIntParameter("Total real workspace   ", lenrw );
    this->setIntParameter("Total integer workspace", leniw );
  }
  snjac_( &inform, &neF, &n, usrfun,
	  iAfun, jAvar, &lenA, &neA, A,
	  iGfun, jGvar, &lenG, &neG,
	  x, xlow, xupp, &mincw, &miniw, &minrw,
	  cw, &lencw, iw, &leniw, rw, &lenrw,
	  cw, &lencw, iw, &leniw, rw, &lenrw,
	  8*500, 8*500 );
  //snjac_ will generate fortran style arrays.
  fortranStyleAG = 1;
}

int SnoptProblem::snmema
   ( integer &amincw, integer &aminiw, integer &aminrw)
{
  int memoryGuess = 0;
  integer nxname = 1; integer nfname = 1;
  if ( neA < 0 ) {
    neA = n*neF;
    memoryGuess = 1;
  }
  if ( neG < 0 ) {
    neG = n*neF;
    memoryGuess = 1;
  }
  snmema_( &inform, &neF, &n, &nxname, &nfname, &neA, &neG,
	   &amincw, &aminiw, &aminrw, cw, &lencw, iw, &leniw,
	   rw, &lenrw, 8*500 );
  return memoryGuess;
}

void SnoptProblem::init()
{
  initCalled = 1;
  sninit_( &iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500 );
}

void SnoptProblem::setParameter( char *stropt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snset_( stropt, &iPrt, &iSum, &inform, cw, &lencw, iw, &leniw,
	  rw, &lenrw, stropt_len, 8*500 );
}

void SnoptProblem::getParameter( char *stroptin, char *stroptout )
{
  assert( initCalled == 1 );

  integer stroptin_len  = strlen(stroptin);
  integer stroptout_len = strlen(stroptout);
  sngetc_( stroptin, stroptout, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stroptin_len, stroptout_len, 8*500 );
}

void SnoptProblem::setIntParameter( char *stropt, integer opt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snseti_( stropt, &opt, &iPrt, &iSum, &inform,
	   cw, &lencw, iw, &leniw, rw, &lenrw, stropt_len, 8*500 );
}

void SnoptProblem::getIntParameter( char *stropt, integer &opt )
{
  assert( initCalled == 1 );
  integer stropt_len = strlen(stropt);
  sngeti_( stropt, &opt, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stropt_len, 8*500 );
}

void SnoptProblem::setRealParameter( char *stropt, doublereal opt )
{
  assert( initCalled == 1 );

  integer iPrt   = 0; // suppresses printing
  integer iSum   = 0;
  integer stropt_len = strlen(stropt);
  snsetr_( stropt, &opt, &iPrt, &iSum, &inform,
	   cw, &lencw, iw, &leniw, rw, &lenrw, stropt_len, 8*500 );
}

void SnoptProblem::getRealParameter( char *stropt, doublereal &opt )
{
  assert( initCalled == 1 );
  integer stropt_len = strlen(stropt);
  sngetr_( stropt, &opt, &inform, cw, &lencw, iw, &leniw,
	   rw, &lenrw, stropt_len, 8*500 );
}

void SnoptProblem::solve( integer starttype )
{
  assert( initCalled == 1 );
  //Ensures all user data initialized.
  userDataSet();
  //Unlike snjac_ we also need neA and neG to be set.
  if ( neA == -1 | neG == -1 ) {
    cerr << "Warning: neA and neG must be set before calling"
	 << "SnoptProblem::solve()\n";
    exit(1);
  }
  integer npname = strlen(Prob);
  integer nS, nInf;
  doublereal sInf;
  this->increment(); //Convert array entries to Fortran style
  this->setMemory();
  snopta_( &starttype, &neF, &n, &nxnames,
	   &nFnames,
	   &ObjAdd, &ObjRow, Prob,
	   usrfun, iAfun, jAvar, &lenA, &neA, A,
	   iGfun, jGvar, &lenG, &neG,
	   xlow, xupp, xnames, Flow,
	   Fupp, Fnames, x, xstate,
	   xmul, F, Fstate, Fmul,
	   &inform, &mincw, &miniw, &minrw, &nS, &nInf, &sInf,
	   cw, &lencw, iw, &leniw, rw, &lenrw,
	   cw, &lencw, iw, &leniw, rw, &lenrw,
	   npname, 8*(nxnames), 8*(nFnames), 8*500, 8*500);
  this->decrement();  //Convert array entries to C style
}

void SnoptProblem::setPrintFile(  char aprintname[] )
{
  assert( initCalled = 1 );
  if (iPrint != 0 ) {
    snoptclose_( &iPrint );
  }
  iPrint = 15;
  strcpy( printname, aprintname );  prnt_len = strlen(printname);
  snoptopenappend_( &iPrint, printname,   &inform, prnt_len );
  this->setIntParameter("Print file", iPrint);
}

void SnoptProblem::setSpecFile( char aspecname[] )
{
  assert( initCalled == 1 );
  if (iSpecs != 0 ) {
    snoptclose_( &iSpecs );
  }
  iSpecs = 4;
  strcpy( specname, aspecname );    spec_len = strlen(specname);
  snoptfilewrapper_( specname, &iSpecs, &inform, cw, &lencw,
		     iw, &leniw, rw, &lenrw, spec_len, 8*lencw);
  if( inform != 101 ){
    printf("Warning: unable to find specs file %s \n", specname);
  }
}

void SnoptProblem::setProblemSize( integer an, integer aneF )
{
  //  checkSet = checkSet+2;
  n   = an;
  neF = aneF;
}

void SnoptProblem::setObjective( integer aObjRow, doublereal aObjAdd )
{
  //checkSet = checkSet+2;
  ObjRow = aObjRow;
  ObjAdd = aObjAdd;
}

void SnoptProblem::setA( integer alenA, integer *aiAfun,
			 integer *ajAvar, doublereal *aA )
{
  //checkSet = checkSet+4;
  lenA  = alenA;
  iAfun = aiAfun;
  jAvar = ajAvar;
  A     = aA;
}

void SnoptProblem::setNeA( integer aneA )
{
  //checkSet = checkSet+1;
  neA = aneA;
}
void SnoptProblem::setG( integer alenG, integer *aiGfun,
			 integer *ajGvar)
{
  //checkSet = checkSet+3;
  lenG  = alenG;
  iGfun = aiGfun;
  jGvar = ajGvar;
}

void SnoptProblem::setNeG( integer aneG )
{
  //checkSet = checkSet+1;
  neG = aneG;
}

void SnoptProblem::setX( doublereal *ax, doublereal *axlow, doublereal *axupp,
			 doublereal *axmul, integer *axstate )
{
  //checkSet = checkSet+5;
  x      = ax;
  xlow   = axlow;
  xupp   = axupp;
  xmul   = axmul;
  xstate = axstate;
}

void SnoptProblem::setF( doublereal *aF, doublereal *aFlow, doublereal *aFupp,
			 doublereal *aFmul, integer *aFstate )
{
  //checkSet = checkSet+5;
  F      = aF;
  Flow   = aFlow;
  Fupp   = aFupp;
  Fmul   = aFmul;
  Fstate = aFstate;
}

void SnoptProblem::setXNames( char *axnames, integer anxnames )
{
  //checkSet = checkSet+2;
  xnames  = axnames;
  nxnames = anxnames;
}

void SnoptProblem::setFNames( char *aFnames, integer anFnames )
{
  //checkSet = checkSet+2;
  Fnames  = aFnames;
  nFnames = anFnames;
}

void SnoptProblem::setProbName( char *aProb )
{
  //checkSet = checkSet+1;
  sprintf(Prob, "%8s", aProb );
}

void SnoptProblem::setUserFun( My_fp ausrfun )
{
  //checkSet = checkSet+1;
  usrfun = ausrfun;
}



// ------------------------------ reliability part ---------------------------

int
SnoptProblem::findDesignPoint()
{


//---------------------define variables ------------------------------------------------------


//	 this->theReliabilityDomain = passedReliabilityDomain;

	int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
	
	int j;
//	int k;
	int zeroFlag;
	Vector dummy(numberOfRandomVariables);
	dummy.Zero();

	reliability_x = dummy;
	u = dummy; 
	int i;
	for(i=0;i<n;i++) x[i]=u(i);

//	Vector u_old(numberOfRandomVariables);
	uSecondLast = dummy;
//	Vector uNew(numberOfRandomVariables);

	alpha = dummy;
	gamma = dummy;
	alphaSecondLast= dummy;
	double gFunctionValue = 1.0;
//	double gFunctionValue_old = 1.0;
	Vector gradientOfgFunction(numberOfRandomVariables);
	gradientInStandardNormalSpace = dummy;
//	Vector gradientInStandardNormalSpace_old(numberOfRandomVariables);
	double normOfGradient =0;
//	double stepSize;
	
	Matrix jacobian_x_u(numberOfRandomVariables,numberOfRandomVariables);
	int evaluationInStepSize = 0;
	int result;
	
	theGFunEvaluator->initializeNumberOfEvaluations();  //	theGFunEvaluator->numberOfEvaluations = 0;

	
	// Prepare output file to store the search points
	
	if (this->outputFile ==0)
		outputFile = new ofstream( fileNamePrint, ios::out );


	if (printFlag == 0) {
		*outputFile << "The user has not specified to store any search points." << endln;
		*outputFile << "This is just a dummy file. " << endln;
	}

	
//------------------------ start point ------------------------------------------------	
	if (startAtOrigin)
	  u.Zero();
	else {

		// Get starting point
	  theReliabilityDomain->getStartPoint(reliability_x);


		// Transform starting point into standard normal space
		/*
		result = theProbabilityTransformation->set_x(reliability_x);  //return theProbabilityTransformation->(*x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not set x in the xu-transformation." << endln;
			return -1;
		}


		result = theProbabilityTransformation->transform_x_to_u();
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from x to u." << endln;
			return -1;
		}
		u = theProbabilityTransformation->get_u();
		*/
		result = theProbabilityTransformation->transform_x_to_u(reliability_x, u);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from x to u." << endln;
			return -1;
		}

		for(int i=0; i<n; i++) x[i] = u(i);
	}  //if (startPoint != 0) 

	
	
	//--------------- get Gfirst ---------------------------------------

		result = theGFunEvaluator->runGFunAnalysis(reliability_x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not run analysis to evaluate limit-state function. " << endln;
			return -1;
		}
		result = theGFunEvaluator->evaluateG(reliability_x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not tokenize limit-state function. " << endln;
			return -1;
		}
		
		gFunctionValue = theGFunEvaluator->getG();
		
		Gfirst = gFunctionValue;



// --------------------------------------------------------------------- 
//                                                                     |
	integer Cold = 0, Basis = 1, Warm = 2;                    //       |
    this->solve( Cold );                                      //       |
	opserr<<"SNOPTProblem::solve(Cold) is called! ...done"<<endln;  // | 
//                                                                     |      
// --------------- after solving the problem ---------------------------
//	this->getHessian();

	for(i=0; i<n; i++) u(i) = x[i];
	
	// Transform from u to x space

	/*
	result = theProbabilityTransformation->set_u(u);
	if (result < 0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " could not set u in the xu-transformation." << endln;
	  return -1;
	}
	
	result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
	if (result < 0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " could not transform from u to x and compute Jacobian." << endln;
	  return -1;
	}
	*/
	//reliability_x = theProbabilityTransformation->get_x();
	//jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();

	result = theProbabilityTransformation->transform_u_to_x(u, reliability_x);	
	if (result < 0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " could not transform from u to x and compute Jacobian." << endln;
	  return -1;
	}	
	result = theProbabilityTransformation->getJacobian_x_to_u(reliability_x, jacobian_x_u);	
	if (result < 0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " could not transform from u to x and compute Jacobian." << endln;
	  return -1;
	}	

	
	// -------------------------- compute alpha ----------------------------
	
	gradientOfgFunction = theGradGEvaluator->getGradG();
	
	
	// Check if all components of the vector is zero
	zeroFlag = 0;
	for (j=0; j<gradientOfgFunction.Size(); j++) {
	  if (gradientOfgFunction[j] != 0.0) {
	    zeroFlag = 1;
	  }
	}
	if (zeroFlag == 0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " all components of the gradient vector is zero. " << endln;
	  return -1;
	}
	
	
	// Gradient in standard normal space
	//	gradientInStandardNormalSpace_old = gradientInStandardNormalSpace;
	gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;
	
	
	
	// Compute the norm of the gradient in standard normal space
	normOfGradient = gradientInStandardNormalSpace.Norm();
	
	
	// Check that the norm is not zero
	if (normOfGradient == 0.0) {
	  opserr << "SNOPT::doTheActualSearch() - " << endln
		 << " the norm of the gradient is zero. " << endln;
	  return -1;
	}
	
	
	// Compute alpha-vector
	alpha = gradientInStandardNormalSpace *  ( (-1.0) / normOfGradient );
	
	
	// Print the design point to file, if desired
	int iii;
	if (printFlag != 0) {
	  outputFile->precision(16);
	  if (printFlag == 3) {
	    
	    outputFile->setf(ios::scientific, ios::floatfield);
	    *outputFile<<endln;
	    for (iii=0; iii<reliability_x.Size(); iii++) {
	      *outputFile<<reliability_x(iii)<<endln;
	    } //for
	  } //if
	  else if (printFlag == 4) {
	    //					static ofstream *outputFile( fileNamePrint, ios::out );
	    outputFile->setf(ios::scientific, ios::floatfield);
	    *outputFile<<endln;
	    for (iii=0; iii<u.Size(); iii++) {
	      *outputFile<<u(iii)<<endln;
	    } //for
	  } //else
	} //if
	
	
	// -------------------- Compute the gamma vector ----------------------------
	
	MatrixOperations theMatrixOperations(jacobian_x_u);
	
	result = theMatrixOperations.computeTranspose();
	if (result < 0) {
	  opserr << "snopt::doTheActualSearch() - " << endln
		 << " could not compute transpose of jacobian matrix. " << endln;
	  return -1;
	}
	Matrix transposeOfJacobian_x_u = theMatrixOperations.getTranspose();
	
	Matrix jacobianProduct = jacobian_x_u * transposeOfJacobian_x_u;
	
	Matrix D_prime(numberOfRandomVariables,numberOfRandomVariables);
	for (j=0; j<numberOfRandomVariables; j++) {
				D_prime(j,j) = sqrt(jacobianProduct(j,j));
	}
	
	Matrix jacobian_u_x(u.Size(), reliability_x.Size());
	theProbabilityTransformation->getJacobian_u_to_x(u, jacobian_u_x);
	
	Vector tempProduct = jacobian_u_x ^ alpha;
	
	gamma = D_prime ^ tempProduct;
	
	gFunctionValue = theGFunEvaluator->getG();
	Glast = gFunctionValue;
	/*			if (fabs(gFunctionValue) > 1.0e-4){ 
	  
	opserr<<"SNOPT can not find design point, final g is "<<Glast<<endln;
	return -1; } // not converge}
	*/
	numberOfEvaluations = theGFunEvaluator->getNumberOfEvaluations();
	outputFile->close();
	return 1;
	//	}
	
	
};



const Vector &
SnoptProblem::get_x()
{
	
	return reliability_x;
}

const Vector &
SnoptProblem::get_u()
{
	for (int i=0;i<n;i++) u(i)=x[i];
	return u;
}

const Vector &
SnoptProblem::get_alpha()
{
  return alpha;
}

const Vector &
SnoptProblem::get_gamma()
{
  static Vector(res);
  res =gamma*(1.0/gamma.Norm());
  return res;
}

int
SnoptProblem::getNumberOfSteps()
{
  return (ii-1);
}

const Vector &
SnoptProblem::getSecondLast_u()
{

  opserr<<"SnoptProblem::getSecondLast_u() is called!  "<<endln;
  return uSecondLast;
	
}

const Vector &
SnoptProblem::getSecondLast_alpha()
{
  opserr<<"SnoptProblem::getSecondLast_alpha() is called! "<<endln;
  return alphaSecondLast;
  
}

const Vector &
SnoptProblem::getLastSearchDirection()
{	
  return gradientInStandardNormalSpace;
}

double
SnoptProblem::getFirstGFunValue()
{
  return Gfirst;
}

double
SnoptProblem::getLastGFunValue()
{
  return Glast;
}


const Vector &
SnoptProblem::getGradientInStandardNormalSpace()
{
  return gradientInStandardNormalSpace;
}



int
SnoptProblem::getNumberOfEvaluations()
{
	return numberOfEvaluations;
}


// Quan and Michele

int SnoptProblem::setStartPt(Vector * pStartPt)
{
  //startPoint->addVector(0.0,(*pStartPt),1.0);
	return 0;
}




int SnoptProblem::setUSecondLast(const Vector &theU)
{
	if (theU.Size() != uSecondLast.Size()) {
		opserr<<"Warning in SnoptProblem::setUSecondLast, different vector size ! "<<endln;;
		opserr<<"theU.Size()"<<theU.Size()<<endln;
		opserr<<"uSecondLast.Size()"<<uSecondLast.Size()<<endln;

	//	return -1;
	}

	uSecondLast.addVector(0.0, theU, 1.0);
	return 1;
}


int SnoptProblem::setU(const Vector &theU)
{
	if (theU.Size() != u.Size()) {
		opserr<<"Warning in SnoptProblem::setU, different vector size ! "<<endln;
		//return -1;
	}

	u.addVector(0.0, theU, 1.0);
	return 1;
}



int SnoptProblem::setAlphaSecondLast(const Vector &theAlpha){

	if (theAlpha.Size() != u.Size()) {
		opserr<<"Warning in SnoptProblem::theAlpha, different vector size ! "<<endln;
		opserr<<"theAlphaSize()"<<theAlpha.Size()<<endln;
		opserr<<"u.Size()"<<u.Size()<<endln;

//		return -1;
	}

	alphaSecondLast.addVector(0.0, theAlpha, 1.0);
	return 1;
}

 






// Quan
double * SnoptProblem::getHessian(){



//* ----- Quan Gu 
//!     ------------------------------------------------------------------
//!     Retrieve the approximate Hessian from memory.
//!     The option 'Hessian full memory' must have been selected.
//!     ------------------------------------------------------------------
 




    integer errors, lenh=10000;
	doublereal * h__ = new doublereal[lenh]; 
	integer nnh;
	
    snreth_(&errors, &lenh, h__, &nnh, cw, &lencw, iw, &leniw, rw, &lenrw, 8*500);


	// Open output file and start writing to it
	ofstream outputFile3( "ApproximateHessian.out", ios::out );



//*  --
  outputFile3 << "=================================================\n";
  outputFile3 << "Approximate Hessian is: " <<  endl;

     int jthcol = 1;

	 
/*      jthcol = 1
      do   j = 1, nnH
         jthcol = jthcol + j - 1
         write (iSumm, 1000) (H(k), k=jthcol,jthcol+j-1)
      end do  */

     int j;
	 for(j = 1; j<=nnh; j++){
         jthcol = jthcol + j - 1;
		 outputFile3 << endl;
	     for (int k=jthcol; k<=jthcol+j-1;k++)
		    outputFile3<< h__[k-1]<<"     ";

	 }
	  outputFile3 << endl;
      outputFile3<< " ========== solution x is: ================="  <<endl; 
	  for (j=0;j<n; j++)
		  outputFile3<< x[j]<<"   ";

	  outputFile3<<endl;

      outputFile3 << " \n === lagrangian multiplier is: ================\n"  ; 
	  outputFile3 << Fmul[1] <<endl;


      


       
  outputFile3 << "=================================================\n";
  

  // Clean up
  outputFile3.close();

  return h__;
}


