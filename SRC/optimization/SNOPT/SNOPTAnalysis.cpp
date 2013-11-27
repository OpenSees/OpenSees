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
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include "snopt.h"
#include "SNOPTAnalysis.h"

using namespace std;



// -------------------------------------------------------------------------------------------------
//// ---- Note: everything here is 'scaled'. So  need scale back when deal with designVariables.   |
// -------------------------------------------------------------------------------------------------

SNOPTAnalysis::~SNOPTAnalysis()
{
	if (fileNamePrint !=0) delete [] fileNamePrint;
	if (iAfun !=0) delete [] iAfun;
	if (jAvar !=0) delete [] jAvar;
	if (iGfun !=0) delete [] iGfun;
	if (jGvar !=0) delete [] jGvar;
	if (gradient !=0) delete [] gradient;
	if (scales !=0) delete [] scales;
	if (x !=0) delete [] x;
	if (xlow !=0) delete [] xlow;
	if (xupp !=0) delete [] xupp;
	if (xmul !=0) delete [] xmul;
	if (xstate !=0) delete [] xstate;
	if (xnames !=0) delete [] xnames;
	if (Fnames !=0) delete [] Fnames;
	if (tclFileToRun !=0) delete [] tclFileToRun;



}


SNOPTAnalysis::SNOPTAnalysis(int passedMaxNumberOfIterations, 
					int pPrintFlag,
					char *pFileNamePrint,
					char * probType, 
					OptimizationDomain * passedOptimizationDomain,
					Tcl_Interp *passedTclInterp,
					char * pTclFileToRun
					): SNOPTCLASS()

{

// --------------------- step 1. SNOPTCLASS::SNOPTCLASS() auto called -------------------------

/*//  iSpecs(0), iSumm(6), iPrint(0), initCalled(0)
initCalled=0;
iSpecs = 0;
iSumm=6;
iPrint=0;
	  
// ------ please note in this construction, only set memory and init value to 0, do not call objectiveFunction or constraintFunction to set value !!
	  




  init2zero();

  // Nothing incremented yet
  fortranStyleObj = 0;
  fortranStyleAG  = 0;

  iSpecs =  4;
  iSumm  = 6;
  //  iPrint = 15;

  // Create temporary memory for the call to sninit_.
  // Lengths must all be >= 500.
  lencw  = 500;
  leniw  = 500;
  lenrw  = 500;
  this->alloc( 500, 500, 500 );

  // sninit_ "undefines" the optional parameters

  this->init();

*/
  
// ---------------------- step 2 save optimization pointer ------------------------------------

 // if probType =="optimization" then
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	printFlag						= pPrintFlag;
	numberOfEvaluations =0;
	fileNamePrint = new char[256];
	if (printFlag != 0) {
		strcpy(fileNamePrint,pFileNamePrint);
	}
	else {
		strcpy(fileNamePrint,"SNOPPoints.out");
	}


// ---------------------- step 3 snopt construction for optimization -------------------------------
	this->theOptimizationDomain = passedOptimizationDomain;
	this->theConstraintFunction=theOptimizationDomain->getConstraintFunctionPtr(1);
	this->theObjectiveFunction=theOptimizationDomain->getObjectiveFunctionPtr(1);
	this->theTclInterp = passedTclInterp;

 
//--------------------- step 3.1 for snopt memory alloc -------------------------------------
  // Allocate and initialize;
	int numOfMaxIter = this->maxNumberOfIterations;  // this number need to be input
	
	
// --- Quan control  from here....	
	this->n     =  theOptimizationDomain->getNumberOfDesignVariables();   // n need to be got from DesignVariables x1 x2 x3 ..
	if (theConstraintFunction ==0) this->neF =1;   // no constraint function
	
	else
		this->neF   =  this->theConstraintFunction->getNumOfConstraintFunctions()+1;// in this case, only 1  objectivefunction

// ------------------- A --------------
///*  Quan 1
	Vector * myA;
	Matrix * myB;
    
	myA  = theObjectiveFunction->getLinearAdditionPtr();   // may or may not needed
	
	if (theConstraintFunction !=0)
		myB  = theConstraintFunction->getLinearAdditionPtr();
	else myB=0;

	if ((myA !=0) || (myB !=0)) {
		opserr<< "linearAddition is not zero"<<endln;

// ---- cheap check for format --
		if (((myA !=0) &&(myA->Size() !=n)) ||((myB !=0) &&(myB->noCols() !=n))){
			opserr<<"fatal: A dimention is not compatible with design variables number!"<<endln;
		    exit(-1);
		}


		
		if (myB !=0)
			lenA = n + myB->noCols() * myB->noRows();
		else  lenA = myA->Size() ; //n
		
		neA = lenA; // correct?



		this->iAfun = new integer[lenA];
		this->jAvar = new integer[lenA];
		this->A  = new doublereal[lenA];
	
	  // ---- fill A --
		int i;
		for (i=0; i<n; i++)
			if (myA !=0) A[i] = (*myA)[i];
			else A[i] = 0;

		if (myB !=0) {
			for (int k=0; k<myB->noRows(); k++){
			  for (int j=0; j<n; j++)
				  {A[i] = (*myB)(k,j); i++;}
			} 
		} //if myB !=0
	}//if
	else{



	// Quan, for snopt interface only
	// YONG, DEBUG HERE.
		lenA=this->n*this->neF;
		iAfun = new integer[lenA];
		jAvar = new integer[lenA];
		A  = new doublereal[lenA];
		this->neA=0;
	}//if ((myA !=0) || (myB !=0))







  


// ------------------ G -----------------  

  this->lenG   = n*neF; //// for optimization only

  this->needG=false;
  if (theObjectiveFunction->isGradientProvided())
	  this->needG = true;
  if (theConstraintFunction !=0){
	  if (theConstraintFunction->isGradientProvided())
		this->needG = true;
  }

  
//  if (needG){
	  this->iGfun = new integer[lenG];
	  this->jGvar = new integer[lenG];
	  gradient = new double [lenG];
	  for(int i=0; i<lenG;i++)	  gradient[i]=0;
//  }
//  else gradient=0;

// --------------- x and F -------------------
	
  this->scales = new doublereal[n];

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

//--------------------- step 3.2 set x and problem bound for optimization --------------------------------------

  this-> ObjRow = 0;
  this-> ObjAdd = 0;


  // --- Set the x and its upper and lower bounds.
  
  DesignVariable * theDesignVariable;
  for (int i=0;i<n;i++){
   theDesignVariable=theOptimizationDomain->getDesignVariablePtr(i+1);


   // ------------- scale here --------------------------------
	scales[i] = theDesignVariable->getScale();

	x[i]      = theDesignVariable->getValue() / scales[i];
	xlow[i]   = theDesignVariable->getLowerBound()/ scales[i];
	xupp[i]   = theDesignVariable->getUpperBound()/ scales[i];
	xstate[i] = 0; //x[i]/ scales[i];
	xmul[i]   = theDesignVariable->getXMultiplier() / scales[i]; //0
  } //for
 
  // --- Set the F and its upper and lower bounds ----------- . 
  
  
  Flow[0] = theObjectiveFunction->getLowerBound();
  Fupp[0] = theObjectiveFunction->getUpperBound();


  Fstate[0] = theObjectiveFunction->getState();
  Fmul[0]   = theObjectiveFunction->getMultiplier();



  	  Vector * lowerBound = 0;
	  Vector * upperBound = 0;
	  Vector * multiplier =0; 
	  Vector * state = 0; 
  if (theConstraintFunction !=0){
	   lowerBound = theConstraintFunction->getLowerBoundPtr();
	   upperBound = theConstraintFunction->getUpperBoundPtr();
	   multiplier = theConstraintFunction->getMultiplierPtr(); 
	   state = theConstraintFunction->getStatePtr(); 
  }
  for(int i=1; i< neF; i++){
	  if (lowerBound !=0) 		Flow[i] = (*lowerBound)[i-1];
	  else                      Flow[i] = -1e20;;

	  if (upperBound !=0) 	Fupp[i] = (*upperBound)[i-1];
	  else                  Fupp[i] = 1e20;

	  if (multiplier !=0) 	Fstate[i] = (*multiplier)[i-1];
	  else                  Fstate[i] = 0; 

	  if (state !=0)     Fmul[i]   = (*state)[i-1];
	  else               Fmul[i]   = 0;
  }
  
//--------------------- step 3.4 set G sequence -----------------------------------------

   
// should be: (0,0) (0,1) (0,2) (0,3)..  (1,0) (1,1) (1,2) (1,3)..


  if (needG){
	  neG=0;
	  for (int i=0; i<neF;i++){
		  for(int j =0; j<n; j++){
			 iGfun[neG] = i;
			 jGvar[neG] = j;
			 neG++;
			}
	  }
  }

// ==== Quan oct 24 for debug ungradient case


//--------------------- step 3.5 set others  -----------------------------------------

   //this->setSpecFile    ( "sntoya.spc" ); // uncommented by YONG
   //this->setIntParameter( "Derivative option", 1 );// uncommented by YONG
   //this->setIntParameter( "Major Iteration limit", maxNumberOfIterations );// uncommented by YONG
    this->setProbName   ( probType );
	this->ii=0;
	if (pTclFileToRun !=0){
		tclFileToRun = new char[30];
		strcpy( tclFileToRun, pTclFileToRun);  // to run this part only if structure analysis. otherwise tclFileToRun=0.
	}
	else tclFileToRun = 0;


	temp = 0;



}



int SNOPTAnalysis::runOptAnalysis (OptimizationDomain *passedOptimizationDomain){

	this->theOptimizationDomain = passedOptimizationDomain;	
//	int result;
	
	// Prepare output file to store the search points
	ofstream outputFile2( fileNamePrint, ios::out );
	if (printFlag == 0) {
		outputFile2 << "The user has not specified to store any search points." << endln;
		outputFile2 << "This is just a dummy file. " << endln;
	}



// 1. ----------------------------------------------------------------- 
//                                                                     |
	if (!needG)  {                                            //       |
	  this->usrfun= ( My_fp ) toyOptusrf_ ;  // Sets the usrfun that supplies F. |
	  this->setPrintFile   ( "Toy0.out" );                    //       |
	  computeJac ();                                          //       |
      setIntParameter( "Derivative option", 0 );	          //       |
	}                                                         //       |
	else {                                                    //       |

	  this->usrfun= ( My_fp) toyOptusrfg_ ; // Sets the usrfun that supplies F,G|
	  this->setPrintFile   ( "Toy1.out" );                    //       |
      setIntParameter( "Derivative option", 1 );	          //       |
    }                                                         //       |
	integer Cold = 0, Basis = 1, Warm = 2;                    //       |
    this->solve( Cold );                                      //       |
	opserr<<"SNOPTProblem::solve(Cold) is called! ...done"<<endln;  // | 
// --------------------------------------------------------------------- 
	
	      
// 2. -------------- after solving the problem  judge whether this point is a local minimum, if so, output into file ---
	
	       

			numberOfEvaluations = theObjectiveFunction->getNumberOfEvaluations();

// Print the design point to file, if desired
		int i;
		if (printFlag != 0) {
				outputFile2.setf(ios::scientific, ios::floatfield);
				outputFile2<<"The Minimum point X is found as followings:"<<endln;
				for (i=0; i<n; i++) {
						outputFile2<<x[i]*scales[i]<<endln;  //scale back
				} //for
		} //if

	
			return 1;

  
	};


int
SNOPTAnalysis::getNumberOfSteps()
{
	return (ii-1);
}


int
SNOPTAnalysis::getNumberOfEvaluations()
{
	return numberOfEvaluations;
}



double * SNOPTAnalysis::getScaledX(){
	return x;   // here is scaled value! mean: x_i/x_i_0
}

double * SNOPTAnalysis::getUnScaledX(){
	if (temp ==0) temp = new double [n];

	for(int i=0; i<n;i++)
		temp[i]=x[i]*scales[i];

	return temp;   // here is UnScaled value!  mean: x_i, in physical space,(user input in tcl)
}

int SNOPTAnalysis::updateXFG(double * newScaledX){

//--------------- get objective function and gradient from ---------------------------------------
	 
// 1. --- reset
// But, how to revertToStart() for both resp and sens , 'reset' need extend to all ele and mat.

	
	if (tclFileToRun != 0) {     // structure analysis
		char theRevertToStartCommand[10] = "reset";
		Tcl_Eval( theTclInterp, theRevertToStartCommand );
	}
	
// 2. ---- update x
   for(int i=0; i<n;i++){
	  
	  DesignVariable * theDesignVariable = theOptimizationDomain->getDesignVariablePtr(i+1);


    //	if(newScaledX[i] ==0){
	//	int aaa=1;
	//}
      theDesignVariable->update(newScaledX[i]*scales[i]);//========== scale here =========================
 
	} //for 

/*   	opserr<<"x is:";
	opserr.setPrecision(16);
	for ( i=0; i<n; i++)
		opserr<<newScaledX[i]*scales[i]<<" , ";
	opserr<<endln;
*/
	// quan debug
//	if (fabs(newScaledX[0]*scales[0]-33.89955589282231)< 1e-14) {
//	opserr<<"bug here"<<endln;

	
//	}



// 3. ---  clear all the recorder file ----------------------------------

   // maybe need change recorder class.

   
// 4, --- here run opensees again ?  yes! file name is in "tclFileToRun" ----------
	// a. run "wipeAnalysis"
    // b. run tclFileToRun'

    // ========================= -----------------------------------------------------------------------------------
	// ========== should add something in the following part !! before calling updateSFG, should pass whether needG
    // ========== when needF but not needG, should automatically remove (skip) sensitiity part. do parse analysis
    // ========== guquan Oct 12 2005 --------------------------------------------------------------------------------


   if (tclFileToRun != 0) {     // structure analysis

/*   user to write this in tcl file 
	   char theWipeAnalysisCommand[15] = "wipeAnalysis";
		if (Tcl_Eval( theTclInterp, theWipeAnalysisCommand ) != TCL_OK ) {
			opserr<<"Fatal: can not wipeAnalysis in SNOPTAnalysis"<<endln;
			exit(-1);
		}
*/

			
/*	
	if( Tcl_EvalFile(theTclInterp, tclFileToRun) !=TCL_OK ) { // -- divergence case ----
			opserr<<"Fatal: can not run tcl file  "<<tclFileToRun<< " in SNOPTAnalysis"<<endln;
			//exit(-1);
			
			return -1;
		}   
 */ 
		double result1 = -1.0;
		char theAnalyzeCommand[50];
		sprintf(theAnalyzeCommand,"[source %s]",tclFileToRun);
		Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result1);
	
		if ( result1 <0) {    // need to change back
			opserr<<"Fatal: can not run tcl file  "<<tclFileToRun<< " in SNOPTAnalysis"<<endln;
			return -1; 
		}


      
    }




// 5. --- call objectiveFunction->update(); and constraintFunction->update()

	theObjectiveFunction->update();
	if (theConstraintFunction !=0)
		theConstraintFunction->update();

// 6. --- F = (objective & constraint)->getValue(); G= (objective & constraint) -> getGradient();
    
	F[0] = theObjectiveFunction->getValue(); //ok 
	Vector * temp = 0;

	if (theConstraintFunction !=0) {
	   temp = theConstraintFunction->getValuePtr();
	   for(int i=1; i<neF; i++)  F[i] = (*temp)[i-1];
	}
	

	if (needG){
		if (theObjectiveFunction->isGradientProvided()) {
			Vector * temp1 = theObjectiveFunction->getGradientPtr();
			for(int i=0; i<n; i++){
				gradient[i] =  (*temp1)[i]*scales[i];   //========= scale here, change to du/dSita, instead of du/dx
			}
		}
		else{ // do nothing		
			};

		if (theConstraintFunction !=0){
			if (theConstraintFunction->isGradientProvided()) {
				Matrix * temp2 = theConstraintFunction->getGradientPtr();  
				int i=n;

				for(int j=0; j<neF-1; j++){
					for(int k=0; k<n; k++){
						gradient[i]= (*temp2)(j,k)*scales[k];   // ====== scale here, change to du/dSita, instead of du/dx
						i++;
					}
				}
			}
			else{ // do nothing		
				};

		} // if constraint exist
	
	} //if needG
	

	return 0;
}