#include <stdio.h>
#include <string.h>
#include "snopt.h"
#include "ToyFunction.h"
#include "SnoptProblem.h"

extern SnoptProblem * theSNOPT;// = new testaa();

int toyusrf_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
	     integer    iu[],    integer *leniu,
	     doublereal ru[],    integer *lenru )
{

double gFunctionValue = 1.0;
Vector gradientOfgFunction(*n);
double normOfGradient =0;
int result,i;
Matrix jacobian_x_u(*n,*n);


	F[0] =0;
	for (i=0; i< *n; i++) F[0] += 0.5*x[i]*x[i];  //F0=1/2*u'*u

//	F[1] = (x[0]-3)*(x[0]-3) + (x[1]-6)*(x[1]-6)+x[2]*x[2]-6;

//---------------- first transfer x(u) to reliability_x----------


	theSNOPT->setUSecondLast(theSNOPT->u);
	for (i=0; i<*n; i++) theSNOPT->u(i)=x[i];

	/*
	result = theSNOPT->theProbabilityTransformation->set_u(theSNOPT->u);
	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
			<< " could not set u in the xu-transformation." << endln;
		return -1;
	}

	result = theSNOPT->theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "SNopt::doTheActualSearch() - " << endln
				<< " could not transform from u to x and compute Jacobian." << endln;
			return -1;
		}
	theSNOPT->reliability_x = theSNOPT->theProbabilityTransformation->get_x();
	jacobian_x_u = theSNOPT->theProbabilityTransformation->getJacobian_x_u();
	*/
	theSNOPT->theProbabilityTransformation->transform_u_to_x(theSNOPT->u, theSNOPT->reliability_x);
	theSNOPT->theProbabilityTransformation->getJacobian_x_to_u(theSNOPT->reliability_x, jacobian_x_u);


	result = theSNOPT->theGFunEvaluator->runGFunAnalysis(theSNOPT->reliability_x);

	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
		<< " could not run analysis to evaluate limit-state function. " << endln;
		//return -1;
		(*Status) = -1;
	}
	//else

	result = theSNOPT->theGFunEvaluator->evaluateG(theSNOPT->reliability_x);
	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
			<< " could not tokenize limit-state function. " << endln;
		return -1;
	}
//	gFunctionValue_old = gFunctionValue;
	gFunctionValue = theSNOPT->theGFunEvaluator->getG();
//------------------------------------------
	F[1]=gFunctionValue;
 
	if ((theSNOPT->outputFile !=0)&&(theSNOPT->printFlag ==2))
	{ 
		*(theSNOPT->outputFile)<<endln;
		 for(int i=0; i<*n;i++)
			*(theSNOPT->outputFile)<<theSNOPT->u[i]<<endln;
	      	*(theSNOPT->outputFile)<<endln;
	}
	else if ((theSNOPT->outputFile !=0)&&(theSNOPT->printFlag ==1))
	{ 
		*(theSNOPT->outputFile)<<endln;
 		 for(int i=0; i<*n;i++)
			*(theSNOPT->outputFile)<<theSNOPT->reliability_x[i]<<endln;
	      	*(theSNOPT->outputFile)<<endln;
	}

	return 0;
}

int toyusrfg_( integer    *Status, integer *n,    doublereal x[],
	       integer    *needF,  integer *neF,  doublereal F[],
	       integer    *needG,  integer *neG,  doublereal G[],
	       char       *cu,     integer *lencu,
	       integer    iu[],    integer *leniu,
	       doublereal ru[],    integer *lenru )
{
  //==================================================================
  // Computes the nonlinear objective and constraint terms for the toy
  // problem featured in the SnoptA users guide.
  // neF = 3, n = 2.
  //
  //   Minimize     1/2*x'*x
  //
  //   subject to   g(x)=0 i.e., (x1-3)^2+(x2-6)^2+x3^2=0
  //
  // The triples (g(k),iGfun(k),jGvar(k)), k = 1:neG, define
  // the sparsity pattern and values of the nonlinear elements
  // of the Jacobian.
  //==================================================================



double gFunctionValue = 1.0;
Vector gradientOfgFunction(*n);
double normOfGradient =0;
int result,i;
Matrix jacobian_x_u(*n,*n);

if( *needF > 0 ) {
	F[0] =0;
	for (i=0; i< *n; i++) F[0] += 0.5*x[i]*x[i];  //F0=1/2*u'*u

/* Quan debug only----.
  opserr.precision(16);
  opserr<<"========= x[] is:"<<endln<<endln<<endln<<endln;
  for (i=0; i< *n; i++) opserr<< x[i]<< "  ";
  opserr<<endln;
// ----Quan debug only.
*/

//	F[1] = (x[0]-3)*(x[0]-3) + (x[1]-6)*(x[1]-6)+x[2]*x[2]-6;

//---------------- first transfer x(u) to reliability_x----------
	theSNOPT->setUSecondLast(theSNOPT->u);
	for (i=0; i<*n; i++) theSNOPT->u(i)=x[i];

	/*
	result = theSNOPT->theProbabilityTransformation->set_u(theSNOPT->u);
	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
			<< " could not set u in the xu-transformation." << endln;
		return -1;
	}

	result = theSNOPT->theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "SNopt::doTheActualSearch() - " << endln
				<< " could not transform from u to x and compute Jacobian." << endln;
			return -1;
		}
	theSNOPT->reliability_x = theSNOPT->theProbabilityTransformation->get_x();
	jacobian_x_u = theSNOPT->theProbabilityTransformation->getJacobian_x_u();
	*/
	theSNOPT->theProbabilityTransformation->transform_u_to_x(theSNOPT->u, theSNOPT->reliability_x);
	theSNOPT->theProbabilityTransformation->getJacobian_x_to_u(theSNOPT->reliability_x, jacobian_x_u);



	//*Quan debug -----
	opserr<<"-------------in physical space, x is:----------------"<<endln;
	opserr<<theSNOPT->reliability_x<<endln;
	opserr<<"-----------------------------------------------------"<<endln;
//*/

	result = theSNOPT->theGFunEvaluator->runGFunAnalysis(theSNOPT->reliability_x);


	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
		<< " could not run analysis to evaluate limit-state function. " << endln;
		//return -1;
		(*Status) = -1;
	}
	result = theSNOPT->theGFunEvaluator->evaluateG(theSNOPT->reliability_x);
	if (result < 0) {
		opserr << "SNopt::doTheActualSearch() - " << endln
			<< " could not tokenize limit-state function. " << endln;
		return -1;
	}
//	gFunctionValue_old = gFunctionValue;
	gFunctionValue = theSNOPT->theGFunEvaluator->getG();
//------------------------------------------
	F[1]=gFunctionValue;

	// Quan debug only----.
    //opserr.precision(16);
	//opserr<<"F is: " <<F[1]<<endln;

} //if *needF





//------------------ grad of [F;G] ---------------------------------------------

  
  if( *needG > 0 ){
	  if (*needF<=0) {opserr<<"needG but no needF !!  in snopt userfunction calling! "<<endln;/* exit(-1);*/}
	  
	  
	  *neG = 0;
	  for (i=0; i<*n; i++) {G[*neG] = x[i];*neG = *neG + 1;}



//----------------- grad of G ------------------------------
		// Gradient in original space
		result = theSNOPT->theGradGEvaluator->computeGradG(gFunctionValue,theSNOPT->reliability_x);
		if (result < 0) {
			opserr << "SNOPT::doTheActualSearch() - " << endln
				<< " could not compute gradients of the limit-state function. " << endln;
			return -1;
		}
		gradientOfgFunction = theSNOPT->theGradGEvaluator->getGradG();

 /*Quan debug only ----.
	  opserr<<"gradient G is:"<<endln;
	  for (i=0; i<(*n); i++) 	opserr <<gradientOfgFunction[i]<<", ";
	  opserr<<endln;
*/// ---	
//

  

		// Check if all components of the vector is zero
		int zeroFlag = 0;
		for (int j=0; j<gradientOfgFunction.Size(); j++) {

			if (gradientOfgFunction[j] != 0.0) {
				zeroFlag = 1;
			}
		}
		if (zeroFlag == 0) {
			opserr << "SNopt::doTheActualSearch() - " << endln
				<< " all components of the gradient vector is zero. " << endln;
			return -1;
		}
	   

		// Gradient in standard normal space
		//gradientInStandardNormalSpace_old = gradientInStandardNormalSpace;
		theSNOPT->gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;


		// Compute the norm of the gradient in standard normal space
		normOfGradient = theSNOPT->gradientInStandardNormalSpace.Norm();


		// Check that the norm is not zero
		if (normOfGradient == 0.0) {
			opserr << "SNopt::doTheActualSearch() - " << endln
				<< " the norm of the gradient is zero. " << endln;
			return -1;
		}
		
	//	opserr.setPrecision(16);
	//	opserr<<"Grad G is : ";

      for (i=0; i<*n; i++) {G[*neG] = theSNOPT->gradientInStandardNormalSpace[i];   *neG = *neG + 1;}

	  theSNOPT->setAlphaSecondLast(theSNOPT->get_alpha());
	  theSNOPT->alpha.addVector(0.0,theSNOPT->gradientInStandardNormalSpace, -1.0/normOfGradient);



/*    // iGfun[*neG] = 1
    // jGvar[*neG] = 1
    *neG = 0;
    G[*neG] = x[0];

    // iGfun[*neG] = 1
    // jGvar[*neG] = 2
    *neG = *neG + 1;
    G[*neG] = x[1];

    // iGfun[*neG] = 1
    // jGvar[*neG] = 3
    *neG = *neG + 1;
    G[*neG] = x[2];

    // iGfun[*neG] = 2
    // jGvar[*neG] = 1
    *neG = *neG + 1;
	//---------------
    G[*neG] = 2*x[0]-6;

    // iGfun[*neG] = 2
    // jGvar[*neG] = 2
    *neG = *neG + 1;
    G[*neG] = 2*x[1]-12;

    // iGfun[*neG] = 2
    // jGvar[*neG] = 3
    *neG = *neG + 1;
    G[*neG] = 2*x[2];
    *neG = *neG + 1;  */
  }

  
	if ((theSNOPT->outputFile !=0)&&(theSNOPT->printFlag ==2))
	{ 
		*(theSNOPT->outputFile)<<endln;
		 for(int i=0; i<*n;i++)
			*(theSNOPT->outputFile)<<theSNOPT->u[i]<<endln;
	      	*(theSNOPT->outputFile)<<endln;
	}
	else if ((theSNOPT->outputFile !=0)&&(theSNOPT->printFlag ==1))
	{ 
		*(theSNOPT->outputFile)<<endln;
 		 for(int i=0; i<*n;i++)
			*(theSNOPT->outputFile)<<theSNOPT->reliability_x[i]<<endln;
	      	*(theSNOPT->outputFile)<<endln;
	}

  return 0;
}
