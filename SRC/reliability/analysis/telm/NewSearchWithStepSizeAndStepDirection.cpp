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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.5 $
// $Date: 2008-10-22 16:41:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewSearchWithStepSizeAndStepDirection.cpp,v $
                                                                     


#include <NewSearchWithStepSizeAndStepDirection.h>
#include <FindDesignPointAlgorithm.h>
#include <ReliabilityDomain.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <MatrixOperations.h>
#include <HessianEvaluator.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>
#include <GammaRV.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;



NewSearchWithStepSizeAndStepDirection::NewSearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					GradientEvaluator *passedGradGEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					HessianEvaluator *passedHessianEvaluator,
					ReliabilityConvergenceCheck *passedReliabilityConvergenceCheck,
					bool pStartAtOrigin,
					int pprintFlag,
					char *pFileNamePrint)
  :FindDesignPointAlgorithm(), theReliabilityDomain(passedReliabilityDomain)
{
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	theGFunEvaluator				= passedGFunEvaluator;
	theGradGEvaluator				= passedGradGEvaluator;
	theStepSizeRule					= passedStepSizeRule;
	theSearchDirection				= passedSearchDirection;
	theProbabilityTransformation	= passedProbabilityTransformation;
	theHessianEvaluator             = passedHessianEvaluator;
	theReliabilityConvergenceCheck  = passedReliabilityConvergenceCheck;
	startAtOrigin = pStartAtOrigin;
	printFlag						= pprintFlag;
	numberOfEvaluations =0;
	numberOfSensAna =0;
	if (printFlag != 0) {
		strcpy(fileNamePrint,pFileNamePrint);
	}
	else {
		strcpy(fileNamePrint,"searchpoints.out");
	}

	xinit=0;
	x=0;
	u=0;
	u_old=0;
	uSecondLast=0;
	uThirdLast=0;
	alpha=0;
	alpha_old=0;
	alphaSecondLast=0;
	alphaThirdLast=0;
	gradientInStandardNormalSpace=0;
	gradientInStandardNormalSpace_old=0;
	gamma=0;
	uSecondLast=0;
	searchDirection=0;
	gradientOfgFunction=0;
	jacobian_x_u =0;

	x_min=0;
	u_min=0;
	alpha_min=0;
	gradientInStandardNormalSpace_min=0;
	gradientOfgFunction_min=0;
	jacobian_x_u_min=0;

	output.open("Search.txt", ios::out);
	outputTMP.open( "finddespoint.txt", ios::out );

}

NewSearchWithStepSizeAndStepDirection::~NewSearchWithStepSizeAndStepDirection()
{

}


int
NewSearchWithStepSizeAndStepDirection::findDesignPoint()
{

  	check1_init=0.0;
	check2_init=0.0;
	check1_conv=0.0;
	check2_conv=0.0;

	// Set the reliability domain (a data member of this class)
	//theReliabilityDomain = passedReliabilityDomain;
	// Declaration of data used in the algorithm
	int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
	int j;
	int zeroFlag;
	Vector dummy(numberOfRandomVariables);
	if(x!=0){ delete x; x=0; }
	x = new Vector(numberOfRandomVariables);
	if(u!=0){ delete u; u=0; }
	u = new Vector(numberOfRandomVariables);
	if(u_old!=0){ delete u_old; u_old=0; }
	u_old = new Vector(numberOfRandomVariables);
	if(searchDirection!=0){ delete searchDirection; searchDirection=0; }
	searchDirection = new Vector(numberOfRandomVariables);
	if(uSecondLast!=0){ delete uSecondLast; uSecondLast=0; }
	uSecondLast = new Vector(numberOfRandomVariables);
	if(uThirdLast!=0){ delete uThirdLast; uThirdLast=0; }
	uThirdLast = new Vector(numberOfRandomVariables);
//	Vector uNew(numberOfRandomVariables);
	if(alpha!=0){ delete alpha; alpha=0; }
	alpha = new Vector(numberOfRandomVariables);
	if(alpha_old!=0){ delete alpha_old; alpha_old=0; }
	alpha_old = new Vector(numberOfRandomVariables);
	if(gamma!=0){ delete gamma; gamma=0; }
	gamma = new Vector(numberOfRandomVariables);
	if(alphaSecondLast!=0){ delete alphaSecondLast; alphaSecondLast=0; }
	alphaSecondLast = new Vector(numberOfRandomVariables);
	if(alphaThirdLast!=0){ delete alphaThirdLast; alphaThirdLast=0; }
	alphaThirdLast = new Vector(numberOfRandomVariables);
	gFunctionValue = 1.0;
	gFunctionValue_old = 1.0;
	if(gradientOfgFunction!=0){ delete gradientOfgFunction; gradientOfgFunction=0; }
	gradientOfgFunction = new Vector(numberOfRandomVariables);
	if(gradientInStandardNormalSpace!=0){ delete gradientInStandardNormalSpace; gradientInStandardNormalSpace=0; }
	gradientInStandardNormalSpace = new Vector(numberOfRandomVariables);
	if(gradientInStandardNormalSpace_old!=0){ delete gradientInStandardNormalSpace_old; gradientInStandardNormalSpace_old=0; }
	gradientInStandardNormalSpace_old = new Vector(numberOfRandomVariables);
	normOfGradient =0;
	stepSize=0.0;
	if(jacobian_x_u!=0){ delete jacobian_x_u; jacobian_x_u=0; }
	jacobian_x_u = new Matrix(numberOfRandomVariables,numberOfRandomVariables);

	if(x_min!=0){ delete x_min; x_min=0; }
	x_min = new Vector(numberOfRandomVariables);
	if(u_min!=0){ delete u_min; u_min=0; }
	u_min = new Vector(numberOfRandomVariables);
	if(alpha_min!=0){ delete alpha_min; alpha_min=0; }
	alpha_min = new Vector(numberOfRandomVariables);
	if(gradientInStandardNormalSpace_min!=0){ delete gradientInStandardNormalSpace_min; gradientInStandardNormalSpace_min=0; }
	gradientInStandardNormalSpace_min = new Vector(numberOfRandomVariables);
	if(gradientOfgFunction_min!=0){ delete gradientOfgFunction_min; gradientOfgFunction_min=0; }
	gradientOfgFunction_min = new Vector(numberOfRandomVariables);
	if(jacobian_x_u_min!=0){ delete jacobian_x_u_min; jacobian_x_u_min=0; }
	jacobian_x_u_min = new Matrix(numberOfRandomVariables,numberOfRandomVariables);

	
	int evaluationInStepSize = 0;
	int result;

	theGFunEvaluator->initializeNumberOfEvaluations();
	theGradGEvaluator->initializeNumberOfEvaluations();

 	fdf=theGradGEvaluator->getfinitedifference();

	// Prepare output file to store the search points
	ofstream outputFile2( fileNamePrint, ios::out );
	if (printFlag == 0) {
		outputFile2 << "The user has not specified to store any search points." << endln;
		outputFile2 << "This is just a dummy file. " << endln;
	}

	if (startAtOrigin)
	  u->Zero();
	else {
	  //theReliabilityDomain->getStartPoint(*x);
	  for (int i = 0; i < numberOfRandomVariables; i++) {
	    RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
	    (*x)(i) = theRV->getStartValue();
	  }

	  result = theProbabilityTransformation->transform_x_to_u(*u);
	  if (result < 0) { 
	    errorMessage_xtou(); delxinit();  return -1; 
	  }
	}

	convergenceAchieved=false;

	/// work area ///
	Vector** uSave;
	uSave= new Vector*[maxNumberOfIterations+1];
	for ( int i=0; i<=maxNumberOfIterations; i++)	{
		uSave[i] =new Vector(numberOfRandomVariables);
	}
	Matrix product(maxNumberOfIterations+1,maxNumberOfIterations+1);
	Matrix distance(maxNumberOfIterations+1,maxNumberOfIterations+1);
	Vector ulengthSave(maxNumberOfIterations+1);
	Vector gValue(maxNumberOfIterations+1);
	Vector vcheck1(maxNumberOfIterations+1);
	Vector vcheck2(maxNumberOfIterations+1);
	Vector gdelta(maxNumberOfIterations+1);
	Vector chkSave(maxNumberOfIterations+1);
	Vector beta1Save(maxNumberOfIterations+1);
	Vector beta2Save(maxNumberOfIterations+1);
	Vector beta3Save(maxNumberOfIterations+1);
	Vector duch1Save(maxNumberOfIterations+1);
	Vector duch2Save(maxNumberOfIterations+1);
	Vector duch3Save(maxNumberOfIterations+1);
	(*uSave[0])=(*u);
	ulengthSave(0)=(*u).Norm();
	/// work area ///

	(*u_old)=(*u);
	(*uSecondLast)=(*u);
	(*uThirdLast)=(*u);
	(*alphaThirdLast)=(*uThirdLast);
	(*alphaThirdLast).Normalize();

	ampchk=1000.0;
	chk_min=10000.0;

	// Loop to find design point
	iter = 1;
	resFlag=-1;
	while ( iter <= maxNumberOfIterations )
	{
		// Transform from u to x space
		//result = theProbabilityTransformation->set_u(*u);
		//if (result < 0) { errorMessage_setu(); delxinit(); return -1; }
		//result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		//if (result < 0) { errorMessage_utox(); delxinit(); return -1; }
		//(*x) = theProbabilityTransformation->get_x();
		//(*jacobian_x_u) = theProbabilityTransformation->getJacobian_x_u();
		result = theProbabilityTransformation->transform_u_to_x(*u, *x);
		if (result < 0) { 
			errorMessage_xtou(); delxinit();  return -1; 
		}
		result = theProbabilityTransformation->getJacobian_x_to_u(*jacobian_x_u);

		// Possibly print the point to output file
		int iii;
		if (printFlag != 0) {
			if (printFlag == 1) {
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<(*x).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*x)(iii)<<endln;
				}
			}
			else if (printFlag == 2) {
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<(*u).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*u)(iii)<<endln;
				}
			}
		}
		// Evaluate limit-state function unless it has been done in 
		// a trial step by the "stepSizeAlgorithm"
		if (evaluationInStepSize == 0) {
			/////// modified by K Fujimura//////////////////
			if(!fdf) theGFunEvaluator->activateSensitivty();
			else theGFunEvaluator->inactivateSensitivty();
			/////// modified by K Fujimura//////////////////
			result = theGFunEvaluator->runAnalysis();
			if (result < 0) { errorMessage_gfun(); delxinit(); return -1; }
			gFunctionValue_old = gFunctionValue;
            gFunctionValue = theGFunEvaluator->evaluateExpression();
		}
		// Set scale parameter
		if (iter == 1)	{
			Gfirst = gFunctionValue;
			opserr << " Limit-state function value at start point, g=" << gFunctionValue << endln;
			opserr << " STEP #0: ";
			theReliabilityConvergenceCheck->setScaleValue(gFunctionValue);
		}

		/// work ///
		gValue(iter-1)=gFunctionValue;
		/// work ///

		// Gradient in original space
		result = theGradGEvaluator->computeGradient(gFunctionValue);
		if (result < 0) { errorMessage_compGradg(); delxinit(); return -1; }
		(*gradientOfgFunction) = theGradGEvaluator->getGradient();

		int noutput=(*gradientOfgFunction).Size();
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
        output<<" step "<<	iter << "\n";
		for(int ijk=0; ijk<noutput; ijk++){
			output << setw(30) << setprecision(10) <<(*gradientOfgFunction)(ijk);
			output << "\n";
		}
		output.flush();

		// Check if all components of the vector is zero
	  	zeroFlag = 0;
		for (j=0; j<(*gradientOfgFunction).Size(); j++) {
			if ((*gradientOfgFunction)[j] != 0.0) {
				zeroFlag = 1;
			}
			if ( zeroFlag == 1 ) break; 
		}
		if (zeroFlag == 0) {errorMessage_checkGradg();delxinit(); return -1;}
		// Gradient in standard normal space
		(*gradientInStandardNormalSpace_old) = (*gradientInStandardNormalSpace);
		(*gradientInStandardNormalSpace) = (*jacobian_x_u)^(*gradientOfgFunction);
		// Compute the norm of the gradient in standard normal space
		normOfGradient = (*gradientInStandardNormalSpace).Norm();
		// Check that the norm is not zero
		if (normOfGradient == 0.0) {errorMessage_zeroGradg();delxinit(); return -1;}
		// Compute alpha-vector
		(*alpha) = (*gradientInStandardNormalSpace) *  ( (-1.0) / normOfGradient );
		// Check convergence
		chkfunc=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
		reschk = theReliabilityConvergenceCheck->check(*u,gFunctionValue,*gradientInStandardNormalSpace);
		chk1=theReliabilityConvergenceCheck->getCheck1();
		chk2=theReliabilityConvergenceCheck->getCheck2();

		/// work ///
		gdelta(iter-1)=fabs(gFunctionValue/normOfGradient);
		vcheck1(iter-1)=chk1;
		vcheck2(iter-1)=chk2;
		chkSave(iter-1)=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
		/// work ///
	
	 	if(iter==1){
			check1_init=chk1;
			check2_init=chk2;
		}

		if (reschk > 0)  {
			// Inform the user of the happy news!
		  	opserr << "Design point found!" << endln;
			check1_conv=chk1;
  			check2_conv=chk2;
			convergenceAchieved=true;
			resFlag=1;
		}else if(reschk ==-1){
			// only direction is in convergence
			// perform line search //
  			opserr << " criterion 2 is satisfied try line search "<< endln;
			int ifound=this->lineSearch();
			if(ifound==2){
				/// work ///
				chkfunc=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
				gdelta(iter-1)=fabs(gFunctionValue/normOfGradient);
				vcheck1(iter-1)=chk1;
				vcheck2(iter-1)=chk2;
				chkSave(iter-1)=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
				/// work ///
				if(convergenceAchieved)resFlag=1;
			}
		}

		if(convergenceAchieved) break;

		// check the change of location if iter 
		if(iter>=4){
			int istop=1;
			double ulenght0=(*uThirdLast).Norm();
			double beta0=(*alphaThirdLast)^(*uThirdLast);
			double beta1=(*alphaSecondLast)^(*uSecondLast);
			double beta2=(*alpha_old)^(*u_old);
			double beta3=(*alpha)^(*u);
			double beta_base=0.0;
			beta_base=beta3;
			double btchange1=fabs(beta0-beta_base)/fabs(beta_base);
			double btchange2=fabs(beta1-beta_base)/fabs(beta_base);
			double btchange3=fabs(beta2-beta_base)/fabs(beta_base);
			beta1Save(iter-1)=btchange1;
			beta2Save(iter-1)=btchange2;
			beta3Save(iter-1)=btchange3;
			if(btchange1>0.01) istop=0;
			if(btchange2>0.01) istop=0;
		 	if(btchange3>0.01) istop=0;
			double ulength=(*u).Norm();
			Vector udif=(*uThirdLast)-(*u);
			double duch3=udif.Norm()/ulength; 
			udif=(*uSecondLast)-(*u);
			double duch2=udif.Norm()/ulength; 
			udif=(*u_old)-(*u);
			double duch1=udif.Norm()/ulength; 
			duch1Save(iter-1)=duch1;
			duch2Save(iter-1)=duch2;
			duch3Save(iter-1)=duch3;
			if(duch1>0.01) istop=0;
			if(duch2>0.01) istop=0;
			if(duch3>0.01) istop=0;
			if(chkfunc<=chk_min){
				(*u_min)=(*u);
				(*x_min)=(*x);
				(*alpha_min)=(*alpha);
				(*gradientInStandardNormalSpace_min)=(*gradientInStandardNormalSpace);
				(*gradientOfgFunction_min)=(*gradientOfgFunction);
				(*jacobian_x_u_min)=(*jacobian_x_u);
				chk_min=chkfunc;
				chk1_min=chk1;
				chk2_min=chk2;
				gFunctionValue_min=gFunctionValue;
				normOfGradient_min=normOfGradient;
			}
			if(istop==1) {
				if(reschk!=-2){
  					opserr << " line search "<< endln;
					int ifound=this->lineSearch();
					if(ifound==2){
						/// work ///
						gdelta(iter-1)=fabs(gFunctionValue/normOfGradient);
						vcheck1(iter-1)=chk1;
						vcheck2(iter-1)=chk2;
						chkSave(iter-1)=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
						/// work ///
						if(convergenceAchieved)	resFlag=1;
					}
				}
				convergenceAchieved=true;
				if(resFlag!=1) {
					resFlag=2;
					check1_conv=chk1;
					check2_conv=chk2;
					chkfunc=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
					if(chkfunc<=chk_min){
						(*u_min)=(*u);
						(*x_min)=(*x);
						(*alpha_min)=(*alpha);
						(*gradientInStandardNormalSpace_min)=(*gradientInStandardNormalSpace);
						(*gradientOfgFunction_min)=(*gradientOfgFunction);
						(*jacobian_x_u_min)=(*jacobian_x_u);
						chk_min=chkfunc;
						chk1_min=chk1;
						chk2_min=chk2;
						gFunctionValue_min=gFunctionValue;
						normOfGradient_min=normOfGradient;
					}
				}
			}
		}	

		if(convergenceAchieved) break;
	
		// Store 'u' and 'alpha' at the second last iteration point
		if(iter>=4)	{
			(*uThirdLast)=(*uSecondLast);
			(*alphaThirdLast)=(*alphaSecondLast);
		}
		if(iter>=3)	{
			(*uSecondLast)=(*u_old);
			(*alphaSecondLast)=(*alpha_old);
		}
		(*u_old) = (*u);
		(*alpha_old) = (*alpha);

		// Let user know that we have to take a new step
		opserr << " STEP #" << iter <<": ";

		// Update Hessian approximation, if any
		if (  (theHessianEvaluator != 0) && (iter != 1)  ) {
			//theHessianApproximation->updateHessianApproximation(*u_old,gFunctionValue_old,
			//		*gradientInStandardNormalSpace_old,stepSize,*searchDirection,
            //        gFunctionValue,*gradientInStandardNormalSpace);
            theHessianEvaluator->computeHessian();
		}


		// Determine search direction
		result = theSearchDirection->computeSearchDirection(iter,*u, gFunctionValue, 
                                            *gradientInStandardNormalSpace );
		if (result < 0) {
			opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute search direction. " << endln;
			if(xinit!=0){
				delete xinit;
				xinit=0;
			}
			return -1;
		}
		(*searchDirection) = theSearchDirection->getSearchDirection();


		// Determine step size
		result = theStepSizeRule->computeStepSize(
			*u, *gradientInStandardNormalSpace, gFunctionValue, *searchDirection, iter, reschk);
		if (result < 0) {  // (something went wrong)
			opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute step size. " << endln;
			if(xinit!=0){
				delete xinit;
				xinit=0;
			}
			return -1;
		}
		else if (result == 0) {  // (nothing was evaluated in step size)
			evaluationInStepSize = 0;
		}
		else if (result == 1) {  // (the gfun was evaluated)
			evaluationInStepSize = 0;
////////////////////////////////////////////////////////////////////////
//////////// Kazuya Fujimura/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//			evaluationInStepSize = 1;
			gFunctionValue_old = gFunctionValue;
            // KRM -- step size no longer storing gFun value, get it from algorithm itself
			//gFunctionValue = theStepSizeRule->getGFunValue();
		}
		stepSize = theStepSizeRule->getStepSize();
		int numReduction = theStepSizeRule->getNumReductions();
		outputTMP<<"stepSize="<<stepSize<<"\n";
		outputTMP<<"numReduction="<<numReduction<<"\n";
		outputTMP<<"\n";
		outputTMP.flush();


		// Determine new iteration point (take the step)
		(*u) = (*u_old) + (*searchDirection)*stepSize;

		(*uSave[iter])=(*u);
		ulengthSave(iter)=(*uSave[iter]).Norm();
		// Increment the loop parameter
		double ulengthi=(*uSave[iter]).Norm();
		for(int jjj=0;jjj<=iter-1;jjj++){
			double ulengthj=(*uSave[jjj]).Norm();
			product(iter,jjj)=1.0-((*uSave[iter])^(*uSave[jjj]))/(ulengthi*ulengthj);
			Vector diff=(*uSave[iter])-(*uSave[jjj]);
			distance(iter,jjj)=diff.Norm();
		}

		iter++;
	}
	// iteration loop end //
	
	if(resFlag!=1&&resFlag!=2) {
		// NON CONVERGENCE CASE //
		resFlag=3;
		if(reschk!=-2){
			// line Search //
			int ifound=this->lineSearch();
			gdelta(iter-2)=fabs(gFunctionValue/normOfGradient);
			vcheck1(iter-2)=chk1;
			vcheck2(iter-2)=chk2;
			chkSave(iter-2)=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
			if(convergenceAchieved)	resFlag=1;
		}
		if(resFlag!=1){
			convergenceAchieved=true;
			check1_conv=chk1;
			check2_conv=chk2;
			chkfunc=(*u).Norm()+ampchk*fabs(gFunctionValue/normOfGradient);
			if(chkfunc<=chk_min){
				(*u_min)=(*u);
				(*x_min)=(*x);
				(*alpha_min)=(*alpha);
				(*gradientInStandardNormalSpace_min)=(*gradientInStandardNormalSpace);
				(*gradientOfgFunction_min)=(*gradientOfgFunction);
				(*jacobian_x_u_min)=(*jacobian_x_u);
				chk_min=chkfunc;
				chk1_min=chk1;
				chk2_min=chk2;
				gFunctionValue_min=gFunctionValue;
				normOfGradient_min=normOfGradient;
			}
		}
	}

	if(resFlag!=1){
	// non convergence //
		convergenceAchieved=true;
		(*u)=(*u_min);
		(*x)=(*x_min);
		(*alpha)=(*alpha_min);
		(*gradientInStandardNormalSpace)=(*gradientInStandardNormalSpace_min);
		(*gradientOfgFunction)=(*gradientOfgFunction_min);
		(*jacobian_x_u)=(*jacobian_x_u_min);
		chkfunc=chk_min;
		check1_conv=chk1_min;
		check2_conv=chk2_min;
		gFunctionValue=gFunctionValue_min;
		normOfGradient=normOfGradient_min;
		outputTMP<<"Convergence is not achieved\n";
		outputTMP<<"point with chkfunc  "<<chkfunc<<" is selected\n";
		outputTMP<<"\n";
	}

	outputTMP<<"Gfunc and uLength\n";
	outputTMP<<"\n";
	outputTMP<<"         gValue";
	outputTMP<<"        ulength";
	outputTMP<<"         gdelta";
	outputTMP<<"        chkfunc";
	outputTMP<<"           chk1";
	outputTMP<<"           chk2";
	outputTMP<<"            bt1";
	outputTMP<<"            bt2";
	outputTMP<<"            bt3";
	outputTMP<<"            du1";
	outputTMP<<"            du2";
	outputTMP<<"            du3";
	outputTMP<<"\n";
  	for(int iii=0;iii<=iter;iii++){
		outputTMP << setw(15) << setprecision(5) << gValue(iii);
		outputTMP << setw(15) << setprecision(5) << ulengthSave(iii);
		outputTMP << setw(15) << setprecision(5) << gdelta(iii);
		outputTMP << setw(15) << setprecision(5) << chkSave(iii);
		outputTMP << setw(15) << setprecision(5) << vcheck1(iii);
		outputTMP << setw(15) << setprecision(5) << vcheck2(iii);
		outputTMP << setw(15) << setprecision(5) << beta1Save(iii);
		outputTMP << setw(15) << setprecision(5) << beta2Save(iii);
		outputTMP << setw(15) << setprecision(5) << beta3Save(iii);
		outputTMP << setw(15) << setprecision(5) << duch1Save(iii);
		outputTMP << setw(15) << setprecision(5) << duch2Save(iii);
		outputTMP << setw(15) << setprecision(5) << duch3Save(iii);
		outputTMP<<"\n";
	}
	outputTMP<<"\n";
	outputTMP<<"innner product matrix\n";
	outputTMP<<"\n";
	for(int iii=1;iii<=iter;iii++){
		for(int jjj=0;jjj<=iii-1;jjj++){
			outputTMP << setw(15) << setprecision(5) << product(iii,jjj);
		}
		outputTMP<<"\n";
	}
	outputTMP<<"distance matrix\n";
	outputTMP<<"\n";
	for(int iii=1;iii<=iter;iii++){
		for(int jjj=0;jjj<=iii-1;jjj++){
			outputTMP << setw(15) << setprecision(5) << distance(iii,jjj);
		}
		outputTMP<<"\n";
	}
	outputTMP.flush();

	for (int i=0; i<=maxNumberOfIterations; i++)	{
		delete uSave[i]; uSave[i]=0;
	}
  	delete uSave; uSave=0;


	if(convergenceAchieved){
		int iii;
		if (printFlag != 0) {
			if (printFlag == 3) {
				/*
				result = theProbabilityTransformation->set_u(*u);
				if (result < 0) {
					opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
						<< " could not set u in the xu-transformation." << endln;
					if(xinit!=0){
						delete xinit;
						xinit=0;
					}
					return -1;
				}
				result = theProbabilityTransformation->transform_u_to_x();
				if (result < 0) {
					opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
							<< " could not transform from u to x." << endln;
					if(xinit!=0){
						delete xinit;
						xinit=0;
					}
					return -1;
				}
				(*x) = theProbabilityTransformation->get_x();
				*/
				result = theProbabilityTransformation->transform_u_to_x(*u, *x);
				if (result < 0) { 
					errorMessage_xtou(); delxinit();  return -1; 
				}
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<(*x).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*x)(iii)<<endln;
				}
			}
			else if (printFlag == 4) {
				static ofstream outputFile2( fileNamePrint, ios::out );
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<(*u).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*u)(iii)<<endln;
				}
			}
		}
		/*
		// Compute the gamma vector
		MatrixOperations theMatrixOperations(*jacobian_x_u);
		result = theMatrixOperations.computeTranspose();
		if (result < 0) {
			opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute transpose of jacobian matrix. " << endln;
			if(xinit!=0){
				delete xinit;
				xinit=0;
			}
			return -1;
		}
		Matrix transposeOfJacobian_x_u = theMatrixOperations.getTranspose();
		Matrix jacobianProduct = (*jacobian_x_u) * transposeOfJacobian_x_u;
		Matrix D_prime(numberOfRandomVariables,numberOfRandomVariables);
		for (j=0; j<numberOfRandomVariables; j++) {
			D_prime(j,j) = sqrt(jacobianProduct(j,j));
		}
		Matrix jacobian_u_x(numberOfRandomVariables,numberOfRandomVariables);
		theProbabilityTransformation->getJacobian_u_to_x(*u, jacobian_u_x);
		Vector tempProduct = jacobian_u_x ^ (*alpha);
		(*gamma) = D_prime ^ tempProduct;
		*/
		
		// Compute the gamma vector
		//const Matrix &jacobian_u_x = theProbabilityTransformation->getJacobian_u_x();
		// Get Jacobian u-space to x-space
		Matrix Jux(numberOfRandomVariables,numberOfRandomVariables);
		Matrix Jxu(numberOfRandomVariables,numberOfRandomVariables);
		result = theProbabilityTransformation->getJacobian_u_to_x(*u, Jux);
		if (result < 0) {
		  opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			 << " could not transform from u to x." << endln;
		  return -1;
		}

		//Vector tempProduct = jacobian_u_x ^ alpha;
		Vector tempProduct(numberOfRandomVariables);
		//tempProduct.addMatrixTransposeVector(0.0, jacobian_u_x, alpha, 1.0);
		tempProduct.addMatrixTransposeVector(0.0, Jux, *alpha, 1.0);

		// Only diagonal elements of (J_xu*J_xu^T) are used
		for (j = 0; j < numberOfRandomVariables; j++) {
		  double sum = 0.0;
		  double jk;
		  for (int k = 0; k < numberOfRandomVariables; k++) {
		    //jk = jacobian_x_u(j,k);
		    jk = Jxu(j,k);
		    sum += jk*jk;
		  }
		  (*gamma)(j) = sqrt(sum) * tempProduct(j);
		}
		
		Glast = gFunctionValue;
		numberOfEvaluations = theGFunEvaluator->getNumberOfEvaluations();
		numberOfSensAna = theGradGEvaluator->getNumberOfEvaluations();
		if(xinit!=0){
			delete xinit;
			xinit=0;
		}
		// Print the design point to file, if desired
		return resFlag;
	}else{
		// Print a message if max number of iterations was reached
		// (Note: in this case the last trial point was never transformed/printed)
		opserr << "Maximum number of iterations was reached before convergence." << endln;

		if (printFlag != 0) {
			if (printFlag == 3) {
				/*
				result = theProbabilityTransformation->set_u(*u);
				if (result < 0) {
					opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
						<< " could not set u in the xu-transformation." << endln;
					if(xinit!=0){
						delete xinit;
						xinit=0;
					}
					return -1;
				}
				result = theProbabilityTransformation->transform_u_to_x();
				if (result < 0) {
					opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
							<< " could not transform from u to x." << endln;
					if(xinit!=0){
						delete xinit;
						xinit=0;
					}
					return -1;
				}
				(*x) = theProbabilityTransformation->get_x();
				*/
				theProbabilityTransformation->transform_u_to_x(*u, *x);
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (int iii=0; iii<(*x).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*x)(iii)<<endln;
				}
			}
			else if (printFlag == 4) {
				static ofstream outputFile2( fileNamePrint, ios::out );
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (int iii=0; iii<(*u).Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<(*u)(iii)<<endln;
				}
			}
		}




	}
	if(xinit!=0){
		delete xinit;
		xinit=0;
	}

	return -1;
}



const Vector&
NewSearchWithStepSizeAndStepDirection::get_x()
{
	return (*x);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::get_u()
{
	return (*u);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::get_alpha()
{
	return (*alpha);
}
double
NewSearchWithStepSizeAndStepDirection::get_beta()
{
	return (*u)^(*alpha);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::get_gamma()
{
	//return (*gamma)*(1.0/(*gamma).Norm());
	if (gamma->Norm() > 1.0)
		(*gamma) = (*gamma) / gamma->Norm();

	return *gamma;
}

int
NewSearchWithStepSizeAndStepDirection::getNumberOfSteps()
{
	return (iter-1);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::getSecondLast_u()
{
	return (*uSecondLast);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::getSecondLast_alpha()
{
	return (*alphaSecondLast);
}

const Vector&
NewSearchWithStepSizeAndStepDirection::getLastSearchDirection()
{
	return (*searchDirection);
}

double
NewSearchWithStepSizeAndStepDirection::getFirstGFunValue()
{
	return Gfirst;
}

double
NewSearchWithStepSizeAndStepDirection::getLastGFunValue()
{
	return Glast;
}


const Vector&
NewSearchWithStepSizeAndStepDirection::getGradientInStandardNormalSpace()
{
	return (*gradientInStandardNormalSpace);
}



int
NewSearchWithStepSizeAndStepDirection::getNumberOfEvaluations()
{
	return numberOfEvaluations;
}
int
NewSearchWithStepSizeAndStepDirection::getNumberOfSensAna()
{
	return numberOfSensAna;
}

void
NewSearchWithStepSizeAndStepDirection::set_x(Vector& xin)
{
	if( xinit !=0){
		delete xinit;
		xinit=0;
	}
	xinit = new Vector(xin);
}
void
NewSearchWithStepSizeAndStepDirection::set_u(Vector& uin)
{
	int result;
	/*
	result = theProbabilityTransformation->set_u(uin);
	if (result < 0) {
		opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set u in the xu-transformation." << endln;
		exit(-1);
	}
	result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
	if (result < 0) {
		opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set u in the xu-transformation." << endln;
		exit(-1);
	}
	(*x) = theProbabilityTransformation->get_x();
	*/
	theProbabilityTransformation->transform_u_to_x(uin, *x);
	if( xinit !=0){
		delete xinit;
		xinit=0;
	}
	xinit = new Vector(*x);
}
void
NewSearchWithStepSizeAndStepDirection::errorMmessage_setx()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set x in the xu-transformation." << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_xtou()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not transform from x to u." << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_setu()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
		<< " could not set u in the xu-transformation." << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_utox()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not transform from u to x and compute Jacobian." << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_gfun()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_evalg()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not tokenize limit-state function. " << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_compGradg()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not compute gradients of the limit-state function. " << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_checkGradg()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " all components of the gradient vector is zero. " << endln;
}
void
NewSearchWithStepSizeAndStepDirection::errorMessage_zeroGradg()
{
	opserr << "NewSearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
		<< " the norm of the gradient is zero. " << endln;
}
void 
NewSearchWithStepSizeAndStepDirection::delxinit()
{
	if(xinit!=0){delete xinit;xinit=0;}
}
double
NewSearchWithStepSizeAndStepDirection::get_check1_init()
{
	return check1_init;
}
double
NewSearchWithStepSizeAndStepDirection::get_check2_init()
{
	return check2_init;
}
double
NewSearchWithStepSizeAndStepDirection::get_check1_conv()
{
	return check1_conv;
}
double
NewSearchWithStepSizeAndStepDirection::get_check2_conv()
{
	return check2_conv;
}
int
NewSearchWithStepSizeAndStepDirection::lineSearch()
{
	theGFunEvaluator->inactivateSensitivty();
	// variables 
	double gpositive = 0.0; 
	double gnegative = 0.0;
	double amppositive = 0.0;
	double ampnegative = 0.0;
	double damp = 0.1;
	Vector uBase = (*u);
	Vector uDelta = (*u)*damp;
	Vector uTest = (*u);
	Vector xTest = (*x);
	double positive = 0.0;
	double gTest = 0.0;
	double ampTest = 0.0;
	if( gFunctionValue > 0 ){
		gpositive = gFunctionValue;
		positive = 1.0;
		opserr << " amppositive "<<amppositive<<endln;
		opserr << " gpositive "<<gpositive<<endln;
	}else{
		gnegative = gFunctionValue;
		positive = -1.0;
		opserr << " ampnegative "<<ampnegative<<endln;
		opserr << " gnegative "<<gnegative<<endln;
	}
	int ifound = 0;
	int result;
	for( int ijk=0; ijk<20; ijk++ ){
		ampTest = positive*damp*(double)(ijk+1);
		uTest=uBase+uDelta*ampTest;
		/*
		result = theProbabilityTransformation->set_u(uTest);
		if (result < 0) { errorMessage_setu(); delxinit(); return -1; }
		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) { errorMessage_utox(); delxinit(); return -1; }
		xTest = theProbabilityTransformation->get_x();
		*/
		theProbabilityTransformation->transform_u_to_x(uTest, xTest);
		result = theGFunEvaluator->runAnalysis();
		if (result < 0) { errorMessage_gfun(); delxinit(); return -1; }
		gTest = theGFunEvaluator->evaluateExpression();
		opserr << " ampTest "<<ampTest<<endln;
		opserr << " gTest "<<gTest<<endln;
		result = theReliabilityConvergenceCheck->checkG(gTest);
		if(result > 0 ){
			ifound = 2;
  			opserr << " lineSearch point found "<<endln;
  			opserr << " ampTest "<<ampTest<<endln;
			(*u)=uTest;
			(*x)=xTest;
		}else{
			if(positive>0.0){
				// look for negative;
				if(gTest < 0.0) {
					ifound = 1;	ampnegative = ampTest; gnegative=gTest;
	  				opserr << " the other end  found "<<endln;
	  				opserr << " ampnegative "<<ampnegative<<endln;
					opserr << " gnegative "<<gnegative<<endln;
				}
			}else{
				// look for positive;
				if(gTest > 0.0) {
					ifound = 1;	amppositive = ampTest; gpositive=gTest;
		  			opserr << " the other end  found "<<endln;
		  			opserr << " amppositive "<<amppositive<<endln;
	  				opserr << " gpositive "<<gpositive<<endln;
				}
			}
		}
		if( ifound !=0) break;
	}

	if( ifound ==1 ){
		// perform line search by secant method //
		// if point is found change ifound from 1 to 2 //
		opserr << " Start Line Search "<<endln;
		opserr << " ampnegative "<<ampnegative<<endln;
		opserr << " gnegative "<<gnegative<<endln;
		opserr << " amppositive "<<amppositive<<endln;
		opserr << " gpositive "<<gpositive<<endln;
		int nlineSearch=0;
		do {
			ampTest = (gpositive*ampnegative-gnegative*amppositive)/(gpositive-gnegative);
			opserr << " ampTest "<<ampTest<<endln;
			uTest=uBase+uDelta*ampTest;
			//result = theProbabilityTransformation->set_u(uTest);
			//if (result < 0) { errorMessage_setu(); delxinit(); return -1; }
			//result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
			//if (result < 0) { errorMessage_utox(); delxinit(); return -1; }
			//xTest = theProbabilityTransformation->get_x();
			result = theProbabilityTransformation->transform_u_to_x(uTest, xTest);
			result = theGFunEvaluator->runAnalysis();
			if (result < 0) { errorMessage_gfun(); delxinit(); return -1; }
			gTest = theGFunEvaluator->evaluateExpression();
			result = theReliabilityConvergenceCheck->checkG(gTest);
			opserr << " gTest "<<gTest<<endln;
			nlineSearch++;
			if( result < 0 ){
				if( gTest > 0 ){
					gpositive = gTest;
					amppositive = ampTest;
				}else{
					gnegative = gTest;
					ampnegative = ampTest;
				}
			}else{
				ifound=2;
				(*u)=uTest;
			}
			if(nlineSearch>=20) break;
		}while(ifound!=2);

		if(ifound==2){
			(*u)=uTest;
			(*x)=xTest;
	  		opserr << " lineSearch point found "<<endln;
	  		opserr << " ampTest "<<ampTest<<endln;
		}else{
			ifound=0;
	  		opserr << " lineSearch point could not be found "<<endln;
	  		opserr << " continue current search  "<<ampTest<<endln;
		}
	}
	if( ifound != 2 ){
		return -4;
	}else{
		// lineserch point is found! check gradient now 
	  	opserr << " now obtain Gradient  "<<endln;
		//(*jacobian_x_u) = theProbabilityTransformation->getJacobian_x_u();
		result = theProbabilityTransformation->getJacobian_x_to_u(*jacobian_x_u);
		if(!fdf) theGFunEvaluator->activateSensitivty();
		else theGFunEvaluator->inactivateSensitivty();
		result = theGFunEvaluator->runAnalysis();
		if (result < 0) { errorMessage_gfun(); delxinit(); return -1; }
		gFunctionValue = theGFunEvaluator->evaluateExpression();
		// Gradient in original space
		result = theGradGEvaluator->computeGradient(gFunctionValue);
		if (result < 0) { errorMessage_compGradg(); delxinit(); return -1; }
		(*gradientOfgFunction) = theGradGEvaluator->getGradient();
//		int noutput=gradientOfgFunction.Size();
//		output.setf(ios::right);
//		output.setf(ios::scientific, ios::floatfield);
//		for(int ijk=0; ijk<noutput; ijk++){
//		output << setw(30) << setprecision(10) <<gradientOfgFunction(ijk);
//		output << "\n";
//		}
//		output.flush();
		// Check if all components of the vector is zero
		int zeroFlag = 0;
		for (int j=0; j<(*gradientOfgFunction).Size(); j++) {
			if ((*gradientOfgFunction)[j] != 0.0) {
				zeroFlag = 1;
			}
			if ( zeroFlag == 1 ) break; 
		}
		if (zeroFlag == 0) {errorMessage_checkGradg();delxinit(); return -1;}
		// Gradient in standard normal space
		(*gradientInStandardNormalSpace) = (*jacobian_x_u) ^ (*gradientOfgFunction);
		// Compute the norm of the gradient in standard normal space
		normOfGradient = (*gradientInStandardNormalSpace).Norm();
		// Check that the norm is not zero
		if (normOfGradient == 0.0) {errorMessage_zeroGradg();delxinit(); return -1;}
		// Compute alpha-vector
		(*alpha) = (*gradientInStandardNormalSpace) *  ( (-1.0) / normOfGradient );
		// Check convergence
		reschk = theReliabilityConvergenceCheck->check(*u,gFunctionValue,*gradientInStandardNormalSpace);
		chk1=theReliabilityConvergenceCheck->getCheck1();
		chk2=theReliabilityConvergenceCheck->getCheck2();
		if (reschk > 0)  {
			check1_conv=chk1;
			check2_conv=chk2;
			// Inform the user of the happy news!
			opserr << "Design point found!" << endln;
			convergenceAchieved=true;
			resFlag=0;
		}
		return ifound;
	}
}


int 
NewSearchWithStepSizeAndStepDirection::setStartPt(Vector *s)
{
  opserr << "NewSearchWithStepSizeAndStepDirection::setStartPt(Vector *s)  - not implemented yet\n";
  return 0;
}
