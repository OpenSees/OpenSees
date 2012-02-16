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
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//

#include <stdio.h> 
#include <math.h>
#include <OrthogonalPlaneSamplingAnalysis.h>
#include <SafeGuardedZeroFindingAlgorithm.h>

 
#include <RandomNumberGenerator.h>
//#include <RandomVariable.h>
#include <NormalRV.h>
#include <fstream>
#include <iostream>

using std::ifstream;
using std::ios;



OrthogonalPlaneSamplingAnalysis::OrthogonalPlaneSamplingAnalysis(   Tcl_Interp *interp,
		                ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						FunctionEvaluator *passedGFunEvaluator,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int passedMaxNumOfIterations,
						double passedTargetCOV,
						int pPrintFlag,
						TCL_Char *passedFileName,
						Vector * pDesignPoint,
						int pAnalysisTypeTag, 
						int zeroFindingType,
						double pFuncTol,
						double pVarTol,
						int pMaxIter,
						double pLittleDt){
	
	numOfGFunEvaluations=0;

	this->theInterp=interp;
	this->theReliabilityDomain=passedReliabilityDomain;
	this->theProbabilityTransformation=passedProbabilityTransformation;
	this->theGFunEvaluator=passedGFunEvaluator;
	this->theRandomNumberGenerator=passedRandomNumberGenerator;
	this->numOfSimulations=passedNumberOfSimulations;
	this->maxNumOfIterations=passedMaxNumOfIterations;
	this->targetCOV=passedTargetCOV;

	this->printFlag=pPrintFlag;
	this->fileName = new char[256];
	strcpy(fileName,passedFileName);

	this->analysisTypeTag=pAnalysisTypeTag;
	this->theDesignPoint=pDesignPoint;
	if (zeroFindingType ==1) { //safeGuardedZeroFinding..
		theAlgorithm = new SafeGuardedZeroFindingAlgorithm(this);
	}
	else{
		opserr<<"unknown zerofinding algorithm type.."<<endln;
		exit(-1);
	}

	if ((analysisTypeTag ==1)||(analysisTypeTag ==2)){
		double * tmp = new double[100];
		double * tmp1 = new double[100];
		theAlgorithm->setX1Pointer(tmp);
		theAlgorithm->setG1Pointer(tmp1);	
		theAlgorithm->set_ii_1(0);
	}
	if (analysisTypeTag ==2){
		double * tmp2 = new double[100];
		double * tmp3 = new double[100];
		theAlgorithm->setX2Pointer(tmp2);
		theAlgorithm->setG2Pointer(tmp3);			
		theAlgorithm->set_ii_2(0);

	}

	// --- need to be passed..
		littleDt = pLittleDt;



	probability=0;
	cov=0;

	uPrime = 0;

	funcTol = pFuncTol;
	varTol = pVarTol;
	maxIter =pMaxIter;
	theDesignPointInUSpace =0;
//	seed=1;



};


OrthogonalPlaneSamplingAnalysis::~OrthogonalPlaneSamplingAnalysis()
{
	if (theAlgorithm !=0)
		delete theAlgorithm;
	if (fileName !=0)
		delete fileName;
	if (uPrime !=0)
		delete uPrime;
	if (theDesignPointInUSpace !=0)
		delete 	theDesignPointInUSpace;


}


int OrthogonalPlaneSamplingAnalysis::analyze(void)
{
	opserr << "OrthogonalPlaneSamplingAnalysis is running ... " << endln;


	ofstream resultsOutputFile( fileName, ios::out | ios::app );
	if ((printFlag ==1)||(printFlag ==2)) resultsOutputFile<< "Failure Probability" <<","<<"c.o.v."<<endln; 

     
	//--- restart option -- 2007 Feb -----------------

///////////
	cov = 999.0;
	double q,qk,qk2,q_square; // probability index. refer Koo
	q=0.0;  // sum of qk
	
	q_square =0.0;  // sum of qk^2
	int k=1;
	seed =1;
	probability =0.0;
////////////




 	if (printFlag ==2) {
	
		// check whether the restart file '_restart.tmp' exist, (this file is wrote by openSees only, not by user)
		 //                      if yes, read data;{ success reading: set values above; otherwise: do noting} 
		 //		                 if no, create file '_restart.tmp'
	   //	data format:    Pf              cov           k             seed           numOfGFunEvaluations
		//  save into data  Pf_restart   cov_restart     k_restart   seed_restart      numOfGFunEvaluations       
		// Possible read data from file if this is a restart simulation
			
		    ifstream inputFile( "_restart.tmp", ios::in );
			if (!inputFile) {
            // file doesn't exist; so it's safe to create it and write to it
    
			}
			else { 
				inputFile >> probability;
				inputFile >> cov;
				inputFile >> k;
				inputFile >> seed;
				inputFile >> numOfGFunEvaluations;

				q=probability*k;
				q_square = (cov*probability)*(cov*probability)*k*(k-1)+ 1/k*q*q;       ;

				k++;
				inputFile.close();
			}
	
	}
    // ------------------------------------------------









	ofstream resultsOutputFile1( "error_points.out", ios::out );
	resultsOutputFile1<< "x:" <<","<<"F(x)"<<endln;


//	if(printFlag ==5 && this->analysisTypeTag ==2) {
		
		char theSecondFileName[50]="SecondSurface_";
		strcat(theSecondFileName, fileName);
		ofstream resultsOutputFile3( theSecondFileName, ios::out );
//	}

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	
	Vector u(numRV);
	Vector x(numRV);
	
	LimitStateFunction *theLimitStateFunction = 0;
 
	double gFunctionValue;


	int result;
	Vector direction(numRV);
	if (uPrime ==0)  	uPrime = new Vector(numRV);
	if (theDesignPointInUSpace ==0)  	
		theDesignPointInUSpace = new Vector(numRV);
	uPrime->Zero();
	Vector projectionPoint(numRV);
	Vector newPoint(numRV);
	double a,b, Fa, Fb;
	double zeroPoint; // defined as distance  from u'_n to the limit state surface along direction..
	double zeroPoint2;

	if (theDesignPoint ==0) {opserr<<"designpoint does not exist!"<<endln; exit(-1);}
	// note designPoint is in physical space ..

	/*
	result = theProbabilityTransformation->set_x(*theDesignPoint);
	if (result < 0) {
		opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
			<< " set the x-vector for xu-transformation. " << endln;
		return -1;
	}
		
	result = theProbabilityTransformation->transform_x_to_u();
	if (result < 0) {
		opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
			<< " transform x to u. " << endln;
		return -1;
	}
	theDesignPointInUSpace->addVector(0.0, theProbabilityTransformation->get_u(),1.0);
	*/

	result = theProbabilityTransformation->transform_x_to_u(*theDesignPointInUSpace);
	if (result < 0) {
	  opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
		 << " transform x to u. " << endln;
	  return -1;
	}

	//opserr.setPrecision(16);
	//opserr<<"theDesignPoint in x space is:"<<*theDesignPoint<<endln;
	//opserr<<"theDesignPoint in u space is:"<<*theDesignPointInUSpace<<endln;


	beta = theDesignPointInUSpace->Norm();  // reliability index
	


	direction.addVector(0.0, *theDesignPointInUSpace ,1.0/beta);  // sampling direction in u space

	bool needCallingZeroFinding;  
	bool FEconvergence;
//	bool contribution;
	char myString[150];


    double gFuncAtProjectionPoint;
	

    static NormalRV aStdNormRV(1,0.0,1.0);


	if (numOfSimulations<2) {
		opserr<<"warning: too less number of simulations, can not do simulation, adjust simulation number to 2.."<<endln;
		numOfSimulations =2;
	}


	bool isFirstSimulation = true;

	while( (k<=numOfSimulations && cov>targetCOV || k<=2) ) {



		// Keep the user posted
		if (printFlag >0) {
			opserr << "Sample #" << k << ":" << endln;
		}




		theAlgorithm->setFEconvergence(true);
		needCallingZeroFinding = true;
		FEconvergence  =  true;
		contribution = true;


		theAlgorithm->set_ii_1(0);  // always use ii_1, x_1, G_1...




		// ---- step 1, simulate independent U_n ~ N(0,1)----
		
		// Create array of standard normal random numbers
		if (isFirstSimulation) {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV,seed);
		}
		else {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV);
		}
		seed = theRandomNumberGenerator->getSeed();
		if (result < 0) {
			opserr << "OrthogonalSamplingAnalysis::analyze() - could not generate" << endln
				<< " random numbers for simulation." << endln;
			return -1;
		}
		u = theRandomNumberGenerator->getGeneratedNumbers();
		
		//opserr.setPrecision(16);
		//opserr<<"original sampling point u is:"<<u<<endln;		

	// ---- step 2,compute U'_(n-1) by project U_n to the plane: alpha'*U_n =0 ---

		double tmp = u ^ direction;
		
		uPrime->addVector(0.0,u,1.0);
		uPrime->addVector(1.0, direction, -tmp);  // projected point on the orthogonal plane passing origin;

		//opserr<<"uPrime is:"<<*uPrime<<endln;

		projectionPoint.addVector(0.0, *uPrime, 1.0);
		projectionPoint.addVector(1.0, *theDesignPointInUSpace,1.0);   // projection to plane passing through design point. 

		//opserr<<"uProjection is:"<<projectionPoint<<endln;

		// Transform into original space
		/*
		result = theProbabilityTransformation->set_u(projectionPoint);
		if (result < 0) {
			opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
				<< " set the u-vector for xu-transformation. " << endln;
			return -1;
		}
		
		result = theProbabilityTransformation->transform_u_to_x();
		if (result < 0) {
			opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
				<< " transform u to x. " << endln;
			return -1;
		}
		x = theProbabilityTransformation->get_x();
		*/

		result = theProbabilityTransformation->transform_u_to_x(projectionPoint, x);
		if (result < 0) {
			opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
				<< " transform u to x. " << endln;
			return -1;
		}


	//	opserr<<"the "<<k<<"th iteration. x is:" <<x<<endln;

	// ----- step 3, u'_n is solved out by G(U'_(n-1), u'_n)=0 -----
		result = theGFunEvaluator->runAnalysis();
		numOfGFunEvaluations++;

		if (result <0) {
			FEconvergence  =  false; 
			gFunctionValue =-1.0;
			theAlgorithm->setFEconvergence(false);
		}
		else { //if (result > 0) {   // converge
			
			// Get value of limit-state function
			gFunctionValue = theGFunEvaluator->evaluateExpression();

			// ------- save data into zeroFinding algorithm ---
			theAlgorithm->saveXG1(beta, gFunctionValue);

//			FEconvergence  =  true; 
	//		theAlgorithm->setFEconvergence(true);

		}

		gFuncAtProjectionPoint = gFunctionValue;


		if (analysisTypeTag == 1) { //---------------
			
			theAlgorithm->setFunctionType(1);

			if (fabs(gFuncAtProjectionPoint)<= funcTol && FEconvergence) {
				zeroPoint = beta;  
				needCallingZeroFinding = false;
			}  // not need zeroFinding
				
			else if (gFuncAtProjectionPoint>funcTol) { //safe domain, must FEconvergence. this point is a and Fa
				
				a = beta;
				Fa = gFuncAtProjectionPoint;
					
				b = a + 2.0;
				newPoint = *uPrime;
				newPoint.addVector(1.0,direction,b);  // projected point on the orthogonal plane passing origin;
					  
				// Transform into original space
				/*
				result = theProbabilityTransformation->set_u(newPoint);
				if (result < 0) {
					opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
					<< " set the u-vector for xu-transformation. " << endln;
					return -1;
				}
					
				result = theProbabilityTransformation->transform_u_to_x();
				if (result < 0) {
					opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
						<< " transform u to x. " << endln;
					return -1;
				}
				x = theProbabilityTransformation->get_x();
				*/
				result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
				if (result < 0) {
					opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
						<< " transform u to x. " << endln;
					return -1;
				}

				result = theGFunEvaluator->runAnalysis();
				numOfGFunEvaluations++;
				
				if(result<0) {  // diverge
					FEconvergence  =  false; 
					theAlgorithm->setFEconvergence(false);
					gFunctionValue =-1.0;
				}
				else { // converge
					gFunctionValue = theGFunEvaluator->evaluateExpression();
					
					// ------- save data into zeroFinding algorithm ---
					theAlgorithm->saveXG1(b, gFunctionValue);

//					FEconvergence  =  true; 
//					theAlgorithm->setFEconvergence(true);

					
				} //else
				
				Fb = gFunctionValue;
				if (Fb>0) {//no contribution
					contribution = false;
					needCallingZeroFinding =false;
				}
			} //else if (gFuncAtProjectionPoint>funcTol)
				
			else if (gFuncAtProjectionPoint< -funcTol) { //unsafe domain, or diverge case

				b= beta; Fb = gFuncAtProjectionPoint;

				// find a, Fa
				double uPrimeNorm=uPrime->Norm();
				if (beta > uPrimeNorm){
					a= pow(beta*beta- uPrimeNorm*uPrimeNorm, 0.5);
					newPoint = *uPrime;
					newPoint.addVector(1.0,direction,a);  // projected point on the orthogonal plane passing origin;
				
					// Transform into original space
					/*
					result = theProbabilityTransformation->set_u(newPoint);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " set the u-vector for xu-transformation. " << endln;
						return -1;
					}
						
					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
					x = theProbabilityTransformation->get_x();
					*/

					result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
				
					result = theGFunEvaluator->runAnalysis();
					numOfGFunEvaluations++;

					if(result<0) {  // diverge
						opserr<<"something wrong with orthogonalPlaneSamplingAnalysis. need save point "<<endln;
						needCallingZeroFinding =false; 
						//gFunctionValue = -1.0;
						contribution = false;
					// save point ....
						resultsOutputFile1<< "a:" <<a<<","<<"F(a): diverge"<<". uPrime:"<<uPrime<<endln;
					
					}
					else { // converge
						gFunctionValue = theGFunEvaluator->evaluateExpression();
						
						// ------- save data into zeroFinding algorithm ---
						theAlgorithm->saveXG1(a, gFunctionValue);

						Fa = gFunctionValue;
						if (fabs(Fa)<funcTol){
							zeroPoint=a; 
							needCallingZeroFinding =false;
						}
						else if (Fa<-funcTol) {
							opserr<<"something wrong with orthogonalPlaneSamplingAnalysis. need save point "<<endln;
							needCallingZeroFinding =false;
							contribution = false;
							// save point ....
							resultsOutputFile1<< "a:" <<a<<","<<"F(a)"<<Fa<<". uPrime:"<<uPrime<<endln;
						

						}	
					} //else converge

						
				} //if (beta > uPrimeNorm)
				else { //if (beta <= uPrimeNorm)
					a=0;
					newPoint = *uPrime;
     				// Transform into original space
					/*
					  result = theProbabilityTransformation->set_u(newPoint);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " set the u-vector for xu-transformation. " << endln;
						return -1;
					}
					
					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
					x = theProbabilityTransformation->get_x();
					*/
					result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
				
					result = theGFunEvaluator->runAnalysis();
					numOfGFunEvaluations++;
					
					if(result<0) {  // diverge, ignore its contribution
					
						contribution = false;
						needCallingZeroFinding =false;
					}
					else { // converge
						gFunctionValue = theGFunEvaluator->evaluateExpression();

						// ------- save data into zeroFinding algorithm ---
						theAlgorithm->saveXG1(a, gFunctionValue);


						Fa = gFunctionValue;
						if (Fa<0) {   // ignore its contribution
							contribution = false;
							needCallingZeroFinding =false;
				
						}	
					} //else converge
				} // else { //if (beta <= uPrimeNorm)
			} //elseif (gFuncAtProjectionPoint<-funcTol) { //unsafe domain, 

				



 		    if (contribution && needCallingZeroFinding){

					// ======== call zeroFinding ==========
					theAlgorithm->setA(a);   theAlgorithm->setFa(Fa);
					theAlgorithm->setB(b);   theAlgorithm->setFb(Fb);

					theAlgorithm->setMaxIterNum(maxNumOfIterations);
					theAlgorithm->setFunctionTolerance(funcTol);
					theAlgorithm->setVariableTolerance(varTol);
					theAlgorithm->setFunctionType(1);

					double x0 = (a+b)/2.0;
					if (theAlgorithm->findZeroPoint(x0)>=0){
						zeroPoint = theAlgorithm->getZeroPoint();
						printf("zeroPoint is %f\n",zeroPoint);
						printf("iteration is %d\n",theAlgorithm->getIterationNumber());

						/* FMK
						if (printFlag ==5){
							// u space
							newPoint=  theProbabilityTransformation->get_u();
							//newPoint.addVector(0.0, *uPrime, 1.0);
							//newPoint.addVector(1.0, direction, zeroPoint);
							for (int i=0; i<numRV;i++){
								resultsOutputFile<< newPoint(i);
								resultsOutputFile<< "\n";
							}
							// x space 
							newPoint=  theProbabilityTransformation->get_x();
							//newPoint.addVector(0.0, *uPrime, 1.0);
							//newPoint.addVector(1.0, direction, zeroPoint);
							for ( i=0; i<numRV;i++){
								resultsOutputFile<< newPoint(i);
								resultsOutputFile<< "\n";
							} //for

							} //if printFlag ==5
						*/
						
			/*			// ===== Quan 2006 Dec. for FE divergence problem only(esp. SFSI problems) ============
						double b = theAlgorithm->getB();
						if (fabs(zeroPoint-b)<=varTol){  //possible diverge case
						
							newPoint = *uPrime;
							newPoint.addVector(1.0,direction,b);

							// Transform into original space
							result = theProbabilityTransformation->set_u(newPoint);
							if (result < 0) {
								opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
								<< " set the u-vector for xu-transformation. " << endln;
								return -1;
							}
								
							result = theProbabilityTransformation->transform_u_to_x();
							if (result < 0) {
								opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
									<< " transform u to x. " << endln;
								return -1;
							}
							x = theProbabilityTransformation->get_x();
							result = theGFunEvaluator->runGFunAnalysis(x);
							//numOfGFunEvaluations++;
							
							if(result<0) {  // diverge. evaluate value of current point 
								result = theGFunEvaluator->evaluateG(x);
								if (result < 0) {
									opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
										<< " tokenize limit-state function. " << endln;
									return -1;
								}
								gFunctionValue = theGFunEvaluator->getG();
							

									
								// divergence caused by other(like soil), but system safe
								if (gFunctionValue > funcTol){ 
									contribution = false;
									opserr << "divergence caused by soil, but system safe. gFunctionValue is:" <<gFunctionValue<<endln;
									opserr << "a is:" <<a<<"b is:"<<b<<endln;
								} 

								// divergence caused by system, system unsafe
							//	else {  
							//		zeroPoint = b;
							//	}
							
							} // if(result<0) {  // diverge
						}  //if     //for FE divergence problem       */
						// ========================================================
					}
					else {//does not converge
						printf("diverge!\n");
						contribution = false;
						// calling other like bisection method to make sure it converge!!
					}
				
			} //if (contribution &&needCallingZeroFinding)

				// compute contribution ...
				// ----- step 4, q=Fai[-h(u'_(n-1))] ----------
			if (contribution) {
				qk = aStdNormRV.getCDFvalue(-zeroPoint);
			} // if contribution.
			else { //	if (!contribution) 
				qk=0.0;
			}


	   } //if analysisTypeTag ==1 ------
	   
/////////////////////////////////////////////////////////////////////////////////////////////////////////
		else if (analysisTypeTag == 2) { //---------------

//=================== compute G1 zeropoint as before ================


			theAlgorithm->set_ii_2(0);
			theAlgorithm->setFunctionType(1); // do zerofinding for 1st zero...
			double G2 ;
			if (FEconvergence){
				G2 = this->getG2FromG1( gFuncAtProjectionPoint, littleDt);
				theAlgorithm->saveXG2(beta,G2);
				if (gFuncAtProjectionPoint<=0 && G2>=0) {contribution = false;}
			}

			if (fabs(gFuncAtProjectionPoint)<= funcTol && FEconvergence && contribution ) {
				zeroPoint = beta;  
				needCallingZeroFinding = false;
			}  // not need zeroFinding
				
			else if (gFuncAtProjectionPoint>funcTol && contribution) { //safe domain, must FEconvergence. this point is a and Fa
				
				a = beta;
				Fa = gFuncAtProjectionPoint;
					
				b = a + 2.0;
				newPoint = *uPrime;
				newPoint.addVector(1.0,direction,b);  // projected point on the orthogonal plane passing origin;
					  
				// Transform into original space
				/*
				result = theProbabilityTransformation->set_u(newPoint);
				if (result < 0) {
					opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
					<< " set the u-vector for xu-transformation. " << endln;
					return -1;
				}
					
				result = theProbabilityTransformation->transform_u_to_x();
				if (result < 0) {
					opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
						<< " transform u to x. " << endln;
					return -1;
				}
				x = theProbabilityTransformation->get_x();
				*/
				int result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
				if (result < 0) {
				  opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
					 << " transform u to x. " << endln;
				  return -1;
				}

				result = theGFunEvaluator->runAnalysis();
				numOfGFunEvaluations++;

				if(result<0) {  // diverge
					FEconvergence  =  false; 
					theAlgorithm->setFEconvergence(false);
					gFunctionValue =-1.0;
				}
				else { // converge
					gFunctionValue = theGFunEvaluator->evaluateExpression();
					theAlgorithm->saveXG1(b,gFunctionValue);
 
					G2 = this->getG2FromG1( gFunctionValue, littleDt);
					theAlgorithm->saveXG2(b,G2);

					if (gFunctionValue<=0 && G2>=0) {contribution = false;}
					
				} //else
				
				Fb = gFunctionValue;
				if (Fb>0) {// too far 
					if (G2 >= 0) {contribution = false;}
					else zeroPoint = 10.0; // far enough
					needCallingZeroFinding =false;
				}
			} //else if (gFuncAtProjectionPoint>funcTol)
				
			else if (gFuncAtProjectionPoint< -funcTol && contribution) { //unsafe domain, or diverge case

				b= beta; Fb = gFuncAtProjectionPoint;

				// find a, Fa
				double uPrimeNorm=uPrime->Norm();
				if (beta > uPrimeNorm){
					a= pow(beta*beta- uPrimeNorm*uPrimeNorm, 0.5);
					newPoint = *uPrime;
					newPoint.addVector(1.0,direction,a);  // projected point on the orthogonal plane passing origin;
				
					// Transform into original space
					/*
					result = theProbabilityTransformation->set_u(newPoint);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " set the u-vector for xu-transformation. " << endln;
						return -1;
					}
						
					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
					x = theProbabilityTransformation->get_x();
					*/
					result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}

					result = theGFunEvaluator->runAnalysis();
					numOfGFunEvaluations++;


					if(result<0) {  // diverge
						opserr<<"==================================================================== "<<endln;
						opserr<<"something wrong with orthogonalPlaneSamplingAnalysis. need save point "<<endln;
						opserr<<"==================================================================== "<<endln;
						resultsOutputFile1<< "a:" <<a<<","<<"F(a): diverge"<<". uPrime:"<<uPrime<<endln;
						
						needCallingZeroFinding =false;  
						contribution = false;
						
					// save point ....
					}
					else { // converge
						gFunctionValue = theGFunEvaluator->evaluateExpression();
						theAlgorithm->saveXG1(a,gFunctionValue);				
						G2 = this->getG2FromG1( gFunctionValue, littleDt);
						theAlgorithm->saveXG2(a,G2);

						if (gFunctionValue<=0 && G2>=0) {contribution = false;}

						Fa = gFunctionValue;
						if (fabs(Fa)<funcTol){
							zeroPoint=a; 
							needCallingZeroFinding =false;
						}
						else if (Fa<-funcTol) {
							opserr<<"==================================================================== "<<endln;
							opserr<<"something wrong with orthogonalPlaneSamplingAnalysis. need save point"<<endln;
							opserr<<"==================================================================== "<<endln;
							resultsOutputFile1<< "a:" <<a<<","<<"F(a)"<<Fa<<". uPrime:"<<uPrime<<endln;
							needCallingZeroFinding =false;
							contribution = false;  //correct
							// save point ....
						}	
					} //else converge

						
				} //if (beta > uPrimeNorm)
				else if (contribution) { //if (beta <= uPrimeNorm)
					a=0;
					newPoint = *uPrime;
					// Transform into original space
					/*
					result = theProbabilityTransformation->set_u(newPoint);
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " set the u-vector for xu-transformation. " << endln;
						return -1;
					}
					
					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
							<< " transform u to x. " << endln;
						return -1;
					}
					x = theProbabilityTransformation->get_x();
					*/
					result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
					if (result < 0) {
					  opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
						 << " transform u to x. " << endln;
					  return -1;
					}

					result = theGFunEvaluator->runAnalysis();
					numOfGFunEvaluations++;
					
					if(result<0) {  // diverge, ignore its contribution
					
						contribution = false;  
						needCallingZeroFinding =false;
					}
					else { // converge
						gFunctionValue = theGFunEvaluator->evaluateExpression();
						theAlgorithm->saveXG1(a,gFunctionValue);					
						G2 = this->getG2FromG1( gFunctionValue, littleDt);
						theAlgorithm->saveXG2(a,G2);
//						if ( G2>=0) {contribution = false;} ?????  Quan Dec 23 2007

						Fa = gFunctionValue;
						if (Fa<0) {   // ignore its contribution
							contribution = false;  
							zeroPoint = 0.0; 
							needCallingZeroFinding =false;
				
						}	
					} //else converge
				} // else { //if (beta <= uPrimeNorm)
			} //elseif (gFuncAtProjectionPoint<-funcTol) { //unsafe domain, 

				



 		    if (contribution && needCallingZeroFinding){
 

					// ======== call zeroFinding ==========
					theAlgorithm->setA(a);   theAlgorithm->setFa(Fa);
					theAlgorithm->setB(b);   theAlgorithm->setFb(Fb);

					theAlgorithm->setMaxIterNum(maxNumOfIterations);
					theAlgorithm->setFunctionTolerance(funcTol);
					theAlgorithm->setVariableTolerance(varTol);
					theAlgorithm->setFunctionType(1);

					double x0 = (a+b)/2.0;
					if (theAlgorithm->findZeroPoint(x0)>=0){
						zeroPoint = theAlgorithm->getZeroPoint();
						printf("zeroPoint is %f\n",zeroPoint);
						printf("iteration is %d\n",theAlgorithm->getIterationNumber());
						/* FMK
						if (printFlag ==5){
							newPoint= theProbabilityTransformation->get_u();
							//newPoint.addVector(0.0, *uPrime, 1.0);
							//newPoint.addVector(1.0, direction, zeroPoint);
							for (int i=0; i<numRV;i++){
								resultsOutputFile<< newPoint(i);
								resultsOutputFile<< "\n";
							}	
														// x space 
							newPoint=  theProbabilityTransformation->get_x();
							//newPoint.addVector(0.0, *uPrime, 1.0);
							//newPoint.addVector(1.0, direction, zeroPoint);
							for ( i=0; i<numRV;i++){
								resultsOutputFile<< newPoint(i);
								resultsOutputFile<< "\n";
							}

						}
						*/
					}
					else {//does not converge
						printf("diverge!\n");
						resultsOutputFile1<<"zero1 diverge!. uPrime:"<<uPrime<<endln;
						contribution = false;
						// calling other like bisection method to make sure it converge!!
					}
						
			} //  if (contribution && needCallingZeroFinding)
			
			if (contribution){
				qk = aStdNormRV.getCDFvalue(-zeroPoint);
			}
			else{	// no contribution
				qk=0.0;
			}
		
// ================================finding the second zero point ==========================================
			if (contribution){
				needCallingZeroFinding = true;
				theAlgorithm->setFEconvergence(true); // default value;

				// --- decide  bound a b----
				int ii_2 = theAlgorithm->get_ii_2();
				double x_2, G_2;

				a = -1.0;
				b = 999;

				for (int i=0; i<ii_2; i++){
					if (needCallingZeroFinding) {
						x_2 = theAlgorithm->getX2(i);
						G_2 = theAlgorithm->getG2(i);
						if (fabs(G_2)<= funcTol) {
							zeroPoint2 = x_2;
							needCallingZeroFinding = false;
						}
						else if (G_2>funcTol && x_2 >a){
							a = x_2;
							Fa = G_2;
						}
						else if (G_2<-funcTol && x_2 < b){
							b = x_2;
							Fb = G_2;
						}
					} //if (needCallingZeroFinding)
				
				}//for

			   
				if (needCallingZeroFinding){

				   if (a == -1.0){
						   
					    a=0;
						newPoint = *uPrime;
    					// Transform into original space
						/*
						result = theProbabilityTransformation->set_u(newPoint);
						if (result < 0) {
							opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
								<< " set the u-vector for xu-transformation. " << endln;
							return -1;
						}
							
						result = theProbabilityTransformation->transform_u_to_x();
						if (result < 0) {
							opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
								<< " transform u to x. " << endln;
							return -1;
						}
						x = theProbabilityTransformation->get_x();
						*/
						result = theProbabilityTransformation->transform_u_to_x(newPoint, x);
						if (result < 0) {
							opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
								<< " transform u to x. " << endln;
							return -1;
						}
						
						result = theGFunEvaluator->runAnalysis();
						numOfGFunEvaluations++;
						
						if(result<0) {  // diverge, ignore its contribution
							
							contribution = false;  //correct
							needCallingZeroFinding =false;
						}
						else { // converge
							gFunctionValue = theGFunEvaluator->evaluateExpression();							
							G2 = this->getG2FromG1( gFunctionValue, littleDt);
							theAlgorithm->saveXG2(a,G2);


							Fa = G2;
							if (Fa<0 && gFunctionValue>0) {   //  unsafe point
								zeroPoint2 = 0.0;    
								needCallingZeroFinding =false;
					
							}	//if (Fa<0)
							else if (gFunctionValue<=0) {contribution = false;}

						


						} //else converge
				   } // if (a==-1.0)

				   if ( b == 999) { 
						opserr<<"=== warning: something wrong "<<endln; 
						resultsOutputFile1<< "something wrong, G2: b does not exit.. uPrime:"<<uPrime<<endln;
						contribution = false;
				   }
						
					

				   if (contribution && needCallingZeroFinding){
					   					// ======== call zeroFinding ==========
						theAlgorithm->setA(a);   theAlgorithm->setFa(Fa);
						theAlgorithm->setB(b);   theAlgorithm->setFb(Fb);
						theAlgorithm->setMaxIterNum(maxNumOfIterations);
						theAlgorithm->setFunctionTolerance(funcTol);
						theAlgorithm->setVariableTolerance(varTol);
						theAlgorithm->setFunctionType(2);

						double x0 = (a+b)/2.0;
						if (theAlgorithm->findZeroPoint(x0)>=0){
							zeroPoint2 = theAlgorithm->getZeroPoint();
							printf("zeroPoint2 is %f\n",zeroPoint2);
							printf("iteration2 is %d\n",theAlgorithm->getIterationNumber());
							if (zeroPoint2 >= zeroPoint) {contribution = false;}  
							/*
							if (printFlag ==5){
								newPoint= theProbabilityTransformation->get_u();
								//newPoint.addVector(0.0, *uPrime, 1.0);
								//newPoint.addVector(1.0, direction, zeroPoint);
								for (int i=0; i<numRV;i++){
									resultsOutputFile3<< newPoint(i);
									resultsOutputFile3<< "\n";
								}
															// x space 
								newPoint=  theProbabilityTransformation->get_x();
								//newPoint.addVector(0.0, *uPrime, 1.0);
								//newPoint.addVector(1.0, direction, zeroPoint);
								for ( i=0; i<numRV;i++){
									resultsOutputFile<< newPoint(i);
									resultsOutputFile<< "\n";
								}

							}
							*/

						}
						else {//does not converge to zeroPoint 2.
							printf("========= diverge!\n something is wrong=========");
							resultsOutputFile1<< "something wrong, G2: can not find zero.. uPrime:"<<uPrime<<endln;
							contribution = false;
							// calling other like bisection method to make sure it converge!!
						}
					
				   } // if (contribution && needCallingZeroFinding)
				   
				}//if (needCallingZeroFinding)
	
			} //if (contribution) 
				// ----- step 4, q=Fai[-h(u'_(n-1))] ----------
				
			if (contribution){
				qk2 = aStdNormRV.getCDFvalue(-zeroPoint2);
				qk = (qk2-qk); 
				if (qk<0) {	resultsOutputFile1<< "something wrong, qk<0.. uPrime:"<<uPrime<<endln; qk=0.0;}
			}

			else {
				qk=0.0;	
			}
} //else if analysisTypeTag ==2 ------------------------


// ====================== post processing ====================
	q += qk;
	q_square += qk*qk;
	probability = q/k; // sample mean

			
	if (k ==1){
		cov = 999;
	}
	else if (probability >0) {
		cov = sqrt((q_square - 1/k*q*q) /(k*(k-1)))/probability;
	}
			
	// ---- inform user
	
	opserr.setPrecision(14);
	
	if (printFlag ==1) 
		opserr<<"\n=====probability:"<<probability<<"    COV:"<<cov<<"======\n"<<endln;
	else if (printFlag ==2) 
		opserr<<"\n=====Pf: "<<probability<<",       cov:"<<cov<<",       k:"<<k<<",     seed:"<<seed<<",    numOfGFunEvaluations: "<<numOfGFunEvaluations<<"======\n"<<endln;


	// ----

	if ((printFlag ==1)||(printFlag ==2)){
		sprintf(myString,"%12.6e     %12.6e ",probability,cov);
		resultsOutputFile << myString << " \n ";
		resultsOutputFile.flush();
	}		

		k++;
		isFirstSimulation = false;
	
		if (printFlag ==2){
 
			// write necessary data into file '_restart.tmp' .... close file
			ofstream resultsOutputFile5( "_restart.tmp");
			resultsOutputFile5<< probability <<endln;
			resultsOutputFile5<< cov         <<endln; 
			resultsOutputFile5<< k-1         <<endln;
			resultsOutputFile5<< seed        <<endln;
			resultsOutputFile5<< numOfGFunEvaluations <<endln;
			
			resultsOutputFile5.flush();
			resultsOutputFile5.close();
		}	

	}//while


	numOfSimulations = k-1; // base class member

	opserr << "Simulation Analysis completed." << endln;
	opserr << "numOfGFunEvaluations = " <<numOfGFunEvaluations<<"     numOfSimulations = "<<numOfSimulations<< endln;

	if ((printFlag ==1)||(printFlag ==2)){
		sprintf(myString,"numOfGFunEvaluations = %d    numOfSimulations = %d ",numOfGFunEvaluations,numOfSimulations);
		resultsOutputFile << myString << " \n ";
		resultsOutputFile.flush();
	}

// ------------- restart option 2007 Feb. -----
	if (printFlag ==2){
 
		// write necessary data into file '_restart.tmp' .... close file
		ofstream resultsOutputFile5( "_restart.tmp");
		resultsOutputFile5<< probability <<endln;
		resultsOutputFile5<< cov         <<endln; 
		resultsOutputFile5<< k-1         <<endln;
		resultsOutputFile5<< seed        <<endln;
		resultsOutputFile5<< numOfGFunEvaluations <<endln;
		
		resultsOutputFile5.flush();
		resultsOutputFile5.close();
	}
 

	resultsOutputFile.close();
	resultsOutputFile1.close();

	if (printFlag ==5 && this->analysisTypeTag ==2){
		resultsOutputFile3.flush();
		resultsOutputFile3.close();
	}




	return 0;
}

double OrthogonalPlaneSamplingAnalysis::getSampledValue(double x)
{
	

	if (! contribution){return 0.0;} // for both analysisType 1 and 2

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	
	Vector newPoint(numRV);
	Vector Xu(numRV);


	double factor = x/beta;
	newPoint.addVector(0.0, *uPrime,1.0);
	newPoint.addVector(1.0, *theDesignPointInUSpace,factor);


	/*
	double result = theProbabilityTransformation->set_u(newPoint);
	if (result < 0) {
		opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
		<< " set the u-vector for xu-transformation. " << endln;
		return -1;
	}
						
	result = theProbabilityTransformation->transform_u_to_x();
	if (result < 0) {
		opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
			<< " transform u to x. " << endln;
		return -1;
	}
	Xu = theProbabilityTransformation->get_x();
	*/
	int result = theProbabilityTransformation->transform_u_to_x(newPoint, Xu);
	if (result < 0) {
		opserr << "OrthogonalSamplingAnalysis::analyze() - could not " << endln
			<< " transform u to x. " << endln;
		return -1;
	}
				
	result = theGFunEvaluator->runAnalysis();
	numOfGFunEvaluations++;
	
	
	double gFunValue, G2;
	if(result<0) {  // diverge
		theAlgorithm->setFEconvergence(false);
		gFunValue = -1.0;
 
	}
	else { // converge
		theAlgorithm->setFEconvergence(true);
		gFunValue = theGFunEvaluator->evaluateExpression();
						
		// ------- save data into zeroFinding algorithm ---
		theAlgorithm->saveXG1(x, gFunValue);

		if (theAlgorithm->getFunctionType() ==2){
			G2 = this->getG2FromG1( gFunValue, littleDt);
			theAlgorithm->saveXG2(x,G2);
					
			if (gFunValue <=0 && G2>=0) {contribution = false;}

			gFunValue = G2;
		}

	}
	return gFunValue;
};



int OrthogonalPlaneSamplingAnalysis::getNumOfGFunEvaluations()
{
	return this->numOfGFunEvaluations;

}

void OrthogonalPlaneSamplingAnalysis::setSeed(int pSeed)
{
	seed=pSeed;
}

int OrthogonalPlaneSamplingAnalysis::getSeed()
{
	return seed;
}

/*int OrthogonalPlaneSamplingAnalysis::getNumOfSimulations()
{
	return numberOfSimulations;
}
*/

double OrthogonalPlaneSamplingAnalysis::getG2FromG1(double gFunctionValue, double littleDeltaT)
{
	
			/*
				LimitStateFunction * theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(1);

				if (theLimitStateFunction == 0) {
					opserr << "ImportanceSamplingAnalysis::analyze() - could not find" << endln
						<< " limit-state function with tag # 1." << endln;
					return -1;
				}

				

				char expressionG1[200];
				strcpy(expressionG1, theLimitStateFunction->getExpression());

				int nodeNumberAndDOFNumber[2];
				int i=0;
				char seps[]   = "}_";
				char *token;

//			    printf( "expressionG1: %s\n", expressionG1 );
 
			    token = strtok( expressionG1, seps );
				token = strtok( NULL, seps ); 

			    while( token != NULL )
				{
			
				  nodeNumberAndDOFNumber[i++]=atoi(token);
				  token = strtok( NULL, seps );
				}

 
				//delete [] expressionG1;

				char expression[200];   
				
				
				double ud=0;
				
				//sprintf(expression,"$ud_%d_%d",nodeNumberAndDOFNumber[0], nodeNumberAndDOFNumber[1]);		
				sprintf(expression,"set udot [nodeVel  %d  %d]",nodeNumberAndDOFNumber[0], nodeNumberAndDOFNumber[1]);		
				if (Tcl_Eval( theInterp, expression) !=TCL_OK ){opserr<<"wrong.."<<endln; exit(-1);};
				
				const char * myStr;
				myStr = Tcl_GetVar(theInterp, "udot",TCL_GLOBAL_ONLY );


				ud = atof(myStr);			

				if (fabs(ud) <1.e-14) {
					opserr<<" ----------------------------------------------------------------------------------" <<endln;
					opserr<<"warning:ImportanceSamplingAnalysis::analyze velocity is 0! \n probabily because you are running static case!"<<endln;
					opserr<<"-----------------------------------------------------------------------------------"<<endln;
				}
				
				// ---------- only for  G= u_lim-u case right now  -----------
                // in general: G2 = G1 + dG/du * udot * dt

				double g2FunctionValue = gFunctionValue-ud*littleDeltaT;


	//			opserr.precision(16);
	//			opserr<<"G1:"<<gFunctionValue<<".  G2: "<<g2FunctionValue<<endln;
				return g2FunctionValue;
*/ // June 2007 ---

  // This needs to be fixed -- MHS 10/7/2011
  //return theGFunEvaluator->getG2(gFunctionValue, littleDeltaT);
  return gFunctionValue; // So it compiles
}

bool OrthogonalPlaneSamplingAnalysis::getContribution()
{
	return contribution;
}
