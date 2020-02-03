// RespSurfaceSimulation.cpp: implementation of the RespSurfaceSimulation class.

//
//////////////////////////////////////////////////////////////////////

#include "RespSurfaceSimulation.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <MatrixOperations.h>

using std::ios;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RespSurfaceSimulation::RespSurfaceSimulation( ReliabilityDomain * pReliabilityDomain,
						Vector * pDesignPoint,
						char * pFileName,
						SurfaceDesign * pSurfaceDesign,
						RandomNumberGenerator * pRandomNumberGenerator,
                        double pTargetCOV,
						int pNumberOfSimulations, bool pIsTimeVariant, double pLittleDtOverDt,double pNormOfGradient)
{
	theReliabilityDomain = pReliabilityDomain;
 	theDesignPoint= new Vector(*pDesignPoint);
    strcpy(fileName,pFileName);
	theSurfaceDesign = pSurfaceDesign;
	theRandomNumberGenerator = pRandomNumberGenerator;
    targetCOV = pTargetCOV;
	numOfSimulations =pNumberOfSimulations;
	isTimeVariant = pIsTimeVariant;
	littleDtOverDt = pLittleDtOverDt;
	normOfGradient = pNormOfGradient;

}

RespSurfaceSimulation::~RespSurfaceSimulation()
{
	delete theDesignPoint;
}

int RespSurfaceSimulation::runSimulationAnalysis()
{

	opserr << "Sampling Analysis is running ... " << endln;

	
	double gFunValue=0, gFunValue2=0, phi=0, q=0, q_bar=0,h=0, variance_of_q_bar=0, cov_of_q_bar=0 ;
	int result=0;
	int k = 1;
	int I=0;
	int seed = 1;
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	Vector u(numRV);
	Vector uPrime(numRV);

 	bool failureHasOccured = false;

	double sum_q=0.0;
	double sum_q_squared=0.0;

	ofstream resultsOutputFile( fileName, ios::out );
	ofstream out("ignoredPoints.out", ios::out);
	
	
	// ======  samplingStdv ======
	double samplingStdv =1.0;

	Matrix covariance(numRV, numRV);
	Matrix chol_covariance(numRV, numRV);
	Matrix inv_covariance(numRV, numRV);

	for (int i=0;  i<numRV;  i++) {
		for (int j=0;  j<numRV;  j++) {
			if (i==j) {
				covariance(i,j) = samplingStdv*samplingStdv;
			}
			else
				covariance(i,j) = 0.0;
		}
	}

	MatrixOperations *theMatrixOperations = new MatrixOperations(covariance);

	result = theMatrixOperations->computeLowerCholesky();
	if (result < 0) {
		opserr << "ImportanceSamplingAnalysis::analyze() - could not compute" << endln
			<< " the Cholesky decomposition of the covariance matrix." << endln;
		return -1;
	}
	chol_covariance = theMatrixOperations->getLowerCholesky();


	// Inverse of covariance matrix
	result = theMatrixOperations->computeInverse();
	if (result < 0) {
		opserr << "ImportanceSamplingAnalysis::analyze() - could not compute" << endln
			<< " the inverse of the covariance matrix." << endln;
		return -1;
	}
	inv_covariance = theMatrixOperations->getInverse();


	// Compute the determinant, knowing that this is a diagonal matrix
	result = theMatrixOperations->computeTrace();
	if (result < 0) {
		opserr << "ImportanceSamplingAnalysis::analyze() - could not compute" << endln
			<< " the trace of the covariance matrix." << endln;
		return -1;
	}
	double det_covariance = theMatrixOperations->getTrace();
	

	// Pre-compute some factors to minimize computations inside simulation loop
	double pi = 3.14159265358979;
	double factor1 = 1.0 / ( pow((2.0*pi),((double)numRV/2.0)) );
	double factor2 = 1.0 / ( pow((2.0*pi),((double)numRV/2.0)) * sqrt(det_covariance) );


// =====



    double  beta  = theDesignPoint->Norm();
	double criterion = 1.0* exp( -0.5 * beta*beta);


// -- translation to get second dp. for time variant case
	
	

	Vector * theDesignPoint2 =0;
	Vector * dp2Prime =0;
	Vector * gradG2=0;


	if (isTimeVariant) {

		if (normOfGradient ==0) {
			opserr<<"Fatal: normOfGradient is 0 in RespSurfaceSimulation::runSimulationAnalysis()"<<endln;
			exit(-1);
		
		}
		theDesignPoint2 = new Vector(theDesignPoint->Size());


		(*theDesignPoint2)(0) = (*theDesignPoint)(0) - littleDtOverDt * ( (*theDesignPoint)(0)-0.0 );

		for (int j=1; j<theDesignPoint2->Size(); j++) {
			(*theDesignPoint2)(j) = (*theDesignPoint)(j) - littleDtOverDt * ( (*theDesignPoint)(j)-(*theDesignPoint)(j-1) );
		}

		//Vector alpha2= theDesignPoint2/theDesignPoint2.Norm();

		dp2Prime = new Vector(*theDesignPoint2);   // second dp in local coord system
		dp2Prime->addVector(1.0, (*theDesignPoint),-1.0);

		gradG2 =  new Vector((*theDesignPoint2)*(-normOfGradient)/theDesignPoint2->Norm()); // gradient G2 = -|gradG|* alpha2;  ?????????????????    Quan
	}

// ---



	double govCov =1.0;

	bool isFirstSimulation = true;
	while( (k<=numOfSimulations && govCov>targetCOV || k<=2) ) {

		

		// Create array of standard normal random numbers
		if (isFirstSimulation) {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV,seed);
		}
		else {
			result = theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numRV);
		}

		seed = theRandomNumberGenerator->getSeed();
		if (result < 0) {
			opserr << "ImportanceSamplingAnalysis::analyze() - could not generate" << endln
				<< " random numbers for simulation." << endln;
			return -1;
		}
		uPrime = theRandomNumberGenerator->getGeneratedNumbers();  // in translated space with origin in u*


		u= *theDesignPoint + chol_covariance *uPrime;

		if (u.Norm()<theDesignPoint->Norm()){
			gFunValue =1.0;
		}
		else {

		   gFunValue = theSurfaceDesign->getFunctionValue(&uPrime); /// BE CAREFUL THAT INPUT IS UPRIME, IN COORD SYSTEM WITH ORIGIN AT U*, WITHOUT ROTATION.

	

		}

		if (! isTimeVariant) {
			if (gFunValue <= 0.0) {
				I = 1;
				failureHasOccured = true;
			}
			else {
				I = 0;
			}
		}
		else {  // time variant


			gFunValue2 = theSurfaceDesign->getFunctionValue2(&uPrime, dp2Prime, gradG2); /// BE CAREFUL THAT INPUT IS UPRIME, IN COORD SYSTEM WITH ORIGIN AT U*, WITHOUT ROTATION.
		//	double useless = theSurfaceDesign->debug(&uPrime, dp2Prime, gradG2); 
			if ((gFunValue > 0.0) &&(gFunValue2 <= 0.0)) {
				I = 1;
				failureHasOccured = true;
			}
			else {
				I = 0;
			}

		}

				// Compute values of joint distributions at the u-point
		// phi = exp( -0.5 * (u ^ u) );
		// h   = exp( -0.5 * (uPrime^uPrime) );

		phi = factor1 * exp( -0.5 * (u ^ u) );
		Vector temp1 = inv_covariance ^ uPrime;
		double temp2 = temp1 ^ uPrime;
		h   = factor2 * exp( -0.5 * temp2 );

		// Update sums
		 q = I * phi / h;
// === check whether it is not good point ===


	/*	 double sum_q_tmp = sum_q + q;
		 double sum_q_squared_tmp = sum_q_squared + q*q;
		 double cov_of_q_bar_tmp = govCov;
		 
		 if (sum_q_tmp > 0.0) {
				double q_bar_tmp = 1.0/(double)k * sum_q_tmp;
				double variance_of_q_bar_tmp = 1.0/(double)k *
					( 1.0/(double)k * sum_q_squared_tmp - (sum_q_tmp /(double)k)*(sum_q_tmp /(double)k));
				if (variance_of_q_bar_tmp  < 0.0) {
					variance_of_q_bar_tmp  = 0.0;
				}
				cov_of_q_bar_tmp  = sqrt(variance_of_q_bar_tmp ) / q_bar_tmp ;
		}
*/
  
	 if ((q > criterion) /*|| ((I==1) &&(k>100)&&(cov_of_q_bar_tmp > 2.0* govCov))*/) {// ignore this point
//			 if ((q > 2*fabs(q)) /*|| ((I==1) &&(k>100)&&(cov_of_q_bar_tmp > 2.0* govCov))*/) {// ignore this point
			opserr<<"warning: bad point "<<k<<" is ignored! valid for ONLY one limit state function! "<<endln;  
//			out<<"===k:"<<k<<",   cov_of_q_bar_tmp: "<<cov_of_q_bar_tmp<<",  govCov:"<<govCov<<endln;
			out<<" point in u space:\n";
			for(int i =0; i< u.Size(); i++)    out<< u(i)<<endln;
	 
			
		 
		 }


		 else {
			sum_q = sum_q + q;
			sum_q_squared = sum_q_squared + q*q;


			if (sum_q > 0.0) {
				q_bar = 1.0/(double)k * sum_q;
				variance_of_q_bar = 1.0/(double)k *
					( 1.0/(double)k * sum_q_squared - (sum_q /(double)k)*(sum_q /(double)k));
				if (variance_of_q_bar  < 0.0) {
					variance_of_q_bar  = 0.0;
				}
				cov_of_q_bar  = sqrt(variance_of_q_bar ) / q_bar ;
			}


			// Compute governing coefficient of variation
			if (failureHasOccured) {
				if (cov_of_q_bar>1.0e-14)
					govCov = cov_of_q_bar;
			}

			resultsOutputFile<<"k:"<<k<<",  Pf: "<<q_bar<<", cov:"<<govCov<<endln;
			
			int m=k%10000;
			if (m==0)
				opserr << "Sample #" << k << ",  " <<",  Pf: "<<q_bar<<", cov:"<<govCov<<endln;

 			k++;
		 }
		isFirstSimulation = false;

	} // while( (k<=numberOfSimulations && govCov>targetCOV || k<=2) )
	k--;

	out.close();	
	resultsOutputFile.close();
	failureProbability = q_bar;
	theCOV = govCov;

	if (isTimeVariant){
		delete theDesignPoint2;
	    delete dp2Prime;
	    delete gradG2;
	
	}


	return 0;
}





double RespSurfaceSimulation::getFailureProbability()
{
	return failureProbability;
}



double RespSurfaceSimulation::getCov()
{
	return theCOV;
}

void RespSurfaceSimulation::setNormOfGradientU(double normU)
{
	normOfGradient = normU;
}
