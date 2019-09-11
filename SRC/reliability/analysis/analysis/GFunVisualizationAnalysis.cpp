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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-10-22 16:41:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/GFunVisualizationAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <GFunVisualizationAnalysis.h>
#include <FunctionEvaluator.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <GradientEvaluator.h>
#include <MeritFunctionCheck.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <RandomVariable.h>
#include <math.h>
#include <fstream>
#include <iomanip>
using std::ios;
using std::ifstream;
#include <ProbabilityTransformation.h>
#include <ReliabilityConvergenceCheck.h>
#include <RootFinding.h>

GFunVisualizationAnalysis::GFunVisualizationAnalysis(
					ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					TCL_Char *passedOutputFileName,
					TCL_Char *passedConvFileName,
					int passedConvResults,
					int passedSpace,
					int passedFunSurf,
					int passedAxes,
					int passedDir)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	theGFunEvaluator = passedGFunEvaluator;
	theProbabilityTransformation = passedProbabilityTransformation;
	theMeritFunctionCheck = 0;
	theGradGEvaluator = 0;
	theReliabilityConvergenceCheck = 0;

	strcpy(outputFileName,passedOutputFileName);
	strcpy(convFileName,passedConvFileName);

	convResults = passedConvResults;
	space = passedSpace;
	funSurf = passedFunSurf;
	axes = passedAxes;
	dir = passedDir;

	nrv = theReliabilityDomain->getNumberOfRandomVariables();

	scaleValue = 1.0;

	if (convResults == 1) {
		convFile.open( convFileName, ios::out );
	}


}

GFunVisualizationAnalysis::~GFunVisualizationAnalysis()
{
	if (convResults == 1) {
		convFile.close( );
	}
}


int
GFunVisualizationAnalysis::analyze(void)
{
	// Meaning of the tags
	// convResults [ 0:no,                   1:yes                 ]
	// space       [ 0:error,                1:X,        2:Y       ]
	// funSurf     [ 0:error,                1:function, 2:surface ]
	// dir         [ 0:(needed for surface), 1:rv,       2:file    ] (pass rvDir or theDirectionVector)
	// axes        [ 0:error,   1:coords1,   2:coords2,  3:file    ] (pass axesVector... or theMatrix+numPts)


	// Alert the user that the visualization analysis has started
	opserr << "Visualization Analysis is running ... " << endln;;


	// Open output file
	ofstream outputFile( outputFileName, ios::out );

	
	// Initial declarations
	int i,j;
	Vector iPoint(nrv);
	Vector fPoint(nrv);
	double result = 0;


	// Determine number of points to be plotted
	int numPoints1 = 1, numPoints2 = 1;
	if (axes==1) {
		numPoints1 = numPts1;
		numPoints2 = 1;
	}
	else if (axes==2) {
		numPoints1 = numPts1;
		numPoints2 = numPts2;
	}
	else if (axes==3) {
		// Should have a warning here if theMatrix isn't set
		numPoints1 = theMatrix.noCols()-1;
		numPoints2 = numLinePts;
	}

	int counter = 0; 
	for (i=1; i<=numPoints1; i++) {
		for (j=1; j<=numPoints2; j++) {

			counter++;
			opserr << counter << " ";

			// Get the current point in relevant space
			if (axes==1 || axes==2) {
				iPoint = this->getCurrentAxes12Point(i,j);
			}
			else if (axes==3) {
				iPoint = this->getCurrentAxes3Point(i,j);
			}

			// Evaluate G or find the surface
			if (funSurf==1) {

				// Transform the point into x-space (1-space) if the user has specified in 2-space
				if (space==2) {
				  /*
					result = theProbabilityTransformation->set_u(thePoint);
					if (result < 0) {
						opserr << "GFunVisualizationAnalysis::analyze() - " << endln
							<< " could not set u in the xu-transformation." << endln;
						return -1;
					}

					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "GFunVisualizationAnalysis::analyze() - " << endln
							<< " could not transform from u to x and compute Jacobian." << endln;
						return -1;
					}
					thePoint = theProbabilityTransformation->get_x();
				  */
				  result = theProbabilityTransformation->transform_u_to_x(iPoint, fPoint);
				  if (result < 0) {
				    opserr << "GFunVisualizationAnalysis::analyze() - " << endln
					   << " could not transform from u to x and compute Jacobian." << endln;
				    return -1;
				  }
				}

				// Evaluate g
				if (j==1) {
					result = this->evaluateGFunction(fPoint,true);
				}
				else {
					result = this->evaluateGFunction(fPoint,false);
				}
			}
			else if (funSurf==2) {

				// Find surface in relevant space
				result = this->findGSurface(iPoint);
			}

			// Print the result (g or distance) to file
			outputFile << result << " ";
		} 
		outputFile << endln;
	}

	// Clean up
	outputFile.close();

	// Print summary of results to screen (more here!!!)
	opserr << endln << "GFunVisualizationAnalysis completed." << endln;

	return 0;
}


Vector
GFunVisualizationAnalysis::getCurrentAxes12Point(int i, int j)
{
	// Initial declarations
	Vector iPoint(nrv);
	Vector fPoint(nrv);
	int result;
    // KRM to compile 4/22/2012
    bool startAtOrigin = false;

	// Find the start point in the space which the user 
	// wants to visualize in. 
	if (startAtOrigin) {
		

		// This indicates the origin in the standard normal space
	  iPoint.Zero();


		// If the user wants to visualize in the x-space; transform it into the x-space
		if (space==1) {
		  /*
			result = theProbabilityTransformation->set_u(thePoint);
			if (result < 0) {
				opserr << "GFunVisualizationAnalysis::analyze() - " << endln
					<< " could not set u in the xu-transformation." << endln;
				return -1;
			}

			result = theProbabilityTransformation->transform_u_to_x();
			if (result < 0) {
				opserr << "GFunVisualizationAnalysis::analyze() - " << endln
					<< " could not transform from u to x and compute Jacobian." << endln;
				return -1;
			}
			thePoint = theProbabilityTransformation->get_x();
		  */
		  result = theProbabilityTransformation->transform_u_to_x(iPoint, fPoint);
		  if (result < 0) {
		    opserr << "GFunVisualizationAnalysis::analyze() - " << endln
			   << " could not transform from u to x and compute Jacobian." << endln;
		    return -1;
		  }
		}
	}
	else {

		// Here the start point is actually given in the original space

	  //theReliabilityDomain->getStartPoint(iPoint);
	  for (int j = 0; j < nrv; j++) {
	    RandomVariable *theParam = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
	    iPoint(j) = theParam->getStartValue();
	  }

		// Transform it into the u-space if that's where the user wants to be
		if (space==2) {
		  /*
			result = theProbabilityTransformation->set_x(thePoint);
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
			thePoint = theProbabilityTransformation->get_u();
		  */
		  result = theProbabilityTransformation->transform_x_to_u(iPoint);
		  if (result < 0) {
		    opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			   << " could not transform from x to u." << endln;
		    return -1;
		  }
		}
	}

	// Now we have the point in the space we want it, 
	// So, set the random variables to be 'ranged'
	fPoint(rv1-1) = from1+(i-1)*interval1; 
	if (axes==2)
		fPoint(rv2-1) = from2+(j-1)*interval2;
	
	return fPoint;
}


Vector
GFunVisualizationAnalysis::getCurrentAxes3Point(int i, int j)
{
	// j is the index that runs along the line
	// i is to say which line we are on


	// Initial declarations
	Vector vector1(nrv);
	Vector vector2(nrv);


	// Extract end vectors
	for (int k=0; k<nrv; k++) {
		vector1(k) = theMatrix(k,i-1);
		vector2(k) = theMatrix(k,i);
	}


	// Compute vector between the points
	Vector vector = vector2-vector1;


	// The point to be returned
	Vector thePoint;
	double vectorFactor = ((double)(j-1))/((double)(numLinePts-1));
	thePoint = ( vector1 + vector*vectorFactor );


	return thePoint;
}


int
GFunVisualizationAnalysis::setDirection(int prvDir)
{
	rvDir = prvDir;

	return 0;
}

int
GFunVisualizationAnalysis::setDirection(Vector pDirectionVector)
{
	theDirectionVector = pDirectionVector;

	return 0;
}
	
int
GFunVisualizationAnalysis::setAxes(Vector axesVector)
{
	// Deschipher the content of the vector
	rv1         = (int)axesVector(0);
	from1       =      axesVector(1);
	double to1  =      axesVector(2);
	numPts1     = (int)axesVector(3);

	double to2 = 0;
	if (axesVector.Size() > 4 ) {
		rv2         = (int)axesVector(4);
		from2       =      axesVector(5);
		to2  =      axesVector(6);
		numPts2     = (int)axesVector(7);
	}

	interval1   = (to1-from1)/(numPts1-1);
	interval2   = (to2-from2)/(numPts2-1);

	return 0;
}

int
GFunVisualizationAnalysis::setAxes(Matrix pMatrix)
{
	theMatrix = pMatrix;

	return 0;
}

int
GFunVisualizationAnalysis::setNumLinePts(int pNumLinePts)
{
	numLinePts = pNumLinePts;

	return 0;
}

int
GFunVisualizationAnalysis::setRootFindingAlgorithm(RootFinding *pRootFinder)
{
	theRootFindingAlgorithm = pRootFinder;

	return 0;
}


int
GFunVisualizationAnalysis::setStartPoint(Vector *pStartPoint)
{
  //theStartPoint = pStartPoint;

	return 0;
}

int
GFunVisualizationAnalysis::setGradGEvaluator(GradientEvaluator *pGradGEvaluator)
{
	theGradGEvaluator = pGradGEvaluator;

	return 0;
}

int
GFunVisualizationAnalysis::setMeritFunctionCheck(MeritFunctionCheck *pMeritFunctionCheck)
{
	theMeritFunctionCheck = pMeritFunctionCheck;

	return 0;
}

int
GFunVisualizationAnalysis::setReliabilityConvergenceCheck(ReliabilityConvergenceCheck *pReliabilityConvergenceCheck)
{
	theReliabilityConvergenceCheck = pReliabilityConvergenceCheck;

	return 0;
}



double
GFunVisualizationAnalysis::findGSurface(Vector thePoint)
{

	// Initial declarations
	int result;
	double g;
	Vector surfacePoint(nrv);
	Vector Direction(nrv);
	Vector theTempPoint(nrv);
	int i;
	Vector distance;
	double scalarDist;
	double a;

	// Find direction; in whichever space user wants
	Direction.Zero();
	if (dir == 1)
		Direction(rvDir-1) = 1.0;
	else if (dir == 2) {
		if (nrv != theDirectionVector.Size()) {
			opserr << "ERROR: There is something wrong with the size of the direction" << endln
				<< " vector of the visualization analysis object." << endln;
		}
		Direction = theDirectionVector;
	}

	// Transform the point into x-space if the user has given it in 2-space
	if (space==2) {
	  /*
		result = theProbabilityTransformation->set_u(thePoint);
		if (result < 0) {
			opserr << "GFunVisualizationAnalysis::analyze() - " << endln
				<< " could not set u in the xu-transformation." << endln;
			return -1;
		}

		result = theProbabilityTransformation->transform_u_to_x();
		if (result < 0) {
			opserr << "GFunVisualizationAnalysis::analyze() - " << endln
				<< " could not transform from u to x and compute Jacobian." << endln;
			return -1;
		}
		theTempPoint = theProbabilityTransformation->get_x();
	  */
	  result = theProbabilityTransformation->transform_u_to_x(thePoint, theTempPoint);
	  if (result < 0) {
	    opserr << "GFunVisualizationAnalysis::analyze() - " << endln
		   << " could not transform from u to x and compute Jacobian." << endln;
	    return -1;
	  }
	}
	else {
		theTempPoint = thePoint;
	}


	// Evaluate limit-state function
	result = theGFunEvaluator->runAnalysis();
	if (result < 0) {
		opserr << "GFunVisualizationAnalysis::analyze() - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
		return -1;
	}
	g = theGFunEvaluator->evaluateExpression();


	// FIND THE POINT IN WHICHEVER SPACE USER HAS SPECIFIED
	surfacePoint = theRootFindingAlgorithm->findLimitStateSurface(space,g,Direction,thePoint);


	// POST-PROCESSING: FIND THE DISTANCE IN THE RELEVANT SPACE

	// Determine scaling factor 'a' in: surfacePoint = thePoint + NewtonDirection*a;
	i=0;
	while (Direction(i)==0.0) {
		i++;
	}
	a = ( surfacePoint(i) - thePoint(i) ) / Direction(i);
	distance = surfacePoint-thePoint;

	// Then the distance is:
	scalarDist = (a/fabs(a)*distance.Norm());


	return scalarDist;
}

double
GFunVisualizationAnalysis::evaluateGFunction(Vector thePoint, bool isFirstPoint)
{
	// Initial declarations
	double g;
	int result;

	Vector uPoint(thePoint);
	int nrv = thePoint.Size();
	Matrix Jxu(nrv, nrv);

	// Evaluate limit-state function
	result = theGFunEvaluator->runAnalysis();
	if (result < 0) {
		opserr << "GFunVisualizationAnalysis::analyze() - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
		return -1;
	}
	g = theGFunEvaluator->evaluateExpression();


	// Possibly compute and print out more comprehensive results
	if (convResults==1) {

		// Initial declarations
		char myString[300];
		Vector u;
		double meritFunctionValue;
		int numCrit = theReliabilityConvergenceCheck->getNumberOfCriteria();
		double criteriumValue;


		// Print a division if this is the first point on the vector between two points
		if (isFirstPoint) {
			double dummy = 0.0;
			sprintf(myString,"%20.14e  %20.14e  ", dummy,dummy);
			convFile << myString;
			sprintf(myString,"%20.14e  ", dummy);
			convFile << myString;
			for (int crit=1; crit<=numCrit; crit++) {
				sprintf(myString,"%20.14e  ", dummy);
				convFile << myString;
			}
			convFile << endln;	
		}

		// Set scale value of convergence criteria
		if (scaleValue == 1.0 && convResults == 1) {
			scaleValue = g;
			theReliabilityConvergenceCheck->setScaleValue(g);
		}


		// Transform the x-point into u-space
		/*
		result = theProbabilityTransformation->set_x(thePoint);
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

		result = theProbabilityTransformation->transform_x_to_u(uPoint);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from x to u." << endln;
			return -1;
		}
	
		/*
		// And back to the x-space again, just to get the jacobian...
		result = theProbabilityTransformation->set_u(u);
		if (result < 0) {
			opserr << "ArmijoStepSizeRule::computeStepSize() - could not set " << endln
				<< " vector u in the xu-transformation. " << endln;
			return -1;
		}
		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "ArmijoStepSizeRule::computeStepSize() - could not  " << endln
				<< " transform u to x. " << endln;
			return -1;
		}
		theProbabilityTransformation->get_x();
		Matrix jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();
		*/
		result = theProbabilityTransformation->getJacobian_x_to_u(Jxu);

		// Print limit-state fnc value and distance to origin to file
		sprintf(myString,"%20.14e  %20.14e  ", g,u.Norm());
		convFile << myString;

		// Compute gradients
		result = theGradGEvaluator->computeGradient(g);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute gradients of the limit-state function. " << endln;
			return -1;
		}
		Vector gradientOfgFunction = theGradGEvaluator->getGradient();
		//gradientOfgFunction = jacobian_x_u ^ gradientOfgFunction;
		gradientOfgFunction = Jxu ^ gradientOfgFunction;

		
		// Evaluate the merit function
		if (isFirstPoint) {
			meritFunctionValue = theMeritFunctionCheck->updateMeritParameters(u, g, gradientOfgFunction);
		}
		meritFunctionValue = theMeritFunctionCheck->getMeritFunctionValue(u, g, gradientOfgFunction);
		sprintf(myString,"%20.14e  ", meritFunctionValue);
		convFile << myString;

		// Get value of convergence criteria
		theReliabilityConvergenceCheck->check(u,g,gradientOfgFunction);
		for (int crit=1; crit<=numCrit; crit++) {
			criteriumValue = theReliabilityConvergenceCheck->getCriteriaValue(crit);
			sprintf(myString,"%20.14e  ", criteriumValue);
			convFile << myString;
		}

		convFile << endln;

	}


	return g;
}

