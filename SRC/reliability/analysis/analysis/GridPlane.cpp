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
 

#include <GridPlane.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////
  
GridPlane::GridPlane(int pTag, Vector * pDesignPt, 
							   PrincipalAxis * firstAxis, PrincipalAxis * secondAxis, Matrix * pRotation, 
							   ProbabilityTransformation * pProbabilityTransformation,
							   FunctionEvaluator * pGFunEvaluator )
{
	numOfPlane = pTag;


	if (pDesignPt !=0)
		theDesignPt=new Vector(*pDesignPt);
	else theDesignPt =0;

	
	theFirstAxis = firstAxis;
	theSecondAxis= secondAxis;
	Rotation = pRotation;
	
    theProbabilityTransformation = pProbabilityTransformation;
	theGFunEvaluator = pGFunEvaluator;
 

	gridValues =0;
	gridValuesG2 =0;
    FEConvergence = true;

 
}

GridPlane::~GridPlane()
{
	if (theDesignPt !=0) delete theDesignPt;
	if (gridValues !=0) delete gridValues;

}




double GridPlane::getValueOnGrid(double p_x, double p_y)
{
    FEConvergence = true;		
	if ( (fabs(p_x) <1.0e-14) && (fabs(p_y) <1.0e-14)) return 0.0;

	Vector u = this->getPtInOrigUSpace(p_x,p_y);
	/*
	theProbabilityTransformation->set_u(u);
	int result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
	if (result < 0) {
		opserr << "Hessian::formHessianByFDM - " << endln
			<< " could not transform from u to x." << endln;
		return -1;
	}
	Vector x = theProbabilityTransformation->get_x();
	*/
	Vector x;
	theProbabilityTransformation->transform_u_to_x(u, x);		
 

	int result = theGFunEvaluator->runAnalysis();
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
		FEConvergence = false;
		return -1;
	}
	double gFunctionValue = theGFunEvaluator->evaluateExpression();

	return gFunctionValue;
}


Vector GridPlane::getPtInOrigUSpace(double x, double y)
{
// U=u* + x*theEigenVector/|theEigenVector|+y*designPt/|designPt|. evaluate by G FunctionEvaluator
// which is same as:  U=u* + Rotation*[x,y]'

	if (theFirstAxis ==0 || theSecondAxis ==0 || theDesignPt==0) { 
		opserr<<"PrincipalPlane::getValueOnGrid wrong. theEigenVector or theDesignPt not exist"<<endln; 
		exit(-1);
	}

	Vector * firstAxis = theFirstAxis->getAxisDirection();
    Vector * secondAxis= theSecondAxis->getAxisDirection();

	double  norm1=firstAxis->Norm();
	double  norm2=secondAxis->Norm();

	static Vector tmp(*theDesignPt);
	tmp.addVector(0.0, *theDesignPt, 1.0);

	
	if (fabs(norm1)<1.0e-16 || fabs(norm2)<1.0e-16){
		opserr<<"PrincipalPlane::getValueOnGrid wrong. theSecondAxis or theFirstAxis is zero"<<endln; 
		exit(-1);
	}

	tmp.addVector(1.0, *firstAxis, x/norm1);
    tmp.addVector(1.0, *secondAxis, y/norm2);
 

	return tmp;

}


int GridPlane::setDesignPoint(Vector *pDesignPt)
{

	theDesignPt = pDesignPt;
	return 0;
}

int GridPlane::setRotation(Matrix *R)
{
	Rotation =R;
	return 0;
}

 

double GridPlane::getBeginOfGridX()
{
   return (theFirstAxis->getExperimentalPointRule())->getBeginOfGrid();
}

double GridPlane::getEndOfGridX()
{
   return (theFirstAxis->getExperimentalPointRule())->getEndOfGrid();
}

int GridPlane::getNumOfGridXPt()
{
	return (theFirstAxis->getExperimentalPointRule())->getNumberOfPoints();
}

double GridPlane::getBeginOfGridY()
{
   return (theSecondAxis->getExperimentalPointRule())->getBeginOfGrid();
}

double GridPlane::getEndOfGridY()
{
   return (theSecondAxis->getExperimentalPointRule())->getEndOfGrid();
}

int GridPlane::getNumOfGridYPt()
{
	return (theSecondAxis->getExperimentalPointRule())->getNumberOfPoints();
}


int GridPlane::cleanGridMatrix()
{
	if (gridValues !=0)
		gridValues->Zero();
	return 0;
}

 
int GridPlane::getNumOfPlane()
{
	return numOfPlane;
}



Matrix * GridPlane::getGridValuesPtr()
{
	return gridValues;
}

Matrix * GridPlane::getRotationPtr()
{
	return Rotation;
}

PrincipalAxis * GridPlane::getAxisPtr(int i)
{
	if (i==1)
	    return theFirstAxis;
	else if (i==2)
		return theSecondAxis;
	else {
		opserr<<"Fatal: GridPlane::getAxisPtr(i) , i is : "<<i<<endln;
		exit(-1);
	} 

	
}

int GridPlane::setNumOfPlane(int num)
{
	numOfPlane=num;
	return 0;
}


 

Vector GridPlane::getPointCoordinate(int i, int j)
{
	static Vector tmp(2);

	tmp(0) = (theFirstAxis->getExperimentalPointRule())->getPointCoordinate(i);
	tmp(1)  = (theSecondAxis->getExperimentalPointRule())->getPointCoordinate(j);

	
	return tmp;
}



int GridPlane::getnumOfAxis(int i){

	if (i==1) return theFirstAxis->getNumOfAxis();
	else if (i==2) return theSecondAxis->getNumOfAxis();

	return -1;

};

double GridPlane::getValueOnGrid(Vector &point)
{
    double x= point(0);
	double y = point(1);
	return getValueOnGrid(x, y);
}


int GridPlane::setAxis( int pAxis, PrincipalAxis * theAxis){

	if (pAxis ==1){
		theFirstAxis = theAxis;
	}
	else if (pAxis ==2){
		theSecondAxis = theAxis;
	}
	else {
	    opserr<<"Fatal: GridPlane::setAxis( ). pAxis="<<pAxis<<endln;
		exit(-1);
	}

	return 0;

 
};

int GridPlane::setGridValue(int i, int j, double v)
{

	if (gridValues ==0){
		int ii = (theFirstAxis->getExperimentalPointRule())->getNumberOfPoints();
	    int jj = (theSecondAxis->getExperimentalPointRule())->getNumberOfPoints();
		gridValues = new Matrix(ii,jj);
		gridValues->Zero();
	
	}

	(*gridValues)(i,j) = v;

	return 0;
}

int GridPlane::setGridValueG2(int i, int j, double v)
{

	if (gridValuesG2 ==0){
		int ii = (theFirstAxis->getExperimentalPointRule())->getNumberOfPoints();
	    int jj = (theSecondAxis->getExperimentalPointRule())->getNumberOfPoints();
		gridValuesG2 = new Matrix(ii,jj);
		gridValuesG2->Zero();
	
	}

	(*gridValuesG2)(i,j) = v;

	return 0;
}

bool GridPlane::isFEConverged()
{
	return FEConvergence;
}

double GridPlane::getSavedValueOnGrid(int i, int j)
{
	if (gridValues ==0) {opserr<<"GridPlane::getSavedValueOnGrid(), gridValues matrix not exist \n"<<endln; exit(-1);}
	return (*gridValues)(i,j) ;
}

double GridPlane::getSavedValueG2OnGrid(int i, int j)
{
	if (gridValuesG2 ==0) {opserr<<"GridPlane::getSavedValueG2OnGrid(), gridValues matrix not exist \n"<<endln; exit(-1);}
	return (*gridValuesG2)(i,j) ;
}


double GridPlane::getValueG2OnGrid(double p_x, double p_y, double valueOfG, double littleDt)
{

	if (!FEConvergence){ return -1.0;}
	double gFunctionValue2;


	if ((fabs(p_x) < 1.0e-14) &&(fabs(p_y) < 1.0e-14)) { // origin
		FEConvergence = true;
		Vector u = this->getPtInOrigUSpace(p_x,p_y);
		/*
		theProbabilityTransformation->set_u(u);
		int result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "Hessian::formHessianByFDM - " << endln
				<< " could not transform from u to x." << endln;
			return -1;
		}
			
		Vector x = theProbabilityTransformation->get_x();
		*/
		Vector x;
		theProbabilityTransformation->transform_u_to_x(u, x);

		int result = theGFunEvaluator->runAnalysis();
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not run analysis to evaluate limit-state function. " << endln;
			FEConvergence = false;
			exit(-1);
		}
		valueOfG = theGFunEvaluator->evaluateExpression();

		// This needs to be fixed -- MHS 10/7/2011
		//gFunctionValue2 = theGFunEvaluator->getG2(valueOfG, littleDt);
		//gFunctionValue2 -= valueOfG;
		gFunctionValue2 = valueOfG; // So it will compile -- MHS
	}
	else {
		// This needs to be fixed -- MHS 10/7/2011
	  //gFunctionValue2 = theGFunEvaluator->getG2(valueOfG, littleDt);
	  gFunctionValue2 = valueOfG; // So it will compile -- MHS
	}

	
	FEConvergence = true;
	return gFunctionValue2;
	
}
 

Matrix * GridPlane::getGridValuesG2Ptr()
{
	return gridValuesG2;
}

