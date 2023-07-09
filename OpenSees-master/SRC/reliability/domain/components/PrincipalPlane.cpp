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

#include <PrincipalPlane.h>
#include <Matrix.h>
#include <math.h>

PrincipalPlane::PrincipalPlane(int pTag, Vector * pDesignPt, 
			       Vector * pEigenVector, Matrix * pRotation, 
			       ProbabilityTransformation * pProbabilityTransformation,
			       FunctionEvaluator * pGFunEvaluator, double pCurv)
{
	if (pDesignPt !=0)
		theDesignPt=new Vector(*pDesignPt);
	else theDesignPt =0;

	if (pEigenVector !=0)
		theEigenVector=new Vector(*pEigenVector);
	else theEigenVector =0;

	if (pRotation !=0)
		Rotation=new Matrix(*pRotation);
	else Rotation =0;

 
	
	numOfPlane = pTag;
    theProbabilityTransformation = pProbabilityTransformation;
	theGFunEvaluator = pGFunEvaluator;
	theCurvature = pCurv;
	FEConvergence = true;
	gridValuesG2 =0;
	gridValues =0;
	//  -- default
	beginOfGridX=-2.0;
	endOfGridX = 2.0;
    numOfGridXPt = 10;

	beginOfGridY=-2.0;
	endOfGridY = 2.0;
    numOfGridYPt = 10;




}

PrincipalPlane::~PrincipalPlane()
{
	if (theDesignPt !=0) delete theDesignPt;
	if (theEigenVector !=0) delete theEigenVector;
	if (Rotation !=0) delete Rotation;
	if (gridValues !=0) delete gridValues;


}

double PrincipalPlane::getValueOnGrid(double p_x, double p_y)
{
		
  FEConvergence = true;
  if ((fabs(p_x) < 1.0e-14) &&(fabs(p_y) < 1.0e-14)) {
    return 0.0; 
  }
  
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
  theProbabilityTransformation->transform_u_to_x(u,x);
  
  // --- output ----
  //	opserr<<"u:"<<u<<endln;
  //	opserr<<"x:"<<x<<endln;
  //	opserr<<"======"<<endln;
  
  // --- compute standard sensitivity --
  //	jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();
  
  // ----remove sensitivity algorithm ---
  
  int  result = theGFunEvaluator->runAnalysis();
  if (result < 0) {
    opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
	   << " could not run analysis to evaluate limit-state function. " << endln;
    FEConvergence = false;
    return -1;
  }
  double gFunctionValue = theGFunEvaluator->evaluateExpression();
  
  return gFunctionValue;
}

Vector 
PrincipalPlane::getPtInOrigUSpace(double x, double y)
{
  // U=u* + x*theEigenVector/|theEigenVector|+y*designPt/|designPt|. evaluate by G FunctionEvaluator
  // which is same as:  U=u* + Rotation*[x,y]'
  
  if (theEigenVector ==0 || theDesignPt==0) { 
    opserr<<"PrincipalPlane::getValueOnGrid wrong. theEigenVector or theDesignPt not exist"<<endln; 
    exit(-1);
  }
  
  double  norm1=theEigenVector->Norm();
  double  norm2=theDesignPt->Norm();
  static Vector tmp(*theDesignPt);
  tmp.addVector(0.0, *theDesignPt, 1.0);
  
  
  if (fabs(norm1)<1.0e-16 || fabs(norm2)<1.0e-16){
    opserr<<"PrincipalPlane::getValueOnGrid wrong. theEigenVector or theDesignPt is zero"<<endln; 
    exit(-1);
  }
  
  tmp.addVector(1.0, *theEigenVector, x/norm1);
  tmp.addVector(1.0, *theDesignPt, y/norm2);
  
  
  return tmp;
}

int PrincipalPlane::setGridXInfo(int numOfGrid, double begin, double end)
{
	numOfGridXPt = numOfGrid;
	beginOfGridX = begin;
	endOfGridX = end;
	return 0;
}

int PrincipalPlane::setGridYInfo(int numOfGrid, double begin, double end)
{
	numOfGridYPt = numOfGrid;
	beginOfGridY = begin;
	endOfGridY = end;
	return 0;
}

int PrincipalPlane::setDesignPoint(Vector *pDesignPt)
{
	if (theDesignPt==0)
		theDesignPt = new Vector(*pDesignPt);
	else 
		theDesignPt->addVector(0.0, *pDesignPt, 1.0);
	return 0;
}

int PrincipalPlane::setRotation(Matrix *R)
{
	if (Rotation==0)
		Rotation = new Matrix(*R);
	else 
		Rotation->addMatrix(0.0, *R, 1.0);

	Rotation = R;
	return 0;
}

int PrincipalPlane::setEigenVector(Vector *pEigenV)
{
	if (theEigenVector ==0)
		theEigenVector = new Vector(*pEigenV);
	else
		theEigenVector->addVector(0.0, *pEigenV, 1.0);

	return 0 ;
}

int PrincipalPlane::setCurvature(double pCurv)
{
	theCurvature = pCurv;
	return 0;
}

int PrincipalPlane::getNumOfGridXPt()
{
	return numOfGridXPt;
}

double PrincipalPlane::getBeginOfGridX()
{
	return beginOfGridX;
}

double PrincipalPlane::getEndOfGridX()
{
	return endOfGridX;
}

int PrincipalPlane::getNumOfGridYPt()
{
	return numOfGridYPt;
}

double PrincipalPlane::getBeginOfGridY()
{
	return beginOfGridY;
}

double PrincipalPlane::getEndOfGridY()
{
	return endOfGridY;
}

int PrincipalPlane::setGridValue(int i, int j, double value)
{
   if (gridValues ==0) 
		gridValues = new Matrix(numOfGridXPt, numOfGridYPt);

	(*gridValues)(i,j)= value;
	return 0;
}

int PrincipalPlane::cleanGridMatrix()
{
	if (gridValues !=0)
		gridValues->Zero();
	return 0;
}

int PrincipalPlane::copyValues(PrincipalPlane *another)
{	
/*	this->beginOfGridX = another->getBeginOfGridX();
	this->beginOfGridY = another->getBeginOfGridY();
	this->endOfGridX = another->getEndOfGridX();
	this->endOfGridY=another->getEndOfGridY();
	this->numOfGridXPt = another->getNumOfGridXPt();
	this->numOfGridYPt = another->getNumOfGridYPt(); 
	this->numOfPlane = another ->getNumOfPlane();  */
	this->theCurvature = another ->getCurvature();
	
/*	if ((this->gridValues !=0) && (another->getGridValuesPtr()!=0))
		this->gridValues->addMatrix(0.0, *(another->getGridValuesPtr()), 1.0);
	else if ((this->gridValues ==0) && (another->getGridValuesPtr()!=0))
		gridValues = new Matrix(*(another->getGridValuesPtr()));
	else if ((this->gridValues !=0) && (another->getGridValuesPtr()==0)) {
		delete gridValues;
		gridValues=0;
	}
*/
	if ((this->Rotation !=0) && (another->getRotationPtr()!=0))
		this->Rotation->addMatrix(0.0, *(another->getRotationPtr()), 1.0);
	else if ((this->Rotation ==0) && (another->getRotationPtr()!=0))
		Rotation = new Matrix(*(another->getRotationPtr()));
	else if ((this->Rotation !=0) && (another->getRotationPtr()==0)) {
		delete Rotation;
		Rotation=0;
	}


	if ((this->theEigenVector !=0) && (another->getEigenVectorPtr()!=0))
		this->theEigenVector->addVector(0.0, *(another->getEigenVectorPtr()), 1.0);
	else if ((this->theEigenVector ==0) && (another->getEigenVectorPtr()!=0))
		theEigenVector = new Vector(*(another->getEigenVectorPtr()));
	else if ((this->theEigenVector !=0) && (another->getEigenVectorPtr()==0)) {
		delete theEigenVector;
		theEigenVector=0;
	}


//	this->Rotation = another->getRotationPtr();
//	this->theDesignPt= another->theDesignPt;
//	this->theEigenVector = another->getEigenVectorPtr();
//	this->theGFunEvaluator = another->theGFunEvaluator;
//	this->theProbabilityTransformation = another->theProbabilityTransformation;
	return 0;
}

int PrincipalPlane::getNumOfPlane()
{
	return numOfPlane;
}

double PrincipalPlane::getCurvature()
{
	return theCurvature;
}

Matrix * PrincipalPlane::getGridValuesPtr()
{
	return gridValues;
}

Matrix * PrincipalPlane::getRotationPtr()
{
	return Rotation;
}

Vector * PrincipalPlane::getEigenVectorPtr()
{
	return theEigenVector;
}

int PrincipalPlane::setNumOfPlane(int num)
{
	numOfPlane=num;
	return 0;
}

double PrincipalPlane::getValueG2OnGrid(double p_x, double p_y, double valueOfG, double littleDt)
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
			opserr << "PrincipalPlane::getValueG2OnGrid - " << endln
			<< " could not transform from u to x." << endln;
			return -1;
		}
			
		Vector x = theProbabilityTransformation->get_x();
		*/
		Vector x;
		int result = theProbabilityTransformation->transform_u_to_x(u, x);
		if (result < 0) {
		  opserr << "PrincipalPlane::getValueG2OnGrid - " << " could not transform from u to x." << endln;
		  return -1;
		}

		result = theGFunEvaluator->runAnalysis();
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
		gFunctionValue2 = valueOfG; // So it compiles
	}
	else {
	  // This needs to be fixed -- MHS 10/7/2011
	  //gFunctionValue2 = theGFunEvaluator->getG2(valueOfG, littleDt);
	  gFunctionValue2 = valueOfG;
	}

	
	FEConvergence = true;
	return gFunctionValue2;

/*	if (!FEConvergence) return -1.0;
	 
	if ((fabs(p_x) < 1.0e-14) &&(fabs(p_y) < 1.0e-14)) {
		FEConvergence = true;
		Vector u = this->getPtInOrigUSpace(p_x,p_y);
		theProbabilityTransformation->set_u(u);
		int result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "Hessian::formHessianByFDM - " << endln
				<< " could not transform from u to x." << endln;
			return -1;
		}

			
		Vector x = theProbabilityTransformation->get_x();
			


		// ----remove sensitivity algorithm ---

		result = theGFunEvaluator->runGFunAnalysis(x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not run analysis to evaluate limit-state function. " << endln;
			FEConvergence = false;
			exit(-1);
		}
		result = theGFunEvaluator->evaluateG(x);
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not tokenize limit-state function. " << endln;
			exit(-1);
		}
		valueOfG = theGFunEvaluator->getG(); 
	}
	double gFunctionValue2 = theGFunEvaluator->getG2(valueOfG, littleDt);
	FEConvergence = true;
	return gFunctionValue2;  */
	
}

int PrincipalPlane::setGridValueG2(int i, int j, double value)
{
   if (gridValuesG2 ==0) 
		gridValuesG2 = new Matrix(numOfGridXPt, numOfGridYPt);

	(*gridValuesG2)(i,j)= value;
	return 0;
}
