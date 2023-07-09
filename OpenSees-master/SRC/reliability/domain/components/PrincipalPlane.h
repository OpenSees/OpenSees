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

#if !defined(GuConte_PRINCIPALPLANE)
#define GuConte_PRINCIPALPLANE

#include <ReliabilityDomainComponent.h>

#include <FunctionEvaluator.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
class PrincipalPlane 
{
public:
	PrincipalPlane(int pTag, Vector * pDesignPt, 
		       Vector * pEigenVector, Matrix * pRotation, 
		       ProbabilityTransformation * pProbabilityTransformation,
		       FunctionEvaluator *, double);
	~PrincipalPlane();

	int setGridValueG2(int i, int j, double value);
	double getValueG2OnGrid(double x,double y,double  valueOfG,double  littleDt);
	int setNumOfPlane( int num);
	Vector * getEigenVectorPtr();
	Matrix * getRotationPtr();
	Matrix * getGridValuesPtr();
	double getCurvature();
	int getNumOfPlane();
	int copyValues( PrincipalPlane * anotherPPlane);
	int cleanGridMatrix();
	int setGridValue(int i, int j, double value);
	double getEndOfGridX();
	double getBeginOfGridX();
	int getNumOfGridXPt();

	double getEndOfGridY();
	double getBeginOfGridY();
	int getNumOfGridYPt();
	
	int setCurvature( double pCurv);
	int setEigenVector( Vector * pEigenV);
	int setRotation(Matrix * R);
	int setDesignPoint (Vector * pDesignPt);
	int setGridYInfo( int numOfGrid, double begin, double end);
	int setGridXInfo( int numOfGrid, double begin, double end);
	Vector getPtInOrigUSpace( double x, double y);
	double getValueOnGrid(double x, double y); 

private:
	bool FEConvergence;
	Matrix * gridValues;
	Matrix * gridValuesG2;
	double theCurvature;
	FunctionEvaluator * theGFunEvaluator;
	ProbabilityTransformation * theProbabilityTransformation;
	Vector * theDesignPt;
	int numOfPlane;

	int numOfGridXPt;
	double endOfGridX;
	double beginOfGridX;

	int numOfGridYPt;
	double endOfGridY;
	double beginOfGridY;

    Matrix * Rotation ;   // u''= R*u 
	Vector * theEigenVector; // the x axils in original U space


};

#endif 
