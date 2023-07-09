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
 

#if !defined AFX_GRIDPLANE_H 
#define AFX_GRIDPLANE_H 

#include <ReliabilityDomainComponent.h>

#include <FunctionEvaluator.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>
#include <ExperimentalPointRule1D.h>
#include <PrincipalAxis.h>

class GridPlane  
{
public:

	int setGridValueG2(int i, int j, double v);
	Matrix * getGridValuesG2Ptr();
	double getValueG2OnGrid(double p_x, double p_y, double gFunValue, double littleDt);


	GridPlane(int pTag, Vector * pDesignPt, 
							   PrincipalAxis * firstAxis, PrincipalAxis * secondAxis, Matrix * pRotation, 
							   ProbabilityTransformation * pProbabilityTransformation,
							   FunctionEvaluator * pGFunEvaluator );

	virtual ~GridPlane();

	int getNumOfPlane();
	int getnumOfAxis(int i);
	PrincipalAxis * getAxisPtr(int i);
	double getValueOnGrid( Vector & point);  // compute the grid value
	Vector getPointCoordinate( int i, int j); // compute the grid value

	double getEndOfGridX();
	double getBeginOfGridX();
	int getNumOfGridXPt();
	double getEndOfGridY();
	double getBeginOfGridY();
	int getNumOfGridYPt();

	Matrix * getRotationPtr();
	Matrix * getGridValuesPtr();
    
	Vector getPtInOrigUSpace( double x, double y);
	double getValueOnGrid(double x, double y); 

	int setNumOfPlane( int num);
	int setAxis( int pAxis, PrincipalAxis * theAxis);
	
	int setRotation( Matrix * R);
	int setDesignPoint (Vector * pDesignPt);
	int setGridValue( int i, int j, double v);



	//int copyValues( GridPlane * anotherPlane);
	int cleanGridMatrix();

	bool isFEConverged();

	double getSavedValueOnGrid(int i, int j);  // does not do computation
	double getSavedValueG2OnGrid(int i, int j);  // does not do computation


private:
	Matrix * gridValuesG2;
	bool FEConvergence;
	Matrix * gridValues;

	FunctionEvaluator * theGFunEvaluator;
	ProbabilityTransformation * theProbabilityTransformation;
	Vector * theDesignPt; //Uspace



    Matrix * Rotation ;   // u''= R*u 


	int numOfPlane;
//	double theCurvature;
	PrincipalAxis * theFirstAxis; // the x axils in original U space
    PrincipalAxis * theSecondAxis; // the y axils in original U space




};

#endif // 
