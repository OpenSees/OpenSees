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
 


#if !defined(AFX_BivariateDecomposition_H__D4014C3F_A17E_4564_A385_0EABA6D3D457__INCLUDED_)
#define AFX_BivariateDecomposition_H__D4014C3F_A17E_4564_A385_0EABA6D3D457__INCLUDED_

#include "SurfaceDesign.h"
#include "PrincipalAxis.h"
#include "GridPlane.h"

class BivariateDecomposition : public SurfaceDesign  
{
public:
	BivariateDecomposition(int numAxis, PrincipalAxis ** pAxis, GridPlane ** thePlane, bool isTimeVariant = false);
	virtual ~BivariateDecomposition();

 
 
public:

	void setGridPlanesPtr(GridPlane ** pPlane);
	Vector ** axisCoeff_2;


  
	void setPincipalAxesPtr(PrincipalAxis ** pAxis);
	//void setGradient( Vector * gradient);
	Vector * Poly( Vector * x);

 
	int fitCurve() ;
    double getFunctionValue(Vector * Point);
	double getFunctionValue2(Vector * point, Vector * dp2prime, Vector * gradG2);
	double debug(Vector * point, Vector * dp2prime, Vector * gradG2);
	char * getType();


private:
	double linearCorrection;
	Vector ** axisCoeff;
	GridPlane ** theGridPlanes;
	PrincipalAxis ** thePrincipalAxes;

	bool isTimeVariant;
	Vector * tmp1;
	Matrix ** coefficients;
	Matrix ** coefficients_2;  // second LLS

    char type[40];
	//Vector * gradientDP;  // U space

	int numAxes;   // = numOfPrincipalAxis +1, since dp direction is accounted. 

};

#endif // !defined(AFX_BivariateDecomposition_H__D4014C3F_A17E_4564_A385_0EABA6D3D457__INCLUDED_)


 
