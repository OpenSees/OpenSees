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
 

#if !defined(AFX_UnivariateDecomposition_H__A0BE5A2C_B670_47A8_8AC9_435AEDA27255__INCLUDED_)
#define AFX_UnivariateDecomposition_H__A0BE5A2C_B670_47A8_8AC9_435AEDA27255__INCLUDED_

 
#include "SurfaceDesign.h"
#include "PrincipalAxis.h"

class UnivariateDecomposition : public SurfaceDesign  
{
public:

	UnivariateDecomposition(int numAxis, /*Vector * gradient,*/ PrincipalAxis ** pAxis, bool isTimeVariant = false);
	virtual ~UnivariateDecomposition();
 
	void setPincipalAxesPtr(PrincipalAxis ** pAxis);
//	void setGradient( Vector * gradient);
//	Vector * Poly( Vector * x);

 
	int fitCurve() ;
    double getFunctionValue(Vector * Point);
	double getFunctionValue2(Vector * point, Vector * dp2prime, Vector * gradG2);
	char * getType();


private:
	double linearCorrection;
	bool isTimeVariant;
	Vector * tmp;
	Vector ** coefficients;
	Vector ** coefficients_2;  // second LLS
	PrincipalAxis ** thePrincipalAxes;
    char type[40];
//	Vector * gradientDP;  // U space

	int numAxes;   // = numOfPrincipalAxis +1, since dp direction is accounted. 
};

#endif // !defined(AFX_UnivariateDecomposition_H__A0BE5A2C_B670_47A8_8AC9_435AEDA27255__INCLUDED_)
