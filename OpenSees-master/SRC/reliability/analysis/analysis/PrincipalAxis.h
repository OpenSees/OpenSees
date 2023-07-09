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

#if !defined(AFX_PRINCIPALAXIS_H)
#define AFX_PRINCIPALAXIS_H
#include "Vector.h"
#include "ExperimentalPointRule1D.h"

class PrincipalAxis  
{
public:
	Vector * Poly(Vector *x);
	Vector * getShapeFuncCoeff( int j);
	int computeShapeFuncCoeff();

	PrincipalAxis( int tag );
	PrincipalAxis(int tag,ExperimentalPointRule1D * ptheExperimentalPointRule);
	PrincipalAxis(int tag,PrincipalAxis * thePrincipalAxis);
	virtual ~PrincipalAxis();

	

	double getValueOnAxis(int i);
	double getValueG2OnAxis(int i);
	double getCurvature();
	void cleanValuesOnAxis();
	Vector * getValuesOnAxis();
	ExperimentalPointRule1D * getExperimentalPointRule();
	int getNumOfAxis();
 
	int copyValues( PrincipalAxis * another);  // copy only curvature and eigenVector, not others( expPtRule, numOfAxis)


	void setCurvature( double pCurv);
	void setValueOnAxis( int i, double value);
	void setValueG2OnAxis( int i, double value);
	void setValuesOnAxis( Vector * value);
	void setExperimentalPointRule(ExperimentalPointRule1D * theExpPtRule);

	void setAxisDirection( Vector * direction);
	void setNumOfAxis(int i);
	Vector * getAxisDirection();



private:
	Vector ** shapeFuncCoeff;   
		
// -- shapeFuncCoeff

	// refer paper by Rahman. "Decomposition methods for structural reliability analysis"
	/* Shape Function at point j is coeff of, Fj = (x-x[1])*(x-x[2])*...*(x-x[j-1])*(x-x[j+1])...(x-x[n])/Denominator
  	   where Denominator = (x[j]-x[1])*(x[j]-x[2])*...*(x[j]-x[j-1])*(x[j]-x[j+1])...(x[j]-x[n])

	   rearrange it to get   Fj = {a[0]*x^(n-1)+a[1]*x^(n-2)+...+a[n-1]*x^0} / Denominator
	   *shapeFuncCoeff[j] = {a[0],  a[1] ...a[n-1]}/Denominator    Note a[0]=1.0    */
	//                      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!



	Vector * valuesG2OnAxis;
 
	Vector * tmp1;
	double curvature;
	Vector * valuesOnAxis;
	int numOfAxis;
	Vector * axisDirection;
	ExperimentalPointRule1D * theExperimentalPointRule;
};

#endif // !defined(AFX_PRINCIPALAXIS_H)
