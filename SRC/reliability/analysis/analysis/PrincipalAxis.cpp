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

#include "PrincipalAxis.h"
#include "UniformExperimentalPointRule1D.h"
#include <string.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PrincipalAxis::PrincipalAxis(int tag){

	numOfAxis =tag;  // from 1
	axisDirection=0;
	theExperimentalPointRule=0;
	valuesOnAxis =0;
	curvature =0.0;
	shapeFuncCoeff = 0;	
	tmp1=0;
}

PrincipalAxis::PrincipalAxis(int tag, ExperimentalPointRule1D * pExperimentalPointRule)
{
	numOfAxis = tag;
	axisDirection=0;

	if (pExperimentalPointRule==0) theExperimentalPointRule=0;

	else if (strcmp(pExperimentalPointRule->getType(),"UniformExperimentalPointRule1D")==0){
		theExperimentalPointRule=new UniformExperimentalPointRule1D( pExperimentalPointRule);
	}
	else {
		opserr<<"Fatal: only UniformExperimentalPointRule1D is implemented now. unrecgnized ExperimentalPointRule!";
		exit(-1);
	}

	curvature =0.0;
	valuesOnAxis = 0;
	valuesG2OnAxis =0;
	shapeFuncCoeff = 0;	
	tmp1=0;
}
	
PrincipalAxis::PrincipalAxis(int tag, PrincipalAxis * thePrincipalAxis){

	numOfAxis = tag;

	if (thePrincipalAxis==0){
		opserr<<"Fatal: thePrincipalAxis=0 in PrincipalAxis::PrincipalAxis"; exit(-1);
	}

	this->numOfAxis = thePrincipalAxis->getNumOfAxis();
	this->axisDirection = new Vector(*(thePrincipalAxis->getAxisDirection()));

	

	ExperimentalPointRule1D * tmp = thePrincipalAxis->getExperimentalPointRule();

	if (tmp==0) 	theExperimentalPointRule=0;

	else if (strcmp(tmp->getType(),"UniformExperimentalPointRule1D")==0){
		theExperimentalPointRule=new UniformExperimentalPointRule1D(tmp);
	}
	else {
		opserr<<"Fatal: only UniformExperimentalPointRule1D is implemented now. unrecgnized ExperimentalPointRule!";
		exit(-1);
	}
	curvature =thePrincipalAxis->getCurvature();	
	valuesOnAxis = 0;
	valuesG2OnAxis =0;
	shapeFuncCoeff = 0;
 	tmp1=0;
};

PrincipalAxis::~PrincipalAxis()
{
	if (axisDirection !=0) {
		delete axisDirection;
		axisDirection=0;
	}

	if (theExperimentalPointRule !=0) delete theExperimentalPointRule;
	if (valuesOnAxis !=0) delete valuesOnAxis;
	if (valuesG2OnAxis !=0) delete valuesG2OnAxis;
	
	if (shapeFuncCoeff !=0) {
		if (shapeFuncCoeff[0] !=0){
			for (int i =0; i< this->getExperimentalPointRule()->getNumberOfPoints(); i++){
				delete shapeFuncCoeff[i];
				shapeFuncCoeff[i]=0;
				}
		}
		delete shapeFuncCoeff;
		shapeFuncCoeff =0;
	
	}
	if (tmp1 !=0) delete tmp1;

}

int PrincipalAxis::getNumOfAxis()
{
	return numOfAxis;
}

Vector * PrincipalAxis::getAxisDirection()
{
	return axisDirection;
}

void PrincipalAxis::setNumOfAxis(int i)
{
	numOfAxis=i;
}

void PrincipalAxis::setAxisDirection(Vector *direction)
{
	if (axisDirection!=0) delete axisDirection;
	axisDirection = new Vector(*direction);
}

ExperimentalPointRule1D * PrincipalAxis::getExperimentalPointRule()
{
	return theExperimentalPointRule;
}

void PrincipalAxis::setExperimentalPointRule(ExperimentalPointRule1D * pExperimentalPointRule)
{
	if (strcmp(pExperimentalPointRule->getType(), "UniformExperimentalPointRule1D")==0){
		theExperimentalPointRule = new UniformExperimentalPointRule1D(pExperimentalPointRule);
		if (valuesOnAxis !=0) delete valuesOnAxis;
	    valuesOnAxis = new Vector(theExperimentalPointRule->getNumberOfPoints());
		valuesOnAxis->Zero();
	}
	else{
		opserr<<"unknown type of ExperimentalPointRule";
		exit(-1);
	}
}

void PrincipalAxis::setValuesOnAxis(Vector *pValue)
{
	if (valuesOnAxis!=0) delete valuesOnAxis;
	valuesOnAxis = new Vector(*pValue);	
}


void PrincipalAxis::setValueOnAxis( int i, double value){

	if ( valuesOnAxis ==0) { valuesOnAxis = new Vector(theExperimentalPointRule->getNumberOfPoints()); valuesOnAxis->Zero();}

	(*valuesOnAxis)(i) = value;


};

void PrincipalAxis::setValueG2OnAxis( int i, double value){

	if ( valuesG2OnAxis ==0) { valuesG2OnAxis = new Vector(theExperimentalPointRule->getNumberOfPoints()); valuesG2OnAxis->Zero();}

	(*valuesG2OnAxis)(i) = value;


};


Vector * PrincipalAxis::getValuesOnAxis()
{
	if (valuesOnAxis ==0)
		opserr<<"Warning: PrincipalAxis::getValuesOnAxis(), valuesOnAxis=0\n";
	
	return valuesOnAxis;
}

void PrincipalAxis::cleanValuesOnAxis()
{
	if (valuesOnAxis !=0) delete valuesOnAxis;
	valuesOnAxis =0;
}

void PrincipalAxis::setCurvature(double pCurv)
{
	curvature = pCurv;
}

double PrincipalAxis::getCurvature()
{
	return curvature;
}

int PrincipalAxis::copyValues(PrincipalAxis *another)
{
	this->curvature = another->getCurvature();
	this->setAxisDirection( another->getAxisDirection()); 
	cleanValuesOnAxis();
	return 0;
}



double PrincipalAxis::getValueOnAxis(int i)
{
	if (valuesOnAxis !=0)
	    return (*valuesOnAxis)(i);
	else{
		opserr<<"Fatal: PrincipalAxis::getValueOnAxis() not axis value exist\n";
		exit(-1);
	}
}


double PrincipalAxis::getValueG2OnAxis(int i)
{
	if (valuesG2OnAxis !=0)
	    return (*valuesG2OnAxis)(i);
	else{
		opserr<<"Fatal: PrincipalAxis::getValueG2OnAxis() not axis value exist\n";
		exit(-1);
	}
}

int PrincipalAxis::computeShapeFuncCoeff()
{
	int numOfPoints = this->getExperimentalPointRule()->getNumberOfPoints();

	if (shapeFuncCoeff !=0) {
		if (shapeFuncCoeff[0] !=0){
			for (int i =0; i< numOfPoints; i++){
				delete shapeFuncCoeff[i];
				shapeFuncCoeff[i]=0;
			}
		}
		delete shapeFuncCoeff;
		shapeFuncCoeff =0;

	}
	
// -- alloc	mem
	shapeFuncCoeff =new Vector * [numOfPoints];
	for (int i =0; i< numOfPoints; i++){
		shapeFuncCoeff[i] =  new Vector(numOfPoints);
		shapeFuncCoeff[i]->Zero();
	}		 

	
// -- compute coeff

	// refer paper by Rahman. "Decomposition methods for structural reliability analysis"
	/* Shape Function at point j is coeff of, Fj = (x-x[1])*(x-x[2])*...*(x-x[j-1])*(x-x[j+1])...(x-x[n])/Denominator
  	   where Denominator = (x[j]-x[1])*(x[j]-x[2])*...*(x[j]-x[j-1])*(x[j]-x[j+1])...(x[j]-x[n])
	   rearrange it to get   Fj = (a[0]*x^(n-1)+a[1]*x^(n-2)+...+a[n-1]*x^0) / Denominator
	   *shapeFuncCoeff[j] = {a[0],  a[1] ...a[n-1]}/Denominator    Note a[0]=1.0    */
	//                      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
	int size = numOfPoints;
	Vector temp(size-1);
	temp.Zero();
	Vector x( *(this->getExperimentalPointRule()->getPointCoordinates()));
	
	double denominator;

	for (int jj =0; jj<size; jj++){ // for each point jj on axis

	  int kk;
	  for(kk=0; kk<jj; kk++) temp(kk) = x(kk);
	  for( kk=jj+1; kk<size; kk++) temp(kk-1) = x(kk);
	  
	  denominator = 1.0;
	  
	  for ( kk=0; kk<size; kk++){
	    if (kk !=jj) denominator *= (x(jj)-x(kk));
	  }
	
//			double theValue = thePrincipalAxes[ii]->getValueOnAxis(jj);


//			opserr<<"theValue:"<<theValue<<"    denom:"<<denominator<<endln;
//			Vector coeff(*Poly(&temp));

		shapeFuncCoeff[jj]->addVector(0.0, *Poly(&temp), 1.0/denominator);



	} // for each point jj on this axis
	

	return 0;
}

Vector * PrincipalAxis::getShapeFuncCoeff(int j)  //j is point
{
	if (shapeFuncCoeff ==0) {
		opserr<<"warning: shapeFuncCoeff is not computed yet, recomputing..."<<endln;
		this->computeShapeFuncCoeff();
	
	}
	return shapeFuncCoeff[j];
}



Vector * PrincipalAxis::Poly(Vector *x)
{


/*	% Expand recursion formula
for j=1:n
    
    for k=2:j+1
        c(k) = c1(k) - e(j)*c1(k-1);
    end
    
     c1 = c;
end

*/

	int n = x->Size();

	if (tmp1 ==0) {
		tmp1 = new Vector(n+1);
	}
	else if (tmp1->Size() != n+1){
		delete tmp1;
		tmp1 = new Vector(n+1);
	}

		
	tmp1->Zero();
	(*tmp1)(0)=1.0;

	Vector  tmp_old(*tmp1);

	for (int j=0; j<=n-1; j++){
		for (int k =1; k<=j+1; k++){  
			(*tmp1)(k) = tmp_old(k) - (*x)(j)*tmp_old(k-1);
		}
		tmp_old.addVector(0.0, *tmp1,1.0);
	}

	
	return tmp1;
}
