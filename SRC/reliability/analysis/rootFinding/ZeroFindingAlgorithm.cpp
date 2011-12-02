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
 

#include "ZeroFindingAlgorithm.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ZeroFindingAlgorithm::ZeroFindingAlgorithm()
{
	theSamplingAnalysis=0;
	functionTol= 1.e-14;
	variableTol=1.e-7;

	a=0;
	b=0;
	zeroPoint=0;
	functionType=-1;
	iterationNum =0;
	maxIterNum = 20;


	ii_1=0;
	ii_2=0;
	x_1=0;
	x_2=0;
	G_1=0;
	G_2=0;
}

ZeroFindingAlgorithm::~ZeroFindingAlgorithm()
{

}



void ZeroFindingAlgorithm::setA(double A){
	a=A;
};

void ZeroFindingAlgorithm::setB(double B){
	b=B;
};

double ZeroFindingAlgorithm::getZeroPoint()
{
	return zeroPoint;
}



	

int ZeroFindingAlgorithm::getIterationNumber()
{
	return iterationNum;
}

void ZeroFindingAlgorithm::setMaxIterNum(int newMaxIterNum)
{
	maxIterNum = newMaxIterNum;
}

void ZeroFindingAlgorithm::setFa(double pFa)
{
	Fa=pFa;
}

void ZeroFindingAlgorithm::setFb(double pFb)
{
	Fb=pFb;
}

double ZeroFindingAlgorithm::getFa()
{
	return Fa;
}

double ZeroFindingAlgorithm::getFb()
{
	return Fb;
}

double ZeroFindingAlgorithm::getA()
{
	return a;
}

double ZeroFindingAlgorithm::getB()
{
	return b;
}

void ZeroFindingAlgorithm::setFEconvergence(bool a)
{
	FEconvergence = a;
}


bool ZeroFindingAlgorithm::getFEconvergence(void){
	return FEconvergence;
}

void ZeroFindingAlgorithm::saveXG1(double pX1, double pG1)
{
	if ((x_1==0)||(G_1==0)) {
		opserr<<"ZeroFindingAlgorithm::saveXG1 wrong"<<endln;
		if (x_1==0) x_1 = new double[100];
		if (G_1==0) G_1 = new double[100];
		ii_1=0;
	}
	
	x_1[ii_1]=pX1;
	G_1[ii_1]=pG1;
	ii_1++;
	
	

}

void ZeroFindingAlgorithm::saveXG2(double pX2, double pG2)
{
	if ((x_2==0)||(G_2==0)) {
		opserr<<"ZeroFindingAlgorithm::saveXG2 wrong"<<endln;
		if (x_2==0) x_2 = new double[100];
		if (G_2==0) G_2 = new double[100];
		ii_2=0;
		
	 
	}
	
	x_2[ii_2]=pX2;
	G_2[ii_2]=pG2;
	ii_2++;
	
	
}

double ZeroFindingAlgorithm::getX2(int i)
{
	if(i>=ii_2) { opserr<<"wrong in ZeroFindingAlgorithm::getX2.  i>=ii_2 "<<endln;return -999;};
	return x_2[i];
}

double ZeroFindingAlgorithm::getG2(int i)
{
	if(i>=ii_2) { opserr<<"wrong in ZeroFindingAlgorithm::getG2.  i>=ii_2 "<<endln; return -999;};
	return G_2[i];
}

void ZeroFindingAlgorithm::setX1Pointer(double * pX1)
{
	this->x_1 = pX1;
}

void ZeroFindingAlgorithm::setG1Pointer(double * pG1)
{
	this->G_1 = pG1;
}


void ZeroFindingAlgorithm::setX2Pointer(double * pX2)
{
	this->x_2 = pX2;
}

void ZeroFindingAlgorithm::setG2Pointer(double * pG2)
{
	this->G_2 = pG2;
}

void ZeroFindingAlgorithm::set_ii_1(int pI)
{
	this->ii_1=pI;
}

void ZeroFindingAlgorithm::set_ii_2(int pI)
{
	this->ii_2=pI;
}

int ZeroFindingAlgorithm::get_ii_1()
{
	return ii_1;
}

int ZeroFindingAlgorithm::get_ii_2()
{
	return ii_2;
}

void ZeroFindingAlgorithm::setFunctionType(int i)
{
	this->functionType = i;
}

int ZeroFindingAlgorithm::getFunctionType()
{
	return this->functionType;
}
