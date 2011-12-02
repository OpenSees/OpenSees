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
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "UniformExperimentalPointRule1D.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

UniformExperimentalPointRule1D::UniformExperimentalPointRule1D()
{
	strcpy(type, "UniformExperimentalPointRule1D");
	// --- default
	begin =-1.0;
	end =1.0;
	num= 11;
	tmp=0;

}

UniformExperimentalPointRule1D::UniformExperimentalPointRule1D(ExperimentalPointRule1D * pExperimentalPointRule1D){

	if (pExperimentalPointRule1D ==0) {
		opserr<<"Fatal: UniformExperimentalPointRule1D::UniformExperimentalPointRule1D(ExperimentalPointRule1D * pExperimentalPointRule1D), pExperimentalPointRule1D=0.";
		exit(-1);
	
	}
	begin = pExperimentalPointRule1D->getBeginOfGrid();
	end = pExperimentalPointRule1D->getEndOfGrid();
	this->num = pExperimentalPointRule1D->getNumberOfPoints();
	strcpy(type, "UniformExperimentalPointRule1D");
	tmp=0;

};

UniformExperimentalPointRule1D::~UniformExperimentalPointRule1D()
{

	if (tmp!=0) {delete tmp; tmp=0;} 
}


double UniformExperimentalPointRule1D::getPointCoordinate(int i){

	if (i>=num) {opserr<<"error in UniformExperimentalPointRule1D::getPointCoord, i is:"<<i<<endln; exit(-1);}
    double tmp= begin+(end-begin)/(num-1)*i;
	return tmp;

};

Vector * UniformExperimentalPointRule1D::getPointCoordinates(){
	
	if (tmp==0) {tmp = new Vector(num);}
	else if (tmp->Size() != num) {delete tmp; tmp = new Vector(num);}

	for( int i =0; i<num; i++)
		(*tmp)(i) = getPointCoordinate(i);
	return tmp;


}

char * UniformExperimentalPointRule1D::getType(){

	return type;


};

int UniformExperimentalPointRule1D::getPointClosestToOrigin()
{
	Vector * temp = this->getPointCoordinates();
	int k=0;
	double value = 1.0;
	for(int i =0; i< num; i++){
		if (fabs((*temp)(i))<value){
			value = fabs((*temp)(i));
			k=i;
		}
	}

	return k;
}
