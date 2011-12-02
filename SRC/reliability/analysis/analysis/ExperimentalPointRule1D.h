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
 

#if !defined AFX_EXPERIMENTALPOINTRULE_H
#define AFX_EXPERIMENTALPOINTRULE_H
#include <Vector.h>

class ExperimentalPointRule1D  
{
public:
	virtual int getPointClosestToOrigin();
	virtual Vector * getPointCoordinates();
	virtual void setInfo( int n, double begin, double end);
	virtual double  getEndOfGrid( ){return end;};
	virtual double getBeginOfGrid( ){return begin;};
	virtual double getPointCoordinate(int i)=0;
	virtual char * getType()=0;
	virtual int getNumberOfPoints(){return num;};

	virtual void setBeginOfGrid(double  x1){begin=x1;};
	virtual void setEndOfGrid(double x2){end=x2;};
	virtual void setNumberOfPoints(int pNum){num=pNum;};
	


	ExperimentalPointRule1D ();
	virtual ~ExperimentalPointRule1D ();

protected:

	double begin;
	double end;
	int num;
	
};

#endif 
