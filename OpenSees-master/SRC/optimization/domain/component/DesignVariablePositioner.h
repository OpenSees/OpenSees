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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//


#ifndef DesignVariablePositioner_h
#define DesignVariablePositioner_h

#include <Parameter.h>
#include <OptimizationDomain.h>

#include <Information.h>
#include <DomainComponent.h>

#include <OptimizationDomainComponent.h>

class DesignVariablePositioner : public OptimizationDomainComponent
{

public:

	DesignVariablePositioner(int tag, OptimizationDomain  * theOptimizationDomain, int DVnumber, DomainComponent *theObject, const char **argv, int argc);
	DesignVariablePositioner(int tag, OptimizationDomain  * theOptimizationDomain, int DVnumber, Parameter &param);
	~DesignVariablePositioner();

	void Print(OPS_Stream &s, int flag =0);

	int getDVNumber(void);
	int update (double newValue); 
 

	int setNewTag(int newTag);
	int setDVNumber(int newRvNumber);

	int sendSelf(int,Channel &) {return 0;}
	int recvSelf(int,Channel &,FEM_ObjectBroker &) {return 0;}


protected:

private:
	int dVNumber;

	Parameter theParameter;
};

#endif






