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
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2008-05-22 20:05:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ParameterPositioner.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ParameterPositioner_h
#define ParameterPositioner_h

#include <ReliabilityDomainComponent.h>
#include <Parameter.h>

class DomainComponent;

class ParameterPositioner : public ReliabilityDomainComponent
{

public:

	ParameterPositioner(int tag, DomainComponent *theObject,
			    const char **argv, int argc);
	ParameterPositioner(int tag, Parameter &param);
	~ParameterPositioner();

	void Print(OPS_Stream &s, int flag =0);

	int update(int newValue); 
	int update(double newValue); 
	int activate(bool active);

	void setGradNumber(int gradNum) {gradNumber = gradNum;}
	int getGradNumber(void) {return gradNumber;}

protected:

private:
	int gradNumber; // 0,...,nparam-1

	Parameter theParameter;
};

#endif
