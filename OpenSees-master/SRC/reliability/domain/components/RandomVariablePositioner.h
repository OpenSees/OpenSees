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
                                                                        
// $Revision: 1.11 $
// $Date: 2010-06-11 15:56:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositioner.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef RandomVariablePositioner_h
#define RandomVariablePositioner_h

#include <ReliabilityDomainComponent.h>
#include <Parameter.h>
class RandomVariable;
class DomainComponent;

class RandomVariablePositioner : public ReliabilityDomainComponent
{

public:

	RandomVariablePositioner(int tag, int RVindex, DomainComponent *theObject,
						const char **argv, int argc,
						RandomVariable *theRVptr);
	RandomVariablePositioner(int tag, int RVindex, Parameter &param, RandomVariable *theRV);
	~RandomVariablePositioner();

	void Print(OPS_Stream &s, int flag =0);

	int update (double newValue); 
	int activate(bool active);

	int setNewTag(int newTag); // Do we really need this anymore? -- MHS
	int setRvIndex(int newRvIndex);
	int getRvIndex(void);
	RandomVariable *getRandomVariable(void) const {return theRV;}

	int getRVTag(void);
	int getParamTag(void);

protected:

private:
	int rvIndex;

	RandomVariable *theRV;

	Parameter theParameter;
};

#endif

