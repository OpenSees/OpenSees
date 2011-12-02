// DesignVariablePositioner.h: interface for the DesignVariablePositioner class.
//
//////////////////////////////////////////////////////////////////////

// ----- Quan Gu Aug 2005 for optimization interface ---

#ifndef DesignVariablePositioner_h
#define DesignVariablePositioner_h

#include <ReliabilityDomainComponent.h>
#include <ReliabilityDomain.h>

#include <Parameter.h>
#include <DomainComponent.h>

// Quan Gu
// this class need to be extended to accomodate the DV is itself a RV case (in RBOD )later. 

class DesignVariablePositioner : public ReliabilityDomainComponent
{

public:

	DesignVariablePositioner(int tag, ReliabilityDomain  * theReliabilityDomain, int DVnumber, DomainComponent *theObject, const char **argv, int argc);
	~DesignVariablePositioner();

	void Print(OPS_Stream &s, int flag =0);

	int getDVNumber(void);
	int update (double newValue); 
//	int activate(bool active);

	int setNewTag(int newTag);
	int setDVNumber(int newRvNumber);

protected:

private:
	int dVNumber;

	Parameter theParameter;
	DomainComponent *theObject;
};

#endif






