/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 6049 $
// $Date: 2019-01-29 $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/truss/ExternalElement.h $


// Written: M. Salehi 
// Created: Jan 2019
//

#ifndef ExternalElement_h
#define ExternalElement_h



#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;
class UniaxialMaterial;
typedef const char *(__stdcall *Ele_GetClassType) (void);
typedef int (__stdcall *Ele_GetNumExternalNodes) (void);
typedef ID* (__stdcall *Ele_GetExternalNodes) (void);
typedef Node **(__stdcall *Ele_GetNodePtrs) (void);
typedef int(__stdcall *Ele_GetNumDOF) (void);
typedef void(__stdcall *Ele_SetDomain) (Domain *theDomain);
typedef int(__stdcall *Ele_CommitState) ();
typedef int(__stdcall *Ele_RevertToLastCommit) ();
typedef int(__stdcall *Ele_RevertToStart) ();
typedef int(__stdcall *Ele_Update) ();
typedef Matrix *(__stdcall *Ele_GetTangentStiff) (void);
typedef Matrix *(__stdcall *Ele_GetInitialStiff) (void);
typedef Matrix *(__stdcall *Ele_GetDamp) (void);
typedef Matrix *(__stdcall *Ele_GetMass) (void);
typedef void(__stdcall *Ele_ZeroLoad) (void);
typedef int(__stdcall *Ele_AddLoad) (ElementalLoad *theLoad, double loadFactor);
typedef int(__stdcall *Ele_AddInertiaLoadToUnbalance) (const Vector &accel);
typedef Vector *(__stdcall *Ele_GetResistingForce) (void);
typedef Vector *(__stdcall *Ele_GetResistingForceIncInertia) (void);
typedef void(__stdcall *Ele_Print) (int flag);
typedef Response *(__stdcall *Ele_SetResponse) (const char **argv, int argc, OPS_Stream &output);
typedef int(__stdcall *Ele_GetResponse) (int responseID, Information &eleInfo);

class ExternalElement : public Element
{
  public:
    ExternalElement(int tag);
    ExternalElement();    
    ~ExternalElement();

	const char *getClassType(void) const { return _Ele_GetClassType(); };

    // public methods to obtain inforrmation about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
	void SetLinks(Ele_GetClassType GetClassType,
		Ele_GetNumExternalNodes GetNumExternalNodes,
		Ele_GetNodePtrs GetNodePtrs,
		Ele_GetExternalNodes GetExternalNodes,
		Ele_GetNumDOF GetNumDOF,
		Ele_SetDomain SetDomain,
		Ele_CommitState CommitState,
		Ele_RevertToLastCommit RevertToLastCommit,
		Ele_RevertToStart RevertToStart,
		Ele_Update Update,
		Ele_GetTangentStiff GetTangentStiff,
		Ele_GetInitialStiff GetInitialStiff,
		Ele_GetDamp GetDamp,
		Ele_GetMass GetMass,
		Ele_ZeroLoad ZeroLoad,
		Ele_AddLoad AddLoad,
		Ele_AddInertiaLoadToUnbalance AddInertiaLoadToUnbalance,
		Ele_GetResistingForce GetResistingForce,
		Ele_GetResistingForceIncInertia GetResistingForceIncInertia,
		Ele_Print Print,
		Ele_GetResponse GetResponse,
		Ele_SetResponse SetResponse)
	{
		this->_Ele_GetClassType = GetClassType;
		this->_Ele_GetNumExternalNodes = GetNumExternalNodes;
		this->_Ele_GetNodePtrs = GetNodePtrs;
		this->_Ele_GetExternalNodes = GetExternalNodes;
		this->_Ele_GetNumDOF = GetNumDOF;
		this->_Ele_SetDomain = SetDomain;
		this->_Ele_CommitState = CommitState;
		this->_Ele_RevertToLastCommit = RevertToLastCommit;
		this->_Ele_RevertToStart = RevertToStart;
		this->_Ele_Update = Update;
		this->_Ele_GetTangentStiff = GetTangentStiff;
		this->_Ele_GetInitialStiff = GetInitialStiff;
		this->_Ele_GetDamp = GetDamp;
		this->_Ele_GetMass = GetMass;
		this->_Ele_ZeroLoad = ZeroLoad;
		this->_Ele_AddLoad = AddLoad;
		this->_Ele_AddInertiaLoadToUnbalance = AddInertiaLoadToUnbalance;
		this->_Ele_GetResistingForce = GetResistingForce;
		this->_Ele_GetResistingForceIncInertia = GetResistingForceIncInertia;
		this->_Ele_Print = Print;
		this->_Ele_GetResponse = GetResponse;
		this->_Ele_SetResponse = SetResponse;
	}
  protected:
    
  private:
	  Ele_GetClassType _Ele_GetClassType;
	  Ele_GetNumExternalNodes _Ele_GetNumExternalNodes;
	  Ele_GetNodePtrs _Ele_GetNodePtrs;
	  Ele_GetExternalNodes _Ele_GetExternalNodes;
	  Ele_GetNumDOF _Ele_GetNumDOF;
	  Ele_SetDomain _Ele_SetDomain;
	  Ele_CommitState _Ele_CommitState;
	  Ele_RevertToLastCommit _Ele_RevertToLastCommit;
	  Ele_RevertToStart _Ele_RevertToStart;
	  Ele_Update _Ele_Update;
	  Ele_GetTangentStiff _Ele_GetTangentStiff;
	  Ele_GetInitialStiff _Ele_GetInitialStiff;
	  Ele_GetDamp _Ele_GetDamp;
	  Ele_GetMass _Ele_GetMass;
	  Ele_ZeroLoad _Ele_ZeroLoad;
	  Ele_AddLoad _Ele_AddLoad;
	  Ele_AddInertiaLoadToUnbalance _Ele_AddInertiaLoadToUnbalance;
	  Ele_GetResistingForce _Ele_GetResistingForce;
	  Ele_GetResistingForceIncInertia _Ele_GetResistingForceIncInertia;
	  Ele_Print _Ele_Print;
	  Ele_SetResponse _Ele_SetResponse;
	  Ele_GetResponse _Ele_GetResponse;
};

#endif
