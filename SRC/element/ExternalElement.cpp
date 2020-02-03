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

#include <ExternalElement.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <Parameter.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

ExternalElement::ExternalElement(int tag)
  :Element(tag,ELE_TAG_ExternalElement)
{
  
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
ExternalElement::ExternalElement()
  :Element(0,ELE_TAG_ExternalElement)
{
  
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
ExternalElement::~ExternalElement()
{
	// free function pointer here
	if (_Ele_GetClassType != 0) free(_Ele_GetClassType);
	if (_Ele_GetNumExternalNodes != 0) free(_Ele_GetNumExternalNodes);
	if (_Ele_GetNodePtrs != 0) free(_Ele_GetNodePtrs);
	if (_Ele_GetExternalNodes != 0) free(_Ele_GetExternalNodes);
	if (_Ele_GetNumDOF != 0) free(_Ele_GetNumDOF);
	if (_Ele_SetDomain != 0) free(_Ele_SetDomain);
	if (_Ele_CommitState != 0) free(_Ele_CommitState);
	if (_Ele_RevertToLastCommit != 0) free(_Ele_RevertToLastCommit);
	if (_Ele_RevertToStart != 0) free(_Ele_RevertToStart);
	if (_Ele_Update != 0) free(_Ele_Update);
	if (_Ele_GetTangentStiff != 0) free(_Ele_GetTangentStiff);
	if (_Ele_GetInitialStiff != 0) free(_Ele_GetInitialStiff);
	if (_Ele_GetDamp != 0) free(_Ele_GetDamp);
	if (_Ele_GetMass != 0) free(_Ele_GetMass);
	if (_Ele_ZeroLoad != 0) free(_Ele_ZeroLoad);
	if (_Ele_AddLoad != 0) free(_Ele_AddLoad);
	if (_Ele_AddInertiaLoadToUnbalance != 0) free(_Ele_AddInertiaLoadToUnbalance);
	if (_Ele_GetResistingForce != 0) free(_Ele_GetResistingForce);
	if (_Ele_GetResistingForceIncInertia != 0) free(_Ele_GetResistingForceIncInertia);
	if (_Ele_Print != 0) free(_Ele_Print);
	if (_Ele_SetResponse != 0) free(_Ele_SetResponse);
	if (_Ele_GetResponse != 0) free(_Ele_GetResponse);
}

int
ExternalElement::getNumExternalNodes(void) const
{
    return  _Ele_GetNumExternalNodes();
}

const ID &
ExternalElement::getExternalNodes(void) 
{
	return *_Ele_GetExternalNodes();
}

Node **
ExternalElement::getNodePtrs(void) 
{
  return _Ele_GetNodePtrs();
}

int
ExternalElement::getNumDOF(void) 
{
	return _Ele_GetNumDOF();	
}

void
ExternalElement::setDomain(Domain *theDomain)
{
	_Ele_SetDomain(theDomain);
}

int
ExternalElement::commitState()
{
  return _Ele_CommitState();
}

int
ExternalElement::revertToLastCommit()
{
	return _Ele_RevertToLastCommit();
}

int
ExternalElement::revertToStart()
{
	return _Ele_RevertToStart();
}

int
ExternalElement::update(void)
{
  return _Ele_Update();
}

const Matrix &
ExternalElement::getTangentStiff(void)
{
    return *_Ele_GetTangentStiff();
}


const Matrix &
ExternalElement::getInitialStiff(void)
{
    return *_Ele_GetInitialStiff();
}


const Matrix &
ExternalElement::getDamp(void)
{
    return *_Ele_GetDamp();
}


const Matrix &
ExternalElement::getMass(void)
{
    return *_Ele_GetMass();
}


void 
ExternalElement::zeroLoad(void)
{
	_Ele_ZeroLoad();
}

int 
ExternalElement::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  return _Ele_AddLoad(theLoad, loadFactor);
}



int 
ExternalElement::addInertiaLoadToUnbalance(const Vector &accel)
{
    return _Ele_AddInertiaLoadToUnbalance(accel);
}

const Vector &
ExternalElement::getResistingForce()
{
    return *_Ele_GetResistingForce();
}

const Vector &
ExternalElement::getResistingForceIncInertia()
{	
    return *_Ele_GetResistingForceIncInertia();
}

int
ExternalElement::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
ExternalElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

int
ExternalElement::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}

void
ExternalElement::Print(OPS_Stream &s, int flag)
{
	_Ele_Print(flag);
}

Response*
ExternalElement::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	return _Ele_SetResponse(argv, argc, output);
}

int 
ExternalElement::getResponse(int responseID, Information &eleInfo)
{
    return _Ele_GetResponse(responseID,eleInfo);
}

int
ExternalElement::setParameter(const char **argv, int argc, Parameter &param)
{
  return 0;
    
}

int
ExternalElement::updateParameter(int parameterID, Information &info)
{
  return 0;
    
}
