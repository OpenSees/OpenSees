#include "stdafx.h"

#include "DomainWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Components;
using namespace OpenSees::Components::LoadPatterns;
using namespace OpenSees::Components::Constraints;
using namespace OpenSees::Components::Loads;

DomainWrapper::DomainWrapper()
{
	_Domain = new Domain();
}

DomainWrapper::~DomainWrapper()
{
	if (_Domain != 0)
	{
		_Domain->~Domain();
		delete _Domain;
	}
}

bool ^ 
DomainWrapper::AddNode(array<NodeWrapper^>^ nodes)
{
	int length = (nodes->Length);
	for (int i = 0; i < length; i++) {
		Node* _node = nodes[i]->_Node;

		if(! _Domain->addNode(_node))
			return false;
		AddNodeEventHandler(nullptr,nullptr);
	}
	
	return true;
}

bool^
DomainWrapper::AddNode(NodeWrapper^ node)
{
	Node* _node = node->_Node;
	return _Domain->addNode(_node);
	//return false;
}

bool ^ DomainWrapper::AddElement(array<ElementWrapper^>^ elements)
{
	int length = (elements->Length);
	for (int i = 0; i < length; i++) {
		Element* _element = elements[i]->_Element;
		if (!_Domain->addElement(_element))
			return false;
	}
	return true;
}

bool^ 
DomainWrapper::AddElement(ElementWrapper^ element)
{
	Element* _element = element->_Element;
	return _Domain->addElement(_element);
}

bool ^ DomainWrapper::AddSP_Constraint(array<SP_ConstraintWrapper^>^ sps)
{
	int length = (sps->Length);
	for (int i = 0; i < length; i++) {
		SP_Constraint* _sP_Constraint = sps[i]->_SP_Constraint;
		if (!_Domain->addSP_Constraint(_sP_Constraint))
			return false;
	}
	return true;
}

bool^
DomainWrapper::AddSP_Constraint(SP_ConstraintWrapper^ sP_Constraint)
{
	SP_Constraint* _sP_Constraint = sP_Constraint->_SP_Constraint;
	return _Domain->addSP_Constraint(_sP_Constraint);
}

bool^
DomainWrapper::AddMP_Constraint(MP_ConstraintWrapper^ mp_Constraint)
{
	return _Domain->addMP_Constraint(mp_Constraint->_MP_Constraint);
}

bool^ 
DomainWrapper::AddMP_Constraint(array<MP_ConstraintWrapper^>^ mps)
{
	int length = (mps->Length);
	for (int i = 0; i < length; i++) {
		if (!_Domain->addMP_Constraint(mps[i]->_MP_Constraint))
			return false;
	}
	return true;
}

bool^
DomainWrapper::AddLoadPattern(LoadPatternWrapper^ loadPattern)
{
	LoadPattern* _loadPattern = loadPattern->_LoadPattern;
	return _Domain->addLoadPattern(_loadPattern);
}

bool^
DomainWrapper::AddNodalLoad(NodalLoadWrapper^ nodalLoad, int loadPatternTag)
{
	NodalLoad* _nodalLoad = nodalLoad->_NodalLoad;
	return _Domain->addNodalLoad(_nodalLoad, loadPatternTag);
}

bool^ 
DomainWrapper::AddNodalLoad(array<NodalLoadWrapper^>^ nodalLoads, int loadPatternTag)
{
	int length = (nodalLoads->Length);
	for (int i = 0; i < length; i++) {
		NodalLoad* nl = nodalLoads[i]->_NodalLoad;
		if (!_Domain->addNodalLoad(nl,loadPatternTag))
			return false;
	}
	return true;
}
