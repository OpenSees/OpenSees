#include "stdafx.h"
#include "Node.h"

#include "NodeWrapper.h"
using namespace OpenSees;
using namespace OpenSees::Components;

NodeWrapper::NodeWrapper()
{
	
}

NodeWrapper::NodeWrapper(int tag, int ndf, double crd1)
{
	this->_Node = new Node(tag, ndf, crd1);
	_TaggedObject = this->_Node;
}

NodeWrapper::NodeWrapper(int tag, int ndf, double crd1, double crd2)
{
	this->_Node = new Node(tag, ndf, crd1,crd2);
	_TaggedObject = this->_Node;
}

NodeWrapper::NodeWrapper(int tag, int ndf, double crd1, double crd2, double crd3)
{
	this->_Node = new Node(tag, ndf, crd1, crd2, crd3);
	_TaggedObject = this->_Node;
}

NodeWrapper::~NodeWrapper()
{
	if (_Node != 0)
		delete _Node;
}




