#include "stdafx.h"
#include "ConstraintWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components::Constraints;

SP_ConstraintWrapper::SP_ConstraintWrapper()
{

}

SP_ConstraintWrapper::SP_ConstraintWrapper(int nodeTag, int ndof, double value, bool isConstant)
{
	_SP_Constraint = new SP_Constraint(nodeTag, ndof, value, isConstant);
}


ImposedMotionSPWrapper::ImposedMotionSPWrapper(int nodeTag, int ndof, int patternTag, int theGroundMotionTag)
{
	_SP_Constraint = new ImposedMotionSP(nodeTag, ndof, patternTag, theGroundMotionTag);
}

ImposedMotionSP1Wrapper::ImposedMotionSP1Wrapper(int nodeTag, int ndof, int patternTag, int theGroundMotionTag)
{
	_SP_Constraint = new ImposedMotionSP1(nodeTag, ndof, patternTag, theGroundMotionTag);
}

MP_ConstraintWrapper::MP_ConstraintWrapper(int nodeRetain,
	int nodeConstr,
	IDWrapper^ constrainedDOF,
	IDWrapper^ retainedDOF,
	int classTag)
{
	_MP_Constraint = new MP_Constraint(nodeRetain, nodeConstr, *constrainedDOF->_ID, *retainedDOF->_ID, classTag);
}

MP_ConstraintWrapper::MP_ConstraintWrapper(int nodeRetain,
	int nodeConstr,
	MatrixWrapper^ constrnt,
	IDWrapper^ constrainedDOF,
	IDWrapper^ retainedDOF)
{
	_MP_Constraint = new MP_Constraint(nodeRetain, nodeConstr, *constrnt->_Matrix, *constrainedDOF->_ID, *retainedDOF->_ID);
}



