//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================
#ifndef FDEPState_CPP
#define FDEPState_CPP

#include "FDEPState.h"

//-------------------------------------------------------------------------------
// Normal Constructor 00
//-------------------------------------------------------------------------------
FDEPState::FDEPState ()
{
	int err;
	err = this->revertToStart();
}

//-------------------------------------------------------------------------------
// Normal Constructor 01
//-------------------------------------------------------------------------------
FDEPState::FDEPState ( const straintensor& xFpInVar,
	               double xStressLikeInVar,
	               double xStrainLikeInVar,
		       const stresstensor& xStressLikeKiVar,
		       const straintensor& xStrainLikeKiVar,
		       //
		       const straintensor& xCommitedFpInVar,
	               double xCommitedStressLikeInVar,
	               double xCommitedStrainLikeInVar,
		       const stresstensor& xCommitedStressLikeKiVar,
		       const straintensor& xCommitedStrainLikeKiVar )
{
	FpInVar = xFpInVar;
	StressLikeInVar = xStressLikeInVar;
	StrainLikeInVar = xStrainLikeInVar;
	StressLikeKiVar = xStressLikeKiVar;
	StrainLikeKiVar = xStrainLikeKiVar;
	//
	CommitedFpInVar = xCommitedFpInVar;
	CommitedStressLikeInVar = xCommitedStressLikeInVar;
	CommitedStrainLikeInVar = xCommitedStrainLikeInVar;
	CommitedStressLikeKiVar = xCommitedStressLikeKiVar;
	CommitedStrainLikeKiVar = xCommitedStrainLikeKiVar;
}

//-------------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------------
FDEPState::~FDEPState ()
{

}

//------------------------------------------------------------------------------
FDEPState *FDEPState::newObj ()
{
    FDEPState *fdeps = new FDEPState ();
    return fdeps;
}

//------------------------------------------------------------------------------
FDEPState::FDEPState( const FDEPState& fds )
{
	setFpInVar(fds.getFpInVar());
	setStressLikeInVar(fds.getStressLikeInVar());
	setStrainLikeInVar(fds.getStrainLikeInVar());
	setStressLikeKiVar(fds.getStressLikeKiVar());
	setStrainLikeKiVar(fds.getStrainLikeKiVar());
	//
	setCommitedFpInVar(fds.getCommitedFpInVar());
	setCommitedStressLikeInVar(fds.getCommitedStressLikeInVar());
	setCommitedStrainLikeInVar(fds.getCommitedStrainLikeInVar());
	setCommitedStressLikeKiVar(fds.getCommitedStressLikeKiVar());
	setCommitedStrainLikeKiVar(fds.getCommitedStrainLikeKiVar());
}

// Get member variables
//------------------------------------------------------------------------------
straintensor FDEPState::getFpInVar() const
{
	return FpInVar;
}

//------------------------------------------------------------------------------
double FDEPState::getStrainLikeInVar() const
{
	return StrainLikeInVar;
}

//------------------------------------------------------------------------------
double FDEPState::getStressLikeInVar() const
{
	return StressLikeInVar;
}

//------------------------------------------------------------------------------
straintensor FDEPState::getStrainLikeKiVar() const
{
	return StrainLikeKiVar;
}

//------------------------------------------------------------------------------
stresstensor FDEPState::getStressLikeKiVar() const
{
	return StressLikeKiVar;
}

//------------------------------------------------------------------------------
straintensor FDEPState::getCommitedFpInVar() const
{
	return CommitedFpInVar;
}

//------------------------------------------------------------------------------
double FDEPState::getCommitedStrainLikeInVar() const
{
	return CommitedStrainLikeInVar;
}

//------------------------------------------------------------------------------
double FDEPState::getCommitedStressLikeInVar() const
{
	return CommitedStressLikeInVar;
}

//------------------------------------------------------------------------------
straintensor FDEPState::getCommitedStrainLikeKiVar() const
{
	return CommitedStrainLikeKiVar;
}

//------------------------------------------------------------------------------
stresstensor FDEPState::getCommitedStressLikeKiVar() const
{
	return CommitedStressLikeKiVar;
}

// Set member variables
//---------------------------------------------------------------------------
void FDEPState::setFpInVar(const straintensor &xFpInVar)
{
	this->FpInVar = xFpInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setStrainLikeInVar(double xStrainLikeInVar)
{
	this->StrainLikeInVar = xStrainLikeInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setStressLikeInVar(double xStressLikeInVar)
{
	this->StressLikeInVar = xStressLikeInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setStrainLikeKiVar(const straintensor& xStrainLikeKiVar)
{
	this->StrainLikeKiVar = xStrainLikeKiVar;
}

//---------------------------------------------------------------------------
void FDEPState::setStressLikeKiVar(const stresstensor& xStressLikeKiVar)
{
	this->StressLikeKiVar = xStressLikeKiVar;
}

//---------------------------------------------------------------------------
void FDEPState::setCommitedFpInVar(const straintensor &xCommitedFpInVar)
{
	this->CommitedFpInVar = xCommitedFpInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setCommitedStrainLikeInVar(double xCommitedStrainLikeInVar)
{
	this->CommitedStrainLikeInVar = xCommitedStrainLikeInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setCommitedStressLikeInVar(double xCommitedStressLikeInVar)
{
	this->CommitedStressLikeInVar = xCommitedStressLikeInVar;
}

//---------------------------------------------------------------------------
void FDEPState::setCommitedStrainLikeKiVar(const straintensor& xCommitedStrainLikeKiVar)
{
	this->CommitedStrainLikeKiVar = xCommitedStrainLikeKiVar;
}

//---------------------------------------------------------------------------
void FDEPState::setCommitedStressLikeKiVar(const stresstensor& xCommitedStressLikeKiVar)
{
	this->CommitedStressLikeKiVar = xCommitedStressLikeKiVar;
}

//----------------------------------------------------------------------
int FDEPState::commitState(void)
{
        CommitedFpInVar = FpInVar;
        CommitedStressLikeInVar = StressLikeInVar;
        CommitedStrainLikeInVar = StrainLikeInVar;
        CommitedStressLikeKiVar = StressLikeKiVar;
        CommitedStrainLikeKiVar = StrainLikeKiVar;

	return 0;
}

//----------------------------------------------------------------------
int FDEPState::revertToLastCommit(void)
{
        FpInVar = CommitedFpInVar;
        StressLikeInVar = CommitedStressLikeInVar;
        StrainLikeInVar = CommitedStrainLikeInVar;
        StressLikeKiVar = CommitedStressLikeKiVar;
        StrainLikeKiVar = CommitedStrainLikeKiVar;

	return 0;
}

//----------------------------------------------------------------------
int FDEPState::revertToStart(void)
{
	tensor tI2("I", 2, def_dim_2);
	tensor t00(2, def_dim_2, 0.0);

	FpInVar = tI2;
	StressLikeInVar = 0.0;
	StrainLikeInVar = 0.0;
	StressLikeKiVar = t00;
	StrainLikeKiVar = t00;

	CommitedFpInVar = tI2;
	CommitedStressLikeInVar = 0.0;
	CommitedStrainLikeInVar = 0.0;
	CommitedStressLikeKiVar = t00;
	CommitedStrainLikeKiVar = t00;

	return 0;
}

# endif







