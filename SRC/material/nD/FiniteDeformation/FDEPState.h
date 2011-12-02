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
#ifndef FDEPState_H
#define FDEPState_H

#include <Tensor.h>
#include <stresst.h>
#include <straint.h>

class FDEPState  {    
    public:
      FDEPState();
      FDEPState ( const straintensor& xFpInVar,
	          double xStressLikeInVar,
	          double xStrainLikeInVar,
		  const stresstensor& xStressLikeKiVar,
		  const straintensor& xStrainLikeKiVar,
		  const straintensor& xCommitedFpInVar,
	          double xCommitedStressLikeInVar,
	          double xCommitedStrainLikeInVar,
	          const stresstensor& xCommitedStressLikeKiVar,
	          const straintensor& xCommitedStrainLikeKiVar);	
      ~FDEPState();
		
	FDEPState *newObj();
	FDEPState( const FDEPState& fdeps );
	
	straintensor getFpInVar() const;
	double getStressLikeInVar() const;
	double getStrainLikeInVar() const;
	stresstensor getStressLikeKiVar() const;
	straintensor getStrainLikeKiVar() const;
	
	straintensor getCommitedFpInVar() const;
	double getCommitedStressLikeInVar() const;
	double getCommitedStrainLikeInVar() const;
	stresstensor getCommitedStressLikeKiVar() const;
	straintensor getCommitedStrainLikeKiVar() const;
	
	void setFpInVar(const straintensor& xFpInVar);
	void setStressLikeInVar(double xStressLikeInVar);
	void setStrainLikeInVar(double xStrainLikeInVar);
	void setStressLikeKiVar(const stresstensor& xStressLikeKiVar);
	void setStrainLikeKiVar(const straintensor& xStrainLikeKiVar);
	
	void setCommitedFpInVar(const straintensor& xCommitedFpInVar);
	void setCommitedStressLikeInVar(double xCommitedStressLikeInVar);
	void setCommitedStrainLikeInVar(double xCommitedStrainLikeInVar);
	void setCommitedStressLikeKiVar(const stresstensor& xCommitedStressLikeKiVar);
	void setCommitedStrainLikeKiVar(const straintensor& xCommitedStrainLikeKiVar);

        int commitState(void) ;
        int revertToLastCommit(void) ;
        int revertToStart(void) ;

    private:
        straintensor FpInVar;                  // Plastic F;
        double StressLikeInVar;                // Scalar stress like evolution variable
        double StrainLikeInVar;                // Scalar strain like evolution variable
        stresstensor StressLikeKiVar;          // Tensor stress like evolution variable
        straintensor StrainLikeKiVar;          // Tensor strain like evolution variable


        straintensor CommitedFpInVar;          // Plastic F;
        double CommitedStressLikeInVar;        // Scalar stress like evolution variable
        double CommitedStrainLikeInVar;        // Scalar strain like evolution variable
        stresstensor CommitedStressLikeKiVar;  // Tensor stress like evolution variable
        straintensor CommitedStrainLikeKiVar;  // Tensor strain like evolution variable

};


#endif
