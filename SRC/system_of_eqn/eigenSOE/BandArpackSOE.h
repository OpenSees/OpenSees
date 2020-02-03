// File: ~/system_of_eqn/eigenSOE/BandArpackSOE.h
//
// Written: Jun Peng 
// Created: February 1999
// Revision: A
//
// Description: This file contains the class definition for BandArpackSOE
// BandArpackSOE is a subclass of EigenSOE. It uses the LAPACK storage
// scheme to store the components of the K, M matrix, which is a full matrix.
// It uses the ARPACK to do eigenvalue analysis.
// ARPACK is an eigen analysis package which was developed by 
// R.B.Lehoucq, D.C.Sorensen and C.Yang at Rice University. ARPACK is a
// collection of FORTRAN77 subroutines designed to solve large scale eigen
// problems. ARPACK is capable of solving large scale non-Hermitian standard 
// and generalized eigen problems. When the matrix <B>K</B> is symmetric, 
// the method is a variant of the Lanczos process called Implicitly Restarted
// Lanczos Method (IRLM).


#ifndef BandArpackSOE_h
#define BandArpackSOE_h

#include <EigenSOE.h>
#include <Vector.h>

class AnalysisModel;
class BandArpackSolver;

class BandArpackSOE : public EigenSOE
{
  public:
    BandArpackSOE(BandArpackSolver &theSolver, AnalysisModel &theModel, 
		  double shift = 0.0);    

    virtual ~BandArpackSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    
   
    virtual void zeroA(void);
    virtual void zeroM(void);
    
    virtual double getShift(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
    friend class BandArpackSolver;

  protected:
    
  private:
    int size, numSuperD, numSubD;    
    double *A;
    int Asize;
    bool factored;
    double shift;
    AnalysisModel *theModel;
};


#endif



