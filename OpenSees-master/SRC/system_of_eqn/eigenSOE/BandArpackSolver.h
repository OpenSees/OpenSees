// Written: Jun Peng
// Created: Feb. 11, 1999
// Revision: A
//
// Description: This file contains the class definition for 
// BandArpackSolver. It solves the BandArpackSOE object by calling
// Arpack routines.


#ifndef BandArpackSolver_h
#define BandArpackSolver_h

#include <EigenSolver.h>
#include <BandArpackSOE.h>

class BandArpackSolver : public EigenSolver
{
  public:
    BandArpackSolver(int numE = 0);    
    virtual ~BandArpackSolver();

    virtual int solve(int numModes, bool generalized, bool findSmallest = true);  

    virtual int setSize(void);
    virtual int setEigenSOE(BandArpackSOE &theSOE);
    
    virtual const Vector &getEigenvector(int mode);
    virtual double getEigenvalue(int mode);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    

  private:
    BandArpackSOE *theSOE;
    int theNev;    
    double *value;
    double *eigenvector;
    int *iPiv;
    int iPivSize;
    Vector *eigenV;
    
    void myMv(int n, double *v, double *result);
    void myCopy(int n, double *v, double *result);
    int getNCV(int n, int nev);
    
};

#endif


