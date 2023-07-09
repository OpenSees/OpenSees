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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: Minjie Zhu

// Description: Quasi-incompressible solver, M and Mp


#ifndef PFEMQuasiLinSOE_h
#define PFEMQuasiLinSOE_h


#include <PFEMLinSOE.h>

class PFEMQuasiSolver;
// class PFEMQuasiSolver_Mumps;

class PFEMQuasiLinSOE : public PFEMLinSOE
{
  public:
    PFEMQuasiLinSOE(PFEMQuasiSolver &theSolver);        
    PFEMQuasiLinSOE();        

    virtual ~PFEMQuasiLinSOE();

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual void zeroA(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    virtual const ID& getDofID()const {return newDofID;}

    friend class PFEMQuasiSolver;
    // friend class PFEMQuasiSolver_Mumps;

private:    

    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize);

private:

    cs* M, *Mp, *Gt;
    ID newDofID;
};

#endif

