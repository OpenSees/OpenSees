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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-08-31 11:36:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMLinSOE.h,v $

// Written: Minjie Zhu
// Created: August 2012
//
// Description: This file contains the class definition for PFEMLinSOE
// PFEMLinSOE is a subclass of SparseGenColLinSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. It solves the equations using the Fractional Step Method in PFEM.
//
// What: "@(#) PFEMLinSOE.h, revA"

#ifndef PFEMLinSOE_h
#define PFEMLinSOE_h


#include <SparseGenColLinSOE.h>
#include <ID.h>
extern "C" {
#include <cs.h>
}

class PFEMLinSOE : public SparseGenColLinSOE
{
  public:
    PFEMLinSOE(SparseGenColLinSolver &theSolver);        
    PFEMLinSOE(int classTag);        
    PFEMLinSOE();        
    PFEMLinSOE(SparseGenColLinSolver &theSolver, int classTag);        


    virtual ~PFEMLinSOE();

    virtual int solve();
    virtual int setLinks(AnalysisModel& model);
    virtual int setSize(Graph& theGraph);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    virtual int setDofIDs();
    virtual cs* setMatIDs();

  protected:
    AnalysisModel* getAnalysisModel();
    
  private:
    AnalysisModel* theModel;
    ID* mID, *pID, *piID, *mIDall;
    ID* Mid, *Mhatid, *Gid, *Gtid, *Lid, *Qtid;
    cs* G, *Gt, *L, *Qt;
};


#endif

