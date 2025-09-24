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

#ifndef SymmGeneralizedEigenSOE_h
#define SymmGeneralizedEigenSOE_h

#include <EigenSOE.h>
#include <Vector.h>

class AnalysisModel;
class SymmGeneralizedEigenSolver;

class SymmGeneralizedEigenSOE : public EigenSOE
{
public:
    SymmGeneralizedEigenSOE(SymmGeneralizedEigenSolver &theSolver,
        AnalysisModel &theModel);

    virtual ~SymmGeneralizedEigenSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    

    virtual void zeroA(void);
    virtual void zeroM(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    friend class SymmGeneralizedEigenSolver;

protected:

private:
    int size;
    double *A;
    int Asize;
    double *M;
    int Msize;
    bool factored;
    AnalysisModel *theModel;
};

#endif
