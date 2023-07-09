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
                                                                        
// $Revision: 1.3 $
// $Date: 2009-05-11 20:56:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.h,v $
                                                                        
// Written: fmk 
// Created: 11/98
//
// Description: This file contains the class definition for 
// UmfpackGenLinSOE. It stores the sparse matrix A in a fashion
// required by the UmfpackLinSolver object.
//
// What: "@(#) UmfpackGenLinSOE.h, revA"


#ifndef UmfpackGenLinSOE_h
#define UmfpackGenLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>
#include <vector>

class UmfpackGenLinSolver;

class UmfpackGenLinSOE : public LinearSOE
{
public:
    UmfpackGenLinSOE(UmfpackGenLinSolver &theSolver);
    UmfpackGenLinSOE();        

    ~UmfpackGenLinSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setUmfpackGenLinSolver(UmfpackGenLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker); 

    friend class UmfpackGenLinSolver;

protected:
    
private:
    Vector X,B;
    std::vector<int> Ap, Ai;
    std::vector<double> Ax;
};


#endif

