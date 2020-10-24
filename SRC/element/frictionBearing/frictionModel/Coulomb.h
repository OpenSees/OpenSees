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

#ifndef Coulomb_h
#define Coulomb_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for the Coulomb
// friction model. In the Coulomb model the friction force is given by
// mu*N, where mu is a constant coefficient of friction and N is a positive
// normal force perpendicular to the sliding surface. If N is negative
// the friction force is zero.

#include "FrictionModel.h"

class Coulomb : public FrictionModel
{
public:
    // constructor
    Coulomb();
    Coulomb(int tag, double mu);
    
    // destructor
    ~Coulomb();
    
    const char *getClassType() const {return "Coulomb";};
    
    // public methods to set and obtain response
    int setTrial(double normalForce, double velocity = 0.0);
    double getFrictionForce();
    double getFrictionCoeff();
    double getDFFrcDNFrc();
    double getDFFrcDVel();
    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    
    FrictionModel *getCopy();
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:

private:
    double mu;  // coefficient of friction (COF)
};

#endif
