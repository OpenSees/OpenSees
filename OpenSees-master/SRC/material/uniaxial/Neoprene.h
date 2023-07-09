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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-12-09 21:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Neoprene.h,v $

// Written: krm
// Created: 12/2005
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterial. Neoprene provides the abstraction
// of an elastic perfectly plastic (tension only) path dependent uniaxial 
// material, with an initial gap offset (force-displacement units)
// For compression only behavior, enter negative gap and ep
// Damage can accumulate through specification of damage = 1 switch,
// otherwise damage = 0
//
//                  gap ep
//                  |<>|>|
//          stress  |    +-----+ fy
//                  |   /E    /
//         ---------+--+-+---+----strain
//                  |       
//                  |   
//                  |      
//
// What: "@(#) Neoprene.h, revA"

#ifndef Neoprene_h
#define Neoprene_h

#include <UniaxialMaterial.h>

class Neoprene : public UniaxialMaterial
{
  public:
    Neoprene(int tag, double E, double gap);
    Neoprene();  
    ~Neoprene();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double commitStrain;
    double trialStrain;
    double E;
    double gap;
    double minElasticYieldStrain;
    double maxElasticYieldStrain;
    
    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
};


#endif
