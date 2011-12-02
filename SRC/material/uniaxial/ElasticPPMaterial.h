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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticPPMaterial.h,v $
                                                                        
                                                                        
#ifndef ElasticPPMaterial_h
#define ElasticPPMaterial_h

// File: ~/material/ElasticPPMaterial.h
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterial. ElasticPPMaterial provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
//                  ep
//                  |>|
//          stress  | +-----+
//                  |/E    /
//             -----+-+---+----strain
//                  |    /
//               ---+---+
//                      |<->|
//                       2ep
//
// What: "@(#) ElasticPPMaterial.h, revA"

#include <UniaxialMaterial.h>

class ElasticPPMaterial : public UniaxialMaterial
{
  public:
    ElasticPPMaterial(int tag, double E, double ep);    
	ElasticPPMaterial();    
    ~ElasticPPMaterial();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getSecant (void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(ostream &s, int flag =0);
    
  protected:
    
  private:
    double commitStrain;
    double trialStrain;
    double E;
    double ep;    
    double maxElasticYieldStrain;
    double minElasticYieldStrain;
};


#endif

