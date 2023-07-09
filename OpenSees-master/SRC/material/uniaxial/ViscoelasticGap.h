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
//                                                                     
// Revision: 1.0
// Date: 05/2019
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscoelasticGap.h
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 05/2019
//
// Description: This file contains the class definition for the
// ViscoelasticGap uniaxialMaterial.
//
// References:
// Goldsmith, W. (1960). "Impact: The Theory and Physical Behavior of Colliding Solids." 
//   E. Arnold: London.
//
// Variables:
// K: linear stiffness
// C: linear damping coefficient
// gap: initial gap (must be a negative value)


#ifndef ViscoelasticGap_h
#define ViscoelasticGap_h

#include <UniaxialMaterial.h>

class ViscoelasticGap : public UniaxialMaterial
{
  public: 
    ViscoelasticGap(int tag, double K, double C, double gap);   
    ViscoelasticGap(); 
    ~ViscoelasticGap();

    const char *getClassType(void) const {return "ViscoelasticGap";};

    int setTrialStrain(double strain, double strainRate); 
	
    double getStrain(void); 
    double getStrainRate(void);
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
    // input variables
    double K; // linear stiffness
    double C; // linear damping coefficient
    double gap; // initial gap distance (must be a negative value)
    
    // state variables
	double commitStrain;
	double commitStrainRate;
	double commitStress;
	double commitTangent;
	double trialStrain;
	double trialStrainRate;
	double trialStress;
	double trialTangent;

	// other variables
	int printFlag; // print flag for impact event
	
};

#endif