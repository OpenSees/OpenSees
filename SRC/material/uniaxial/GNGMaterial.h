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
                                                                        
// 'Grip n Grab' ratcheting, tension-only device material model
// - file edited from EPPGapMaterial file
// Jarrod Cook, University of Canterbury, Christchurch, New Zealand

//	^
//  |
//	|                ________(3)________
//  |               /                  /
//	F              /                  /
//	O             /                  /
//	R            /                  /
//	C          (2)                (4)
//	E          /                  /
//	|         /                  /
//	|        /                  /
//	|__(1)__/     _____(5)_____/
//  --------------DISPLACEMENT---------->
//
//	LOADING
// 		BELOW ENGAGEMENT THRESHOLD (1)
// 		ELASTIC REGION (2)
//		BEYOND YIELD (3)
//	UNLOADING
//		ELASTIC RECOVERY (4) 
// 		BELOW ENGAGEMENT THRESHOLD (5)

#ifndef GNGMaterial_h
#define GNGMaterial_h

#include <UniaxialMaterial.h>
#include <Information.h>

class GNGMaterial : public UniaxialMaterial
{
  public:
  
	//full constructor
    GNGMaterial(int tag, double E, double sigY, double P, double eta);
	
	//null constructor
    GNGMaterial();  
	
	//destructor
    ~GNGMaterial();

    const char *getClassType(void) const {return "GNGMaterial";};

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
    
	Response* setResponse (const char **argv, int argc, OPS_Stream &theOutput);
	int getResponse (int responseID, Information &matInfo);
	
  protected:
    
  private:
  
    double commitStrain;
    double trialStrain;
    double E;
    double sigY;
    double P;
    double eta;
    double epsY;
    double epsE;
	double epsP;
	double sigP;
	double pdemand; //cumulative plastic demand
	int nratchet; //ratchet count

    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
	
};


#endif
