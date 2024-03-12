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
      

#ifndef BoucWenInfill_h
#define BoucWenInfill_h

// Written by Stefano Sirotti (stefano.sirotti@unimore.it)
// Created on January 2022
//
// Description: This file contains the class definition for 
// BoucWenInfill material. BoucWenInfill material provides 
// a hysteretic uniaxial material described by the Bouc-Wen law 
// with stiffness, strength degradation and pinching. 
// Particularly suitable to simulate hysteresis of infill panels.
//
// Reference: Sirotti, S., Pelliciari, M., Di Trapani, F.,
// Briseghella, B., Carlo Marano, G., Nuti, C., & Tarantino, A. M. (2021).
// Development and validation of new bouc�wen data-driven hysteresis model 
// for masonry infilled rc frames. 
// Journal of Engineering Mechanics, 147(11), 04021092.
// 

#include <UniaxialMaterial.h>
#include <Matrix.h>



class BoucWenInfill : public UniaxialMaterial
{
  public:
    BoucWenInfill(int tag,
     double mass,
	 double alpha,
	 double beta0,
	 double eta0,
	 double n,
	 double k,
	 double xy,
	 double deltak,
	 double deltaf,
	 double psi,
	 double Zs,
	 double As,
	 double epsp,
	 double tolerance,
	 int maxNumIter);
	
    BoucWenInfill();
    ~BoucWenInfill();  

    const char *getClassType(void) const {return "BoucWenMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double signum(double);
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);
    
    double getInitialTangent(void);

  protected:
    
  private:

    // Material parameters
    double mass;
	double alpha;
	double beta0;
	double eta0;
	double n;
	double k;
	double xy;
	double deltak;
	double deltaf;
	double psi;
	double Zs;
	double As;
	double epsp;
    
    // History variables (trial and committed)
	double xmaxp;
    double xmax;
	double Tstrain, Cstrain;
    double Tz, Cz;
    double Te, Ce;
    
    // Other variables
    double Tstress, Ttangent;
    
    double tolerance;
    int maxNumIter;
};


#endif



