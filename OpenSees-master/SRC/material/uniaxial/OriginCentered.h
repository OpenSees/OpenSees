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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/OriginCentered.h,v $

//
// Description: This file contains the class definition for 
// OriginCentered. 

#ifndef OriginCentered_h
#define OriginCentered_h

#include <UniaxialMaterial.h>

class OriginCentered : public UniaxialMaterial
{
  public:
    OriginCentered(int tag,
	    double f1, double e1, double f2,
	    double e2, double f3, double e3);
	    
    OriginCentered(void);
    virtual ~OriginCentered();

    const char *getClassType(void) const {return "OriginCentered";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);      
    double getStress(void);
    double getTangent(void);
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
 protected:
    
 private:
    // matpar : STEEL FIXED PROPERTIES
    double f1;  //  = matpar(1)
    double e1;  //  = matpar(2) 
    double f2;   //  = matpar(3) 
    double e2;  //  = matpar(4) 
    double f3; //  = matpar(5)
    double e3; //  = matpar(6) 

	double E1;
	double E2;
	double E3;

    // 
    double TepsMax;
    double TepsMin;
    double TsigMax;
    double TsigMin;
    
	double Tsig;
	double Teps;
	double Ttan;

    double CepsMax;
    double CepsMin;
    double CsigMax;
    double CsigMin;
    
	double Csig;
	double Ceps;
	double Ctan;
};


#endif

