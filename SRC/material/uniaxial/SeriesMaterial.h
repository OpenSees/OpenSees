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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SeriesMaterial.h,v $

#ifndef SeriesMaterial_h
#define SeriesMaterial_h

// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the class definition for 
// SeriesMaterial. SeriesMaterial is an aggregation
// of UniaxialMaterial objects all considered acting in Series.

#include <UniaxialMaterial.h>

class SeriesMaterial : public UniaxialMaterial
{
  public:
    SeriesMaterial(int tag, 
			  int numMaterial, 
			  UniaxialMaterial **theMaterials,
			  int maxIter = 1, double tol = 1.0e-10);
    SeriesMaterial();
    ~SeriesMaterial();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

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
    double Tstrain;
	double Cstrain;
	double Tstress;
	double Cstress;
	double Ttangent;
	double Ctangent;

	int maxIterations;
	double tolerance;

	double *stress;
	double *flex;
	double *strain;

	bool initialFlag;

    int numMaterials;   // the number of UniaxialMaterials in the aggregation
    UniaxialMaterial **theModels; // an array of pointers to the UniaxialMaterials
};


#endif

