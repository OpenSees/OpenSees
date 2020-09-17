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
                                                                       
// Written: Alborz Ghofrani, Pedro Arduino
//			May 2013, University of Washington
                                                                      
// Description: This file contains the implementation for the ManzariDafaliasRO class.

#ifndef ManzariDafaliasRO_h
#define ManzariDafaliasRO_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

#include <Information.h>
//#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

#include "ManzariDafalias.h"

#include <elementAPI.h>

class ManzariDafaliasRO : public ManzariDafalias
{
  public:

    // full constructor
    ManzariDafaliasRO(int tag, int classTag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c, 
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd, 
					double z_max, double cz, double mDen, double kappa = 2.0, int integrationScheme = 2, int tangentType = 2, int JacoType = 1, 
					double TolF = 1.0e-7, double TolR = 1.0e-7);
    // null constructor
    ManzariDafaliasRO();
    // destructor
    ~ManzariDafaliasRO();
 
    NDMaterial *getCopy(const char *type);

	void integrate();
    int commitState(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
	int			getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    void Print(OPS_Stream &s, int flag =0);
	
  protected:

	double m_B;
	double m_a1;
	double m_gamma1;
	double m_kappa;
	
	// internal variables
	Vector mDevEpsSR;	// deviatoric strain at shear reversal point
	Vector mSigmaSR;	// stress at shear reversal point
	double mChi_r;		// Chi_r Ramberg-Osgood parameter
	double mDChi_e;		// change in Chi_e Ramberg-Osgood parameter
	double mEta1;		// eta1 Ramber-Osgood parameter
	bool   mIsFirstShear; // boolean to determine if it's first shearing

	//Member Functions specific for ManzariDafaliasRO model
	void	initialize();
	void	GetElasticModuli(const Vector& sigma, const double& en, const double& en1,
				const Vector& nEStrain, const Vector& cEStrain, double &K, 
				double &G);
	void	GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G);
};

#endif
