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

// Written: Saumyasuchi Das
//          May 2013, University of Canterbury
// Updated: Chris McGann
//          June 2015, Washington State University

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include <StressDensityModel.h>

#ifdef _WIN32

#define sdm3d_ SDM3D

extern "C"  {void sdm3d_(
			 double _stress_current[], 
			 double _strain_current[], 
			 double _strain_next[], 
			 double _model_parameter[], 
			 double _ssl_void_ratio[], 
			 double _ssl_pressure[], 
			 double _hsl_void_ratio[], 
			 double _hsl_pressure[], 
			 double _hard_para_real[], 
			 int    _hard_para_int[],
			 double _anisotropy[],
			 double (*_tangent)[6]);}

#else

extern "C"  {void sdm3d_(
			 double _stress_current[], 
			 double _strain_current[], 
			 double _strain_next[], 
			 double _model_parameter[], 
			 double _ssl_void_ratio[], 
			 double _ssl_pressure[], 
			 double _hsl_void_ratio[], 
			 double _hsl_pressure[], 
			 double _hard_para_real[], 
			 int    _hard_para_int[],
			 double _anisotropy[],
			 double (*_tangent)[6]);}
#endif

class StressDensityModel3D : public StressDensityModel {

public:
	
	//default constructor 
	StressDensityModel3D();

	//initialisation constructor
	StressDensityModel3D(int tag, double constDensity,
                           // SD model  parameters
                           double initialVoidRatio,	double constA, double exponentN, double poissonRatio,	
                           double constAlpha1, double constBeta1, double constAlpha2, double constBeta2,
                           double constAlpha3, double constBeta3, double constDegradation, double constMumin,
                           double constMucyclic, double constDilatancyStrain, double constMumax, double constPatm,
                           // steady state line void ratio
                           double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
                           double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
                           double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
                           double constsslvoidatP10,
                           // hydrostatic state line void ratio
                           double consthslvoid, 
                           // reference pressures 
                           double constP1, double constP2, double constP3, double constP4, double constP5,
                           double constP6, double constP7, double constP8, double constP9, double constP10,
                           // offset of the failure surface
                           double constRxx, double constRyy, double constRzz, double constRxy,
                           double constRyz, double constRzx);

	//destructor
	~StressDensityModel3D();

    NDMaterial* getCopy( ) ;
    const char* getType( ) const ;
    int getOrder( ) const ;

    int setTrialStrain(const Vector &strain_from_element);

    // Unused trialStrain functions
    int setTrialStrain(const Vector &v, const Vector &r);
    
    //send back the strain
    const Vector& getStrain( ) ;

    //send back the stress 
    const Vector& getStress( ) ;

    int commitState(void);

    //send back the tangent 
    const Matrix& getTangent( ) ;
    const Matrix& getInitialTangent( ) ;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  protected:

    Vector stressCurrent,
		   stressNext,
		   strainCurrent,
		   strainNext;

	Matrix initialTangent,
		   currentTangent;

    double hard_para_real[460];

    // parameters sent to FORTRAN 
	double _stress_current[6], 
		   _strain_current[6], 
		   _strain_next[6], 
		   _model_parameter[16], 
		   _ssl_void_ratio[10], 
		   _ssl_pressure[10], 
		   _hsl_void_ratio[10], 
		   _hsl_pressure[10], 
           _hard_para_real[460],
		   _anisotropy[6],
		   _tangent[6][6];
    int _hard_para_int[2];

    // member functions
    void initialise();
    void CalInitialTangent(void);
	void GetCurrentStress(void);
};
