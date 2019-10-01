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

#define sdm2d_ SDM2D

extern "C"  {void sdm2d_(
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
			 double (*_tangent)[3]);}

#else

extern "C"  {void sdm2d_(
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
			 double (*_tangent)[3]);}
#endif

class StressDensityModel2D : public StressDensityModel {

public:
	
	//default constructor 
	StressDensityModel2D();

	//initialisation constructor
	StressDensityModel2D(int tag, double constDensity,
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
                           double constRxx, double constRyy, double constRxy);

	//destructor
	~StressDensityModel2D();

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

    double hard_para_real[250];

    // parameters sent to FORTRAN 
	double _stress_current[3], 
		   _strain_current[3], 
		   _strain_next[3], 
		   _model_parameter[16], 
		   _ssl_void_ratio[10], 
		   _ssl_pressure[10], 
		   _hsl_void_ratio[10], 
		   _hsl_pressure[10], 
           _hard_para_real[250],
		   _tangent[3][3];
    int _hard_para_int[2];

    // member functions
    void initialise();
    void CalInitialTangent(void);
	void GetCurrentStress(void);
};
