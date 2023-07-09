
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Fluid Material
//------------------------------------------


#ifndef AcousticMedium_h
#define AcousticMedium_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#define ND_TAG_AcousticMedium  200001

class AcousticMedium : public NDMaterial
{
  public:
    // Only called by subclasses to pass their tags to NDMaterialModel
    AcousticMedium (int tag, int classTag, double k, double rho, double gamma = 0.0);

    // Called by clients
    AcousticMedium (int tag, double k, double rho, double gamma = 0.0);

    // For parallel processing
    AcousticMedium (void);

    virtual ~AcousticMedium (void);

    virtual const char *getClassType(void) const {return "AcousticMedium";};

    virtual double getRho( ) ;

    virtual int setTrialStrain (const Vector &v);
    virtual int setTrialStrain (const Vector &v, const Vector &r);
    virtual int setTrialStrainIncr (const Vector &v);
    virtual int setTrialStrainIncr (const Vector &v, const Vector &r);
    virtual const Matrix &getTangent (void);
    virtual const Matrix &getInitialTangent (void);
    virtual const Vector &getStress (void);
    virtual const Vector &getStrain (void);
    
    virtual int commitState (void);
    virtual int revertToLastCommit (void);
    virtual int revertToStart (void);
    
    // Create a copy of material parameters AND state variables
    // Called by GenericSectionXD
    virtual NDMaterial *getCopy (void);

    // Create a copy of just the material parameters
    // Called by the continuum elements
    virtual NDMaterial *getCopy (const char *type);

    // Return a string indicating the type of material model
    virtual const char *getType (void) const;

    virtual int getOrder (void) const;
    
    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
   

      Response *setResponse (const char **argv, int argc, OPS_Stream &output);
      int getResponse (int responseID, Information &matInformation);


    void Print(OPS_Stream &s, int flag = 0);

 //   virtual int setParameter(const char **argv, int argc, Parameter &param);
   // virtual int updateParameter(int parameterID, Information &info);

  protected:
    double Kf;	// Bulk modulus
    double rho ; //mass per unit 3D volume
    double Gamma;	// volumetric drag, force per unit volume per velocity

  private:
    static Vector sigma;        // Stress vector
    static Matrix D;            // Elastic constantsVector sigma;
    Vector epsilon;		// Strain vector
    static Matrix DSensitivity;            // Elastic constantsVector sigma;



// ---------------------sensitivity -----------------------
public:
    int            setParameter             (const char **argv, int argc, Parameter &param);
    int            updateParameter          (int parameterID, Information &info);
	int            activateParameter        (int parameterID);
	const Vector & getStressSensitivity     (int gradNumber, bool conditional);
	int            commitSensitivity        (const Vector & strainGradient, int gradNumber, int numGrads);
	const Matrix & getDampTangentSensitivity     (int gradNumber);
	double         getRhoSensitivity        (int gradNumber);
private:	
 
	Matrix *SHVs;
	int parameterID;

};

#endif
