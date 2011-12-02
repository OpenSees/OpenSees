// $Revision: 1.9 $
// $Date: 2003-02-25 23:33:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/FluidSolidPorousMaterial.h,v $
                                                                        
// Written: ZHY
// Created: August 2000
// Revision: A
//
// Description: This file contains the class prototype for FluidSolidPorousMaterial.
//
// What: "@(#) FluidSolidPorousMaterial.h, revA"

#ifndef FluidSolidPorousMaterial_h
#define FluidSolidPorousMaterial_h

#include <NDMaterial.h>
#include <MultiYieldSurface.h>
#include <Matrix.h>
#include <Tensor.h>

class Response;

class FluidSolidPorousMaterial : public NDMaterial
{
  public:
     // Initialization constructor
     FluidSolidPorousMaterial (int tag, int nd, NDMaterial &soilMat,
			       double combinedBulkModul);

     // Default constructor
     FluidSolidPorousMaterial ();

     // Copy constructor
     FluidSolidPorousMaterial (const FluidSolidPorousMaterial &);

     // Destructor: clean up memory storage space.
     virtual ~FluidSolidPorousMaterial ();

     // Sets the values of the trial strain tensor.
     int setTrialStrain (const Vector &strain);

     // Sets the values of the trial strain and strain rate tensors.
     int setTrialStrain(const Vector &v, const Vector &r);

     int setTrialStrainIncr(const Vector &v);
     int setTrialStrainIncr(const Vector &v, const Vector &r);

     // Calculates current tangent stiffness.
     const Matrix &getTangent (void);
     const Matrix &getInitialTangent (void);

     double getRho(void);

     // Calculates the corresponding stress increment (rate), for a given strain increment. 
     const Vector &getStress (void);
     const Vector &getStrain (void);
     const Vector &getCommittedStress (void);
     const Vector &getCommittedStrain (void);
     const Vector &getCommittedPressure (void);

     // Accepts the current trial strain values as being on the solution path, and updates 
     // all model parameters related to stress/strain states. Return 0 on success.
     int commitState (void);

     // Revert the stress/strain states to the last committed states. Return 0 on success.
     int revertToLastCommit (void);

     int revertToStart(void);

     // Return an exact copy of itself.
     NDMaterial *getCopy (void);

     // Return a copy of itself if "code"="PlainStrain" or "ThreeDimensional", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "PlaneStrain" or "ThreeDimensional".
     const char *getType (void) const ;

     // Return 3 or 6.
     int getOrder (void) const ;

     int sendSelf(int commitTag, Channel &theChannel);  
     int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);     
     Response *setResponse (const char **argv, int argc, Information &matInfo);
     int getResponse (int responseID, Information &matInformation);
     void Print(OPS_Stream &s, int flag =0);

     int updateParameter(int responseID, Information &eleInformation);

   protected:

   private:
     int ndm;
     static int loadStage;
     double combinedBulkModulus;
     NDMaterial * theSoilMaterial;
     double trialExcessPressure;
     double currentExcessPressure;
     double trialVolumeStrain;
     double currentVolumeStrain;

     static Vector workV3;
     static Vector workV6;
     static Matrix workM3;
     static Matrix workM6;
};

#endif
