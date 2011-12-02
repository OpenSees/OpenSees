// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
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
#include "MultiYieldSurface.h"
#include <Matrix.h>
#include <Tensor.h>

class Response;

class FluidSolidPorousMaterial : public NDMaterial
{
  public:
     // Initialization constructor
     FluidSolidPorousMaterial (int tag, int nd, NDMaterial * soilMat,
			       double combinedBulkModul, double atm);

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
        
     // Calculates the corresponding stress increment (rate), for a given strain increment. 
     const Vector &getStress (void);
     const Vector &getStrain (void);
     const Vector &getCommittedStress (void);
     const Vector &getCommittedStrain (void);
     double getCommittedPressure (void);

     // Accepts the current trial strain values as being on the solution path, and updates 
     // all model parameters related to stress/strain states. Return 0 on success.
     int commitState (void);

     // Revert the stress/strain states to the last committed states. Return 0 on success.
     int revertToLastCommit (void);

     int revertToStart(void) {return 0;}

     // Return an exact copy of itself.
     NDMaterial *getCopy (void);

     // Return a copy of itself if "code"="FluidSolidPorous", otherwise return null.
     NDMaterial *getCopy (const char *code);

     // Return the string "FluidSolidPorous".
     const char *getType (void) const ;

     // Return 2 (2D material).
     int getOrder (void) const ;

     int sendSelf(int commitTag, Channel &theChannel);  
     int recvSelf(int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);     
     Response *setResponse (char **argv, int argc, Information &matInfo);
     int getResponse (int responseID, Information &matInformation);
     void Print(ostream &s, int flag =0);

     int updateParameter(int responseID, Information &eleInformation);

   protected:

   private:
     int ndm;
     static int loadStage;
     static double AtmoPress;
     static Vector * workV;
     static Matrix * workM;
     double combinedBulkModulus;
     NDMaterial * theSoilMaterial;
     double trialExcessPressure;
     double currentExcessPressure;
     double trialVolumeStrain;
     double currentVolumeStrain;
};

#endif
