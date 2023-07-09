// $Revision: 1.16 $
// $Date: 2010-06-12 18:07:00 $
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
#include <Matrix.h>
#include "soil/T2Vector.h"
class MultiYieldSurface;

class Response;

class FluidSolidPorousMaterial : public NDMaterial
{
  public:
     // Initialization constructor
     FluidSolidPorousMaterial (int tag, int nd, NDMaterial &soilMat,
			       double combinedBulkModul, double atm=101.);

     // Default constructor
     FluidSolidPorousMaterial ();

     // Copy constructor
     FluidSolidPorousMaterial (const FluidSolidPorousMaterial &);

     // Destructor: clean up memory storage space.
     virtual ~FluidSolidPorousMaterial ();

     const char *getClassType(void) const {return "FluidSolidPorousMaterial";};

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
     Response *setResponse (const char **argv, int argc, OPS_Stream &s);
     int getResponse (int responseID, Information &matInformation);
     void Print(OPS_Stream &s, int flag =0);

     int setParameter(const char **argv, int argc, Parameter &param);
     int updateParameter(int responseID, Information &eleInformation);

     // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

   protected:

   private:
     static int matCount;
     static int* ndmx;
     static int* loadStagex;
     static double* combinedBulkModulusx;
     static double pAtm;
     int matN;
     NDMaterial * theSoilMaterial;
     double trialExcessPressure;
     double currentExcessPressure;
     double trialVolumeStrain;
     double currentVolumeStrain;
     double initMaxPress;
     int e2p;

     Vector theSoilCommittedStress;
     Vector theSoilCommittedStrain;

     static Vector workV3;
     static Vector workV6;
     static Matrix workM3;
     static Matrix workM6;
};

#endif

