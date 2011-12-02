//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Pressure dependent elastic isotropic material implementation:
//# CLASS:             PressureDependentElastic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              07July2001
//# UPDATE HISTORY:    22Nov2002 small fixes, formating...
//#
//#
//===============================================================================


#ifndef PressureDependentElastic3D_h
#define PressureDependentElastic3D_h

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <straint.h>
#include <stresst.h>
#include <Tensor.h>
#include <Channel.h>


class PressureDependentElastic3D : public ElasticIsotropicMaterial
{
  public:
    PressureDependentElastic3D (int tag,
                                double E,
                                double nu,
                                double rhop,
                                double expp = 0.6,
                                double pr = 100.0,
                                double pop = 0.5);
    PressureDependentElastic3D ();
    ~PressureDependentElastic3D ();

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent (void);
    const Vector &getStress (void);
    const Vector &getStrain (void);

    int setTrialStrain (const Tensor &v);
    int setTrialStrain (const Tensor &v, const Tensor &r);
    int setTrialStrainIncr (const Tensor &v);
    int setTrialStrainIncr (const Tensor &v, const Tensor &r);
    const Tensor &getTangentTensor (void);
    const stresstensor getStressTensor (void);
    const straintensor getStrainTensor (void);
    const straintensor getPlasticStrainTensor (void);

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);
    void ComputeElasticStiffness(void);

  //Private functions
  private:


  protected:

  private:
    static Vector sigma;  // Stress vector
    static Matrix D;    // Elastic constants
    Vector epsilon;    // Strain vector

    double exp;                 // exponent usually 0.6
    double p_ref;               // Reference pressure, usually atmosphere pressure, i.e. 100kPa
    double p_cutoff;            // Cutoff pressure of this material point

    stresstensor Stress;  // Stress tensor
    Tensor Dt;      // Elastic constants tensor
    straintensor Strain;  // Strain tensor

};
#endif
