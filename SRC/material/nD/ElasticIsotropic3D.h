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
//# PURPOSE:           Elastic Isotropic Material implementation:
//# CLASS:             ElasticIsotropic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              10Oct2000
//# UPDATE HISTORY:    22Nov2002 small fixes 
//#                    Aug2006   Z.Cheng
//#
//===============================================================================
                                                                        
                                                                        
#ifndef ElasticIsotropic3D_h
#define ElasticIsotropic3D_h
	 
#include <ElasticIsotropicMaterial.h>

class ElasticIsotropic3D : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropic3D (int tag, double E, double nu, double rho);
    ElasticIsotropic3D ();
    ~ElasticIsotropic3D ();

    const char *getClassType(void) const {return "ElasticIsotropic3D";};

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
    const stresstensor& getStressTensor (void);
    const straintensor& getStrainTensor (void);
 
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    void Print(OPS_Stream &s, int flag =0);
    
  private:  
    void setInitElasticStiffness(void);

  protected:

  private:
    static Vector sigma;        // Stress vector
    static Matrix D;            // Elastic constantsVector sigma;
    Vector epsilon;		// Strain vector

    static stresstensor Stress;	// Stress tensor    
    static Tensor Dt;           // Elastic constants tensor
    straintensor Strain;	// Strain tensor    
};


#endif


