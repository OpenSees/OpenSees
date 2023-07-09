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
//# PURPOSE:           Elastic Cross Anisotropic Material implementation:
//# CLASS:             ElasticIsotropic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Yi Bian, Boris Jeremic
//#
//#
//# DATE:              10Oct2002
//# UPDATE HISTORY:    March 20, 2003 Re-activated Joey Yang
//#                    Aug2006   Z.Cheng
//#
//===============================================================================

#ifndef ElasticCrossAnisotropic_h
#define ElasticCrossAnisotropic_h

//#include <Channel.h>
//#include <string.h>
//#include <G3Globals.h>

#include <Matrix.h>
#include <Vector.h>
#include <Tensor.h>
#include <stresst.h>
#include <straint.h>

#include <NDMaterial.h>

class ElasticCrossAnisotropic : public NDMaterial
{
public:
  ElasticCrossAnisotropic(int tag,
                          double Ehp,
			  double Evp,
			  double nuhvp,
			  double nuhhp,
			  double Ghvp,
			  double rhop = 0.0);
  ElasticCrossAnisotropic ();
  ~ElasticCrossAnisotropic ();

        const char *getClassType(void) const {return "ElasticCrossAnisotropic";};

        double getrho ();
	double getMatParameter(int MatParameterID);

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
        NDMaterial *getCopy (const char *type);
        const char *getType (void) const;

        void Print(OPS_Stream &s, int flag = 0);

        int sendSelf(int commitTag, Channel &theChannel);
        int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

private:
    static stresstensor Stress;   // Stress tensor
    static Tensor Dt;         // Elastic constants tensor
    straintensor Strain;   // Strain tensor

// all the directions are relative so we call them "horizontal" and "vertical", take that
// horizontal is one plane of anisotropy while vertical is the axes perpendicular to that plane.
  double Eh;             // Eh: Young's modulus in any horizontal direction.
  double Ev;             // Ev: Young's modulus in a vertical direction.
  double nuhv;           // nuhv: Poisson's ratio for strain in the vertical direction due to a horizontal direct stress.
  double nuhh;           // nvhh: Poisoon's ratio for strain in any horizontal direction due to a horizontal direct stress at right angles.
  double Ghv;            // Ghv: Modulus of shear deformation in a vertical plane.
  double rho;            // Mass density
};

#endif

