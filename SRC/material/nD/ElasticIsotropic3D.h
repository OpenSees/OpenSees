//#############################################################################
//# COPYRIGHT (C):     :-))                                                   #
//# PROJECT:           Object Oriented Finite Element Program                 #
//# PURPOSE:                                                                  #
//#                                                                           #
//# CLASS:                                                                    #
//#                                                                           #
//# VERSION:                                                                  #
//# LANGUAGE:          C++																																																				#
//# TARGET OS:         DOS || UNIX || . . .                                   #
//# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                            #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                            #
//# CONTACT:           jeremic@ucdavis.edu                                    #
//#                                                                           #
//# DATE:              Aug, Sept, Oct 2000                                    #
//# UPDATE HISTORY:                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//# SHORT EXPLANATION: 																																																							#
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#############################################################################
                                                                        

//$Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropic3D.h,v $
//$Date: 2000-12-18 10:50:41 $
//$Revision: 1.2 $
                                                                        
#ifndef ElasticIsotropic3D_h
#define ElasticIsotropic3D_h

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <Tensor.h>
#include <stresst.h>
#include <straint.h>

class ElasticIsotropic3D : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropic3D (int tag, double E, double nu);
    ElasticIsotropic3D ();
    ~ElasticIsotropic3D ();

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
    const Tensor &getStressTensor (void);
    const Tensor &getStrainTensor (void);

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker);    
    
    void Print(ostream &s, int flag =0);

  //Private functions
  private:
    void setElasticStiffness(void);


  protected:

  private:
    Vector sigma;		// Stress vector
    Matrix D;			// Elastic constants
    Vector epsilon;		// Strain vector

    stresstensor Stress;	// Stress tensor    
    Tensor Dt;			// Elastic constants tensor
    straintensor Strain;	// Strain tensor    
	     
};


#endif


