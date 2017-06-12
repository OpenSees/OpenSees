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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $


#ifndef PlasticDamageConcretePlaneStressThermal_h
#define PlasticDamageConcretePlaneStressThermal_h

 // Based on PlasticDamageConcretePlaneStress
 
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// Type Definitions
typedef struct {
  double E;
  double nu;
  double ft;
  double fc;
  double Ap;
  double An;
  double Bn;
  double beta;
} structMatData;

typedef struct {
  double sig[3];
  double eps_p[4];
  double rn;
  double rp;
  double dp;
  double dn;
} structCurrent;

typedef struct {
  double eps_p[4];
  double rn;
  double rp;
  double dp;
  double dn;
} structCommitted;

//typedef struct {
//  double sig[3];
//  double km[9];
//  struct2_T Pres;
//  double eps[3];
//  struct3_T Past;
//  double Deps[3];
//} struct1_T;

class PlasticDamageConcretePlaneStressThermal : public NDMaterial
{
  public:
  PlasticDamageConcretePlaneStressThermal(int tag, 
				   double E, 
				   double nu, 
				   double ft,
				   double fc, 
				   double beta = 0.6, 
				   double Ap = 0.5, 
				   double An = 2.0, 
				   double Bn = 0.75);
  PlasticDamageConcretePlaneStressThermal();
  ~PlasticDamageConcretePlaneStressThermal();
  
  const char *getClassType(void) const {return "PlasticDamageConcrete3d";};
  
  int setTrialStrain (const Vector &v);
  int setTrialStrain (const Vector &v, const Vector &r);
  int setTrialStrainIncr (const Vector &v);
  int setTrialStrainIncr (const Vector &v, const Vector &r);
  const Matrix &getTangent (void);
  const Matrix &getInitialTangent (void);
  
  const Vector &getStress (void);
  const Vector &getStrain (void);
  
  //Temperature
  double setThermalTangentAndElongation(double &TempT, double &, double &);//L. Jiang [SIF]

  const Vector& getTempAndElong( void);  ///L. Jiang [SIF]
  
  int commitState (void);
  int revertToLastCommit (void);
  int revertToStart (void);
  
  NDMaterial*getCopy(const char *type);
  NDMaterial *getCopy (void);
  const char *getType (void) const;
  int getOrder (void) const;
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  // parameters
  double E;     // elastic modulus
  double nu;    // Poisson ratio 
  double ft;    // tensile yield strength
  double fc;    // compressive yield strength
  double beta;  // plastic deformation rate
  double Ap;    // damage parameter
  double An;    // damage parameter
  double Bn;    // damage parameter
  
  double E0,fc0,ft0;
  double sig[3];
  double eps[3];
  double Deps[3];

  Matrix Ce;
  Matrix Ce0;
  Matrix CeCommitted;

  //  double Ce[9];  // current tangent
  //  double Ce0[9]; // initial tangent
  //  double CeCommitted[9]; // committed tangent
  
  // current state variables
  //double sig[3];
  double eps_p[4];
  double rn;
  double rp;
  double dp;
  double dn;
  
  // committed state variables
  double Committed_sig[4];
  double Committed_eps[3];
  double Committed_eps_p[4];
  double Committed_rn;
  double Committed_rp;
  double Committed_dp;
  double Committed_dn;

  Vector stress;
  Vector strain;
  Vector Cstress;
  Vector Cstrain;
  Vector TempAndElong;
};

#endif
