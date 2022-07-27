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
** ****************************************************************** */

#ifndef VonPapaDamage_h
#define VonPapaDamage_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <NDMaterial.h>

#include <Information.h>
#include <Parameter.h>
#include <Response.h>

class VonPapaDamage : public NDMaterial {

//-------------------Declarations-------------------------------

public :

  //null constructor
  VonPapaDamage( ) ;

  //full constructor
  VonPapaDamage(int tag,
                double E1,
                double E2,
                double nu12,
                double nu21,
                double G12,
                double rho,
                double Xt,
                double  Xc,
                double  Yt,
                double  Yc,
                double  S,
                double  c1,
                double  c2,
                double  c3,
                double  c4,
                double  c5,
                double  c6,
                double  c7,
                double  c8,
                double  c9,
                double  b) ;


  //destructor
  ~VonPapaDamage( ) ;

  const char *getClassType(void) const {return "VonPapaDamage";};

  //make a clone of this material
  NDMaterial* getCopy( ) ;

  //send back type of material
  const char* getType( ) const ;

  //send back order of strain in vector form
  int getOrder( ) const ;

  //mass per unit volume
  double getRho();

  //get the strain and integrate plasticity equations
  int setTrialStrain( const Vector &strain_from_element) ;

  //unused trial strain functions
  int setTrialStrain( const Vector &v, const Vector &r ) ;
  int setTrialStrainIncr( const Vector &v ) ;
  int setTrialStrainIncr( const Vector &v, const Vector &r ) ;

  //send back the strain
  const Vector& getStrain( ) ;

  //send back the stress
  const Vector& getStress( ) ;

  //send back the tangent
  const Matrix& getTangent( ) ;
  const Matrix& getInitialTangent( ) ;

  //swap history variables
  int commitState( ) ;
  int revertToLastCommit( ) ;
  int revertToStart( ) ;

  //sending and receiving
  int sendSelf(int commitTag, Channel &theChannel) ;
  int recvSelf(int commitTag, Channel &theChannel,
               FEM_ObjectBroker &theBroker ) ;
  //print out material data
  void Print(OPS_Stream &s, int flag = 0) ;

  // Have to implement these
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  // int activateParameter(int parameterID);
  // int setVariable(const char *variable, Information &);
  // int getVariable(const char *variable, Information &);


  Response *setResponse (const char **argv, int argc, OPS_Stream &output);
  int getResponse (int responseID, Information &matInformation);


private :

  const Vector& getDamageState() const;
  void advanceDamageState(int Ncycles);
  void calculateDerDamage(double maxdamage);
  void computeNJUMP(double maxdamage);
  ID getNJUMP(double maxdamage);
  // ID getNJUMP();
  void resetMaxStress();

  //static vectors and matrices
  Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

  //Parametros
  double E1, E2, nu12, nu21, G12, rho;
  double Xt, Xc, Yt, Yc, S;
  double c1, c2, c3, c4, c5, c6, c7, c8, c9, b;

  // Parametros internos
  double D11, D22, D12;
  double deltaSigma1_t, deltaSigma1_c, deltaSigma2_t, deltaSigma2_c, deltaSigma12_t, deltaSigma12_c;
  double bigSigma1t, bigSigma1c, bigSigma2t, bigSigma2c, bigSigma12t, bigSigma12c;
  double strain1_p, strain2_p;

  double dft, dfc, dmt, dmc, dst, dsc;
  double ddft, ddfc, ddmt, ddmc, ddst, ddsc, dstrain1_p, dstrain2_p ;

  static int NVonPapaMaterials;
  static int * NJUMPVEC;
  static int i_current_material_point;
  static int NJUMP;

  ID NJUMP_local;


  double proot_quadraticequ(double a, double b, double c);

  // epsilon plastico

  //index mapping special for plane stress because of
  // condensation on tangent
  void index_map( int matrix_index, int &i, int &j ) ;

} ; //end of VonPapaDamage declarations


#endif
