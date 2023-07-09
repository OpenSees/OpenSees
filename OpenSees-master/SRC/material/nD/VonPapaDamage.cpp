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

// $Revision: 1.7 $
// $Date: 2008-10-20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/VonPapaDamage.cpp,v $

// Written: Ed "C++" Love
//
// VonPapaDamage isotropic hardening material class
//
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi)
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q
//
//  Linear Viscosity
//  gamma = phi / eta  ( if phi > 0 )
//
//  Backward Euler Integration Routine
//  Yield condition enforced at time n+1
//
//  Send strains in following format :
//
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2
//
//  set eta := 0 for rate independent case
//

#include <VonPapaDamage.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Response.h>
#include <MaterialResponse.h>



// Includes for general things (Felipe Elgueta)
#include <algorithm>  // for "max" function
using namespace::std;
#include <stdlib.h>
// #include <cstdlib>
//#include <bits/stdc++.h>
#include <limits.h>

// Vector VonPapaDamage :: strain_vec(3) ;
Vector VonPapaDamage :: stress_vec(3) ;
Matrix VonPapaDamage :: tangent_matrix(3, 3) ;

int VonPapaDamage::NVonPapaMaterials = 0;
int VonPapaDamage::i_current_material_point = 0;
int* VonPapaDamage::NJUMPVEC = 0;
int VonPapaDamage::NJUMP = 0;




void* OPS_VonPapaDamage(void)
{

  opserr << "ndMaterial VonPapaDamage tag E1, E2, nu12, nu21, G12, rho, Xt ,Xc ,Yt ,Yc ,S  ,c1 ,c2 ,c3 ,c4 ,c5 ,c6 ,c7 ,c8 ,c9 , b \n";

  int tag;
  double E1, E2, nu12, nu21, G12, rho;
  double Xt, Xc, Yt, Yc, S;
  double c1, c2, c3, c4, c5, c6, c7, c8, c9, b;


  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 22) {
    opserr << "ndMaterial VonPapaDamage tag E1, E2, nu12, nu21, G12, rho, Xt ,Xc ,Yt ,Yc ,S  ,c1 ,c2 ,c3 ,c4 ,c5 ,c6 ,c7 ,c8 ,c9 , b \n";
    return 0;
  }

  int iData[1];
  double dData[21];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial VonPapaDamage \n";
    return 0;
  }
  tag = iData[0];

  numData = 21;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid double values: nDMaterial VonPapaDamage " << tag << endln;
    return 0;
  }
  E1 = dData[0];
  E2 = dData[1];
  nu12 = dData[2];
  nu21 = dData[3];
  G12 = dData[4];
  rho = dData[5];

  Xt = dData[6];
  Xc = dData[7];
  Yt = dData[8];
  Yc = dData[9];
  S  = dData[10];

  c1 = dData[11];
  c2 = dData[12];
  c3 = dData[13];
  c4 = dData[14];
  c5 = dData[15];
  c6 = dData[16];
  c7 = dData[17];
  c8 = dData[18];
  c9 = dData[19];
  b  = dData[20];

  opserr << "Creating new VonPapaDamage with \n"
         << "tag  = " << tag  << endln
         << "E1   = " << E1   << endln
         << "E2   = " << E2   << endln
         << "nu12 = " << nu12 << endln
         << "nu21 = " << nu21 << endln
         << "G12  = " << G12  << endln
         << "rho  = " << rho  << endln
         << "Xt   = " << Xt   << endln
         << "Xc   = " << Xc   << endln
         << "Yt   = " << Yt   << endln
         << "Yc   = " << Yc   << endln
         << "S    = " << S    << endln
         << "c1   = " << c1   << endln
         << "c2   = " << c2   << endln
         << "c3   = " << c3   << endln
         << "c4   = " << c4   << endln
         << "c5   = " << c5   << endln
         << "c6   = " << c6   << endln
         << "c7   = " << c7   << endln
         << "c8   = " << c8   << endln
         << "c9   = " << c9   << endln
         << "b    = " << b    << endln << endln ;




  NDMaterial *theMaterial =  new VonPapaDamage (tag, E1, E2, nu12, nu21, G12, rho, Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b);

  return theMaterial;
}






//null constructor
VonPapaDamage ::  VonPapaDamage( ) :
  NDMaterial(0, ND_TAG_VonPapaDamage),
  strain_vec(3),
  E1(0),
  E2(0),
  nu12(0),
  nu21(0),
  G12(0),
  rho(0),
  D11(0),
  D22(0),
  D12(0),
  NJUMP_local(0)
{
  NVonPapaMaterials++;
}


//full constructor
VonPapaDamage ::
VonPapaDamage(int tag,
              double E1_,
              double E2_,
              double nu12_,
              double nu21_,
              double G12_,
              double rho_,
              double Xt_,
              double  Xc_,
              double  Yt_,
              double  Yc_,
              double  S_,
              double  c1_,
              double  c2_,
              double  c3_,
              double  c4_,
              double  c5_,
              double  c6_,
              double  c7_,
              double  c8_,
              double  c9_,
              double  b_) :
  NDMaterial(tag, ND_TAG_VonPapaDamage),
  strain_vec(3),
  E1(E1_),
  E2(E2_),
  nu12(nu12_),
  nu21(nu21_),
  G12(G12_),
  rho(rho_),
  Xt(Xt_),
  Xc(Xc_),
  Yt(Yt_),
  Yc(Yc_),
  S(S_),
  c1(c1_),
  c2(c2_),
  c3(c3_),
  c4(c4_),
  c5(c5_),
  c6(c6_),
  c7(c7_),
  c8(c8_),
  c9(c9_),
  b(b_),
  D11(0),
  D22(0),
  D12(0),
  deltaSigma1_t(0),
  deltaSigma1_c(0),
  deltaSigma2_t(0),
  deltaSigma2_c(0),
  deltaSigma12_t(0),
  deltaSigma12_c(0),
  bigSigma1t(0),
  bigSigma1c(0),
  bigSigma2t(0),
  bigSigma2c(0),
  bigSigma12t(0),
  bigSigma12c(0),
  strain1_p(0),
  strain2_p(0),
  dft(0),
  dfc(0),
  dmt(0),
  dmc(0),
  dst(0),
  dsc(0),
  NJUMP_local(3)
{
  D11 = 0;
  D22 = 0;
  D12 = 0;
  deltaSigma1_t = 0;
  deltaSigma1_c = 0;
  deltaSigma2_t = 0;
  deltaSigma2_c = 0;
  deltaSigma12_t = 0;
  deltaSigma12_c = 0;
  bigSigma1t = 0;
  bigSigma1c = 0;
  bigSigma2t = 0;
  bigSigma2c = 0;
  bigSigma12t = 0;
  bigSigma12c = 0;
  strain1_p = 0;
  strain2_p = 0;
  dft = 0;
  dfc = 0;
  dmt = 0;
  dmc = 0;
  dst = 0;
  dsc = 0;

  double ddft, ddfc, ddmt, ddmc, ddst, ddsc, dstrain1_p, dstrain2_p ;

  NVonPapaMaterials++;

}


//destructor
VonPapaDamage :: ~VonPapaDamage( )
{  }


//make a clone of this material
NDMaterial* VonPapaDamage :: getCopy( )
{
  VonPapaDamage  *clone;
  clone = new VonPapaDamage( ) ;   //new instance of this class
  *clone = *this ;                 //asignment to make copy
  return clone ;
}


//send back type of material
const char* VonPapaDamage :: getType( ) const
{
  return "OrthotropicPlaneStress" ;
}


//send back order of strain in vector form
int VonPapaDamage :: getOrder( ) const
{
  return 3 ;
}

//mass per unit volume
double
VonPapaDamage::getRho( )
{
  return rho ;
}

//get the strain and integrate plasticity equations
int VonPapaDamage :: setTrialStrain( const Vector &strain_from_element )
{
  // opserr << "VonPapaDamage :: setTrialStrain called!" << endln;
  // opserr << "strain_from_element = " << strain_from_element << endln;
  strain_vec = strain_from_element;

  // Σ11 = max(abs(sigma11),Σ11);
  // Σ22 = max(abs(sigma22),Σ22);
  // Σ12 = max(abs(sigma12),Σ12);

  // Computing damage varriables
  // advanceDamageState(100);


  return 0 ;
}

void VonPapaDamage :: resetMaxStress()
{

  // opserr << "VonPapaDamage :: resetMaxStress called" << endln;

  // opserr <<'deltaSigma1_t = '<< deltaSigma1_t  << endln;
  // opserr <<'deltaSigma1_c = '<< deltaSigma1_c  << endln;
  // opserr <<'deltaSigma2_t = '<< deltaSigma2_t  << endln;
  // opserr <<'deltaSigma2_c = '<< deltaSigma2_c  << endln;
  // opserr <<'deltaSigma12_t= '<< deltaSigma12_t << endln;
  // opserr <<'deltaSigma12_c= '<< deltaSigma12_c << endln;

  // Resetting max stresses
  deltaSigma1_t  = 0.0;
  deltaSigma1_c  = 0.0;
  deltaSigma2_t  = 0.0;
  deltaSigma2_c  = 0.0;
  deltaSigma12_t = 0.0;
  deltaSigma12_c = 0.0;
  // bigSigma1t = 0;
  // bigSigma1c = 0;
  // bigSigma2t = 0;
  // bigSigma2c = 0;
  // bigSigma12t = 0;
  // bigSigma12c = 0;


  // Resetting permanent strains
  // strain1_p = 0;
  // strain2_p = 0;
}


//unused trial strain functions
int VonPapaDamage :: setTrialStrain( const Vector &v, const Vector &r )
{
  opserr << "VonPapaDamage :: setTrialStrain( const Vector &v, const Vector &r ) -- should not be used! \n";
  return this->setTrialStrain( v ) ;
}

int VonPapaDamage :: setTrialStrainIncr( const Vector &v )
{
  opserr << "VonPapaDamage :: setTrialStrainIncr( const Vector &v ) -- should not be used! \n";
  return -1 ;
}

int VonPapaDamage :: setTrialStrainIncr( const Vector &v, const Vector &r )
{
  opserr << "VonPapaDamage :: setTrialStrainIncr( const Vector &v, const Vector &r ) -- should not be used! \n";
  return this->setTrialStrainIncr(v);
}


//send back the strain
const Vector& VonPapaDamage :: getStrain( )
{
  return strain_vec ;
}


//send back the stress
const Vector& VonPapaDamage :: getStress( )
{

  stress_vec = this->getTangent() * strain_vec;

  return stress_vec ;
}

//send back the tangent
const Matrix& VonPapaDamage :: getTangent( )
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 )
  //


  double den = 1 - nu12 * nu21;

  tangent_matrix(0, 0) = E1 / den * (1 - D11) ;
  tangent_matrix(1, 1) = E2 / den * (1 - D22) ; ;
  tangent_matrix(2, 2) = G12 * (1 - D12) ; ;

  tangent_matrix(0, 1) = nu21 * E1 / den * sqrt(1 - D11) * sqrt(1 - D22) ;
  tangent_matrix(1, 0) = nu12 * E2 / den * sqrt(1 - D11) * sqrt(1 - D22) ;

  tangent_matrix(0, 2) = 0. ;
  tangent_matrix(2, 0) = 0. ;

  tangent_matrix(1, 2) = 0. ;
  tangent_matrix(2, 1) = 0. ;

  return tangent_matrix ;
}


//send back the tangent
const Matrix& VonPapaDamage :: getInitialTangent( )
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 )
  //


  return this->getTangent() ;
}

int
VonPapaDamage::commitState( )
{
  // opserr << "VonPapaDamage::commitState called!" << endln;
  //Calcular los stresses nuevos y guardar los max min...

  // Getting actual strains
  const Vector& epsilon = getStrain();
  // opserr << "epsilon = " << epsilon << endln;

  // Compute the stiffness matrix in the LOCAL MATERIAL coordinate system
  // The elastic constants are updated through the composite fatigue law (ref. Paepegem model)
  const Matrix& Q = getTangent();

  // Getting stresses (should remove the permament strain too)
  static Vector stress(3);   // malloc
  stress.Zero();

  // Compute P.K. II stresses in the current ply in the LOCAL MATERIAL coord. system
  // Gotta be careful with voigt notation in epsilon(2)
  stress(0) = Q(0, 0) * (epsilon(0) - strain1_p) + Q(0, 1) * (epsilon(1) - strain2_p) + Q(0, 2) * epsilon(2);
  stress(1) = Q(1, 0) * (epsilon(0) - strain1_p) + Q(1, 1) * (epsilon(1) - strain2_p) + Q(1, 2) * epsilon(2);
  stress(2) = Q(2, 0) * epsilon(0) +               Q(2, 1) * epsilon(1) +               Q(2, 2) * epsilon(2);

  // if (stress(0) != 0)
  // {
  //   opserr << "stress(0) = " << stress(0) << endln;
  //   opserr << "deltaSigma1_t = " << deltaSigma1_t << endln;
  //   opserr << "deltaSigma1_c = " << deltaSigma1_c << endln << endln;
  // }

  // if (stress(1) != 0)
  // {
  //   opserr << "stress(1) = " << stress(1) << endln;
  //   opserr << "deltaSigma2_t = " << deltaSigma2_t << endln;
  //   opserr << "deltaSigma2_c =" << deltaSigma2_c << endln << endln;
  // }

  // if (stress(2) != 0)
  // {
  //   opserr << "stress(2) = " << stress(2) << endln;
  //   opserr << "deltaSigma12_t = " << deltaSigma12_t << endln;
  //   opserr << "deltaSigma12_c = " << deltaSigma12_c << endln << endln;
  // }


  // // Updating max stresses
  // if (stress(0) >= 0.0)
  // {
  //   deltaSigma1_t = max(stress(0), deltaSigma1_t);
  // } else if (stress(0) < 0.0) {
  //   deltaSigma1_c = min(stress(0), deltaSigma1_c);
  // }

  // if (stress(1) >= 0.0)
  // {
  //   deltaSigma2_t = max(stress(1), deltaSigma2_t);
  // } else if (stress(1) < 0.0) {
  //   deltaSigma2_c = min(stress(1), deltaSigma2_c);
  // }

  // if (stress(2) >= 0.0)
  // {
  //   deltaSigma12_t = max(stress(2), deltaSigma12_t);
  // } else if (stress(2) < 0.0) {
  //   deltaSigma12_c = min(stress(2), deltaSigma12_c);
  // }


  double bigSigma1t_1D, bigSigma1t_2D ,  bigSigma1c_1D , bigSigma1c_2D ,  bigSigma2t_1D , bigSigma2t_2D ,  bigSigma2c_1D , bigSigma2c_2D ,  bigSigma12t_1D , bigSigma12t_2D ,  bigSigma12c_1D , bigSigma12c_2D;
  bigSigma1t_1D = bigSigma1t_2D =  bigSigma1c_1D = bigSigma1c_2D =  bigSigma2t_1D = bigSigma2t_2D =  bigSigma2c_1D = bigSigma2c_2D =  bigSigma12t_1D = bigSigma12t_2D =  bigSigma12c_1D = bigSigma12c_2D = 0.0;




  double aa;
  double bb;
  double cc;
  double proot;

  // D11 = fiber = df
  // D22 = matrix = dm
  // D12 = shear = ds

  // Begin to compute the updated damage variables in tension and in compression (fiber/matrix) and shear

  aa = 1.0 / (Xt * Xc) * pow( stress(0) / (1.0 - D11), 2.0 );
  bb = (1.0 / Xt - 1.0 / Xc) * (stress(0) / (1.0 - D11));
  cc = (1.0 / Yt - 1.0 / Yc) * (stress(1) / (1.0 - D22)) + 1.0 / (Yt * Yc) * pow(stress(1) / (1.0 - D22), 2.0) + 1.0 / pow(S, 2.0) * pow( stress(2) / (1.0 - D12), 2.0) - 1.0;

  proot = proot_quadraticequ(aa, bb, cc);
  // opserr << "proot = " << proot << endln;

  if (stress(0) > deltaSigma1_t)
  {
    deltaSigma1_t = stress(0);
    if (D11 < 1.0)
    {
      bigSigma1t_2D = 1.0 / proot;
      bigSigma1t_1D = stress(0) / ( (1.0 - D11) * Xt );
    }
    else
    {
      bigSigma1t_2D = 1.0;
      bigSigma1t_1D = 1.0;
    }
    bigSigma1t = bigSigma1t_2D / (1.0 + bigSigma1t_2D - bigSigma1t_1D);
    // opserr << "bigSigma1t_1D = " << bigSigma1t_1D << endln;
    // opserr << "bigSigma1t_2D = " << bigSigma1t_2D << endln;
    // opserr << "bigSigma1t = " << bigSigma1t << endln << endln;
  }

  if (stress(0) < deltaSigma1_c)
  {
    deltaSigma1_c = stress(0);
    if (D11 < 1.0)
    {
      bigSigma1c_2D = 1.0 / proot;
      bigSigma1c_1D = -stress(0) / ((1.0 - D11) * Xc);
    }
    else
    {
      bigSigma1c_2D = 1.0;
      bigSigma1c_1D = 1.0;
    }
    bigSigma1c = bigSigma1c_2D / (1.0 + bigSigma1c_2D - bigSigma1c_1D);
    // opserr << "bigSigma1c_1D = " << bigSigma1c_1D << endln;
    // opserr << "bigSigma1c_2D = " << bigSigma1c_2D << endln;
    // opserr << "bigSigma1c = " << bigSigma1c << endln << endln;
  }

  // opserr << "bigSigma1t_1D = " << bigSigma1t_1D << endln;
  // opserr << "bigSigma1c_1D = " << bigSigma1c_1D << endln;
  // opserr << "bigSigma1t_2D = " << bigSigma1t_2D << endln;
  // opserr << "bigSigma1c_2D = " << bigSigma1c_2D << endln;


  // opserr << "bigSigma1t = " << bigSigma1t << endln;
  // opserr << "bigSigma1c = " << bigSigma1c << endln;

  aa = 1.0 / (Yt * Yc) * pow(stress(1) / (1.0 - D22), 2.0);
  bb = (1.0 / Yt - 1.0 / Yc) * (stress(1) / (1.0 - D22));
  cc = (1.0 / Xt - 1.0 / Xc) * (stress(0) / (1.0 - D11)) + 1.0 / (Xt * Xc) * pow(stress(0) / (1.0 - D11), 2.0) + 1.0 / pow(S, 2.0) * pow(stress(2) / (1.0 - D12), 2.0) - 1.0;

  proot = proot_quadraticequ(aa, bb, cc);

  if (stress(1) > deltaSigma2_t)
  {
    deltaSigma2_t = stress(1);
    if (D22 < 1.0)
    {
      bigSigma2t_2D = 1.0 / proot;
      bigSigma2t_1D = stress(1) / ((1.0 - D22) * Yt);
    }
    else
    {
      bigSigma2t_2D = 1.0;
      bigSigma2t_1D = 1.0;
    }
    bigSigma2t = bigSigma2t_2D / (1.0 + bigSigma2t_2D - bigSigma2t_1D);
  }

  if (stress(1) < deltaSigma2_c)
  {
    deltaSigma2_c = stress(1);
    if (D22 < 1.0)
    {
      bigSigma2c_2D = 1.0 / proot;
      bigSigma2c_1D = -stress(1) / ((1.0 - D22) * Yc);
    }
    else
    {
      bigSigma2c_2D = 1.0;
      bigSigma2c_1D = 1.0;
    }
    bigSigma2c = bigSigma2c_2D / (1.0 + bigSigma2c_2D - bigSigma2c_1D);
  }




  // opserr << "bigSigma2t_1D = " << bigSigma2t_1D << endln;
  // opserr << "bigSigma2c_1D = " << bigSigma2c_1D << endln;
  // opserr << "bigSigma2t_2D = " << bigSigma2t_2D << endln;
  // opserr << "bigSigma2c_2D = " << bigSigma2c_2D << endln;

  // opserr << "bigSigma2t = " << bigSigma2t << endln;
  // opserr << "bigSigma2c = " << bigSigma2c << endln;

  aa = 1.0 / pow(S, 2) * pow(stress(2) / (1.0 - D12), 2.0);
  bb = 0.0;
  cc = (1.0 / Xt - 1.0 / Xc) * stress(0) / (1.0 - D11) + (1.0 / Yt - 1.0 / Yc) * stress(1) / (1.0 - D22) + 1.0 / (Xt * Xc) * pow(stress(0) / (1.0 - D11), 2.0) + 1.0 / (Yt * Yc) * pow(stress(1) / (1.0 - D22), 2.0) - 1.0;

  proot = proot_quadraticequ(aa, bb, cc);

  if (stress(2) > deltaSigma12_t)
  {
    deltaSigma12_t = stress(2);
    if (D12 < 1.0)
    {
      bigSigma12t_2D = 1.0 / proot;
      bigSigma12t_1D = abs(stress(2)) / ((1.0 - D12) * S);
    }
    else
    {
      bigSigma12t_2D = 1.0;
      bigSigma12t_1D = 1.0;
    }
    bigSigma12t = bigSigma12t_2D / (1.0 + bigSigma12t_2D - bigSigma12t_1D);
  }

  if (stress(2) < deltaSigma12_c)
  {
    deltaSigma12_c = stress(2);
    if (D12 < 1.0)
    {
      bigSigma12c_2D = 1.0 / proot;
      bigSigma12c_1D = abs(stress(2)) / ((1.0 - D12) * S);
    }
    else
    {
      bigSigma12c_2D = 1.0;
      bigSigma12c_1D = 1.0;
    }
    bigSigma12c = bigSigma12c_2D / (1.0 + bigSigma12c_2D - bigSigma12c_1D);
  }

  // Doing bound checks, not sure if should do
  if (bigSigma1t >= 1.0) bigSigma1t = 0.999;
  if (bigSigma1c >= 1.0) bigSigma1c = 0.999;
  if (bigSigma2t >= 1.0) bigSigma2t = 0.999;
  if (bigSigma2c >= 1.0) bigSigma2c = 0.999;
  if (bigSigma12t >= 1.0) bigSigma12t = 0.999;
  if (bigSigma12c >= 1.0) bigSigma12c = 0.999;

  // opserr << "bigSigma12t_1D = " << bigSigma12t_1D << endln;
  // opserr << "bigSigma12c_1D = " << bigSigma12c_1D << endln;
  // opserr << "bigSigma12t_2D = " << bigSigma12t_2D << endln;
  // opserr << "bigSigma12c_2D = " << bigSigma12c_2D << endln;

  // opserr << "bigSigma12t = " << bigSigma12t << endln;
  // opserr << "bigSigm12c = " << bigSigma12c << endln << endln;

  // opserr << "bigSigma1t = " << bigSigma1t << endln;
  // opserr << "bigSigma1c = " << bigSigma1c << endln;
  // opserr << "bigSigma2t = " << bigSigma2t << endln;
  // opserr << "bigSigma2c = " << bigSigma2c << endln;
  // opserr << "bigSigma12t = " <<  bigSigma12t << endln;
  // opserr << "bigSigma12c = " <<  bigSigma12c << endln;




  return 0;
}

int
VonPapaDamage::revertToLastCommit( ) {


  return 0;
}


int
VonPapaDamage::revertToStart( ) {


  return 0;
}

int
VonPapaDamage::sendSelf(int commitTag, Channel &theChannel)
{


  return -1;
}

int
VonPapaDamage::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{


  return -1;
}


//print out material data
void VonPapaDamage :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "VonPapaDamage : " ;
  s << this->getType( ) << endln ;
  s << "Elastic Modulus 1 =   " << E1        << endln ;
  s << "Elastic Modulus 2 =   " << E2        << endln ;
  s << "Poisson's ratio 12=  " << nu12       << endln ;
  s << "Poisson's ratio 21=  " << nu21       << endln ;
  s << "Shear constant G12=  " << G12       << endln ;
  s << "mass density =        " << rho     << endln ;
  s << endln ;
}


//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3
void
VonPapaDamage :: index_map( int matrix_index, int &i, int &j )
{
  switch ( matrix_index + 1 ) { //add 1 for standard tensor indices

  case 1 :
    i = 1 ;
    j = 1 ;
    break ;

  case 2 :
    i = 2 ;
    j = 2 ;
    break ;

  case 3 :
    i = 1 ;
    j = 2 ;
    break ;

  case 4 :
    i = 3 ;
    j = 3 ;
    break ;

  case 5 :
    i = 2 ;
    j = 3 ;
    break ;

  case 6 :
    i = 3 ;
    j = 1 ;
    break ;


  default :
    i = 1 ;
    j = 1 ;
    break ;

  } //end switch

  i-- ; //subtract 1 for C-indexing
  j-- ;

  return ;
}

void VonPapaDamage :: calculateDerDamage(double max_incr_dam_var)
{
  // Calcular dD11_dN dD22_dN dD12_dN;

  // opserr << "bigSigma1t = " << bigSigma1t << endln;
  // opserr << "bigSigma1c = " << bigSigma1c << endln;
  // opserr << "bigSigma2t = " << bigSigma2t << endln;
  // opserr << "bigSigma2c = " << bigSigma2c << endln;
  // opserr << "bigSigma12t = " << bigSigma12t << endln;
  // opserr << "bigSigm12c = " << bigSigma12c << endln << endln;

  // opserr << "deltaSigma1_t = "<< deltaSigma1_t << endln;
  // opserr << "deltaSigma1_c = "<< deltaSigma1_c << endln;
  // opserr << "deltaSigma2_t = "<< deltaSigma2_t << endln;
  // opserr << "deltaSigma2_c = "<< deltaSigma2_c << endln;
  // opserr << "deltaSigma12_t = "<< deltaSigma12_t << endln;
  // opserr << "deltaSigma12_c = "<< deltaSigma12_c << endln << endln;

  ddft = ddfc = ddmt = ddmc = ddst = ddsc = 0.0;

  // Njump for damage in each damage type
  // int NJUMP_threshold = 5e4;
  // int NJUMP_threshold = 1e6; // Gotta try this
  int NJUMP_threshold = INT_MAX; // Maybe this is too much

  int Njump_dft, Njump_dfc, Njump_dmt, Njump_dmc, Njump_dst, Njump_dsc;
  Njump_dft = Njump_dfc = Njump_dmt = Njump_dmc = Njump_dst = Njump_dsc = 0 ;


  double eps = 1e-4;
  double tmp = 0;

  // dD11/dN in case of positive Sigma11 stress
  if (deltaSigma1_t > 0.0 && bigSigma1t > 0.0)
  {
    if (dft < 1.0)
    {
      tmp = - c2 * dft / (sqrt(bigSigma1t) * (1.0 + pow(dfc, 2) + pow(D12, 2.0)));
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddft = c1 * (1.0 + pow(dfc, 2.0) + pow(D12, 2.0)) * bigSigma1t * exp(tmp) + c3 * dft * pow(bigSigma1t, 2.0) * (1.0 + sqrt(dfc) * exp(c8 * sqrt(dfc))  / (1.0 + exp(-c5 * (bigSigma1t - c7))) ) * (1.0 + exp(c5 * (bigSigma1t - c4)));

      if (max_incr_dam_var / ddft > NJUMP_threshold) {
        Njump_dft = NJUMP_threshold;
      } else {
        Njump_dft = max_incr_dam_var / ddft;
      }
    } else {
      Njump_dft = NJUMP_threshold;
    }
    // opserr << "tmp ddft = " << tmp << endln;
    // opserr << "dft =" << dft << endln;
    // opserr << "dfc =" << dfc << endln;
    // opserr << "D12 =" << D12 << endln;
    // opserr << "bigSigma1t =" << bigSigma1t << endln;

    // opserr << "c1 * (1.0 + pow(dfc, 2) + pow(D12, 2))* bigSigma1t * exp(tmp) =" << c1 * (1.0 + pow(dfc, 2) + pow(D12, 2))* bigSigma1t * exp(tmp) << endln;
    // opserr << "c3 * dft * pow(bigSigma1t, 2) * (1.0 + (sqrt(dfc) * exp(c8 * sqrt(dfc)))  / (1.0 + exp(-c5 * (bigSigma1t - c7))) ) * (1.0 + exp(c5 * (bigSigma1t - c4))) = "<< c3 * dft * pow(bigSigma1t, 2) * (1.0 + (sqrt(dfc) * exp(c8 * sqrt(dfc)))  / (1.0 + exp(-c5 * (bigSigma1t - c7))) ) * (1.0 + exp(c5 * (bigSigma1t - c4)))
  }

  // dD11/dN in case of negative Sigma11 stress
  if (deltaSigma1_c < 0.0 && bigSigma1c > 0.0)
  {
    if (dfc < 1.0)
    {
      tmp = -c2 * dfc / (sqrt(bigSigma1c) * (1.0 + pow(dft, 2.0) + pow(D12, 2.0)));
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddfc = pow( (c1 * (1.0 + pow(dft, 2) + pow(D12, 2.0)) * bigSigma1c * exp(tmp)) , (1.0 + 2.0 * exp( - c6 * dft - D12)) ) + c3 * dfc * pow(bigSigma1c, 2.0) * (1.0 + sqrt(dft) * exp(c8 * sqrt(dft))  / (1.0 + exp( - c5 * (bigSigma1c - c7))) ) * (1.0 + exp(c5 * (bigSigma1c - c4) / 3.0));

      if (max_incr_dam_var / ddfc > NJUMP_threshold) {
        Njump_dfc = NJUMP_threshold;
      } else {
        Njump_dfc = max_incr_dam_var / ddfc;
      }
    } else {
      Njump_dfc = NJUMP_threshold;
    }
    // opserr << "tmp ddfc = " << tmp << endln;

  }

  // dD22/dN in case of positive Sigma22 stress
  if (deltaSigma2_t > 0.0 && bigSigma2t > 0.0)
  {
    if (dmt < 1.0)
    {
      tmp = -c2 * dmt / (sqrt(bigSigma2t) * (1.0 + pow(dmc, 2.0) + pow(D12, 2.0)));
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddmt = c1 * (1.0 + pow(dmc, 2.0) + pow(D12, 2.0)) * bigSigma2t * exp(tmp) + c3 * dmt * pow(bigSigma2t, 2.0) * (1.0 + sqrt(dmc) * exp(c8 * sqrt(dmc)) / (1.0 + exp( - c5 * (bigSigma2t - c7)))) * (1.0 + exp(c5 * (bigSigma2t - c4)));

      if (max_incr_dam_var / ddmt > NJUMP_threshold) {
        Njump_dmt = NJUMP_threshold;
      } else {
        Njump_dmt = max_incr_dam_var / ddmt;
      }
    } else {
      Njump_dmt = NJUMP_threshold;
    }
  }

  // dD22/dN in case of negative Sigma22 stress
  if (deltaSigma2_c < 0.0 && bigSigma2c > 0.0)
  {
    if (dmc < 1.0)
    {
      tmp = -c2 * dmc / (sqrt(bigSigma2c) * (1.0 + pow(dmt, 2.0) + pow(D12, 2.0)));
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddmc = pow(c1 * (1.0 + pow(dmt, 2.0) + pow(D12, 2.0)) * bigSigma2c * exp(tmp) , (1.0 + 2.0 * exp( - c6 * dmt - D12)) ) + c3 * dmc * pow(bigSigma2c, 2.0) * (1.0 + sqrt(dmt) * exp(c8 * sqrt(dmt)) / (1.0 + exp( - c5 * (bigSigma2c - c7)))) * (1.0 + exp(c5 * (bigSigma2c - c4) / 3.0));

      if (max_incr_dam_var / ddmc > NJUMP_threshold) {
        Njump_dmc = NJUMP_threshold;
      } else {
        Njump_dmc = max_incr_dam_var / ddmc;
      }
    } else {
      Njump_dmc = NJUMP_threshold;
    }
  }

  // dD12/dN in case of positive Sigma12 stress
  if (deltaSigma12_t > 0.0 && bigSigma12t > 0.0)
  {
    if (dst < 1.0)
    {
      // tmp = -c2 * dst / (sqrt(bigSigma12t + eps) * 2.0 * (1.0 + pow(dsc, 2))); // don't know what the purpose of eps is
      tmp = -c2 * dst / (sqrt(bigSigma12t + eps) * 2.0 * (1.0 + pow(dsc, 2.0))); // don't know what the purpose of eps is
      // tmp = -c2 * dst / (sqrt(bigSigma12t) * 2.0 * (1.0 + pow(dsc, 2.0))); // don't know what the purpose of eps is
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddst = c1 * (1.0 + pow(dsc, 2.0)) * bigSigma12t * exp(tmp);

      if (max_incr_dam_var / ddst > NJUMP_threshold) {
        Njump_dst = NJUMP_threshold;
      } else {
        Njump_dst = max_incr_dam_var / ddst;
      }
    } else {
      Njump_dst = NJUMP_threshold;
    }
  }

  // dD12/dN in case of negative Sigma12 stress
  if (deltaSigma12_c < 0.0 && bigSigma12c > 0.0)
  {
    if (dsc < 1.0)
    {
      // tmp = -c2 * dsc / (sqrt(bigSigma12c + eps) * 2.0 * (1.0 + pow(dst, 2.0)));
      tmp = -c2 * dsc / (sqrt(bigSigma12c) * 2.0 * (1.0 + pow(dst, 2.0)));
      if (tmp < log(1e-280)) tmp = log(1e-280); // maybe this is a bound check for the exponential or something
      ddsc = c1 * (1.0 + pow(dst, 2.0)) * bigSigma12c * exp(tmp);

      if (max_incr_dam_var / ddsc > NJUMP_threshold) {
        Njump_dsc = NJUMP_threshold;
      } else {
        Njump_dsc = max_incr_dam_var / ddsc;
      }
    } else {
      Njump_dsc = NJUMP_threshold;
    }
    // tmp = -c2 * dsc / (sqrt(bigSigma12c + eps) * 2.0 * (1.0 + pow(dst, 2))); // don't know what the purpose of eps is

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //  Computes the derivatives of the plastic strain !
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  // Getting actual strains
  const Vector& epsilon = getStrain();

  if (deltaSigma1_t > 0.0) {
    dstrain1_p = c9 * epsilon(0) * max(ddst, ddsc);
  } else {
    dstrain1_p = 0.0;
  }

  if (deltaSigma2_t > 0.0) {
    dstrain2_p = c9 * epsilon(1) * max(ddst, ddsc);
  } else {
    dstrain2_p = 0.0;
  }

  // opserr << "Njump_dft = " << Njump_dft << endln;
  // opserr << "ddft = " << ddft << endln;
  // opserr << "Njump_dfc = " << Njump_dfc << endln;
  // opserr << "ddfc = " << ddfc << endln << endln;


  // opserr << "Njump_dmt = " << Njump_dmt << endln;
  // opserr << "ddmt = " << ddmt << endln;
  // opserr << "Njump_dmc = " << Njump_dmc << endln;
  // opserr << "ddmc = " << ddmc << endln << endln;


  // opserr << "Njump_dst = " << Njump_dst << endln;
  // opserr << "ddst = " << ddst << endln;
  // opserr << "Njump_dsc = " << Njump_dsc << endln;
  // opserr << "ddsc = " << ddsc << endln << endln;

  // Calculate NJUMP_local(3) for fiber, membrane and shear. Material can be only in tension or compression
  NJUMP_local.Zero();
  if (Njump_dft == 0) NJUMP_local(0) = Njump_dfc ;
  else NJUMP_local(0) = Njump_dft;

  if (Njump_dmt == 0) NJUMP_local(1) = Njump_dmc ;
  else NJUMP_local(1) = Njump_dmt;

  if (Njump_dst == 0) NJUMP_local(2) = Njump_dsc ;
  else NJUMP_local(2) = Njump_dst;

}

ID VonPapaDamage::getNJUMP(double max_incr_dam_var)
{

  // opserr << "VonPapaDamage::getNJUMP called!" << endln;

  calculateDerDamage(max_incr_dam_var);

  return NJUMP_local;
}

void VonPapaDamage :: advanceDamageState(int Ncycles)
{

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //     Update the fatigue damage variables by using the Forward Euler algorithm        !
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  // opserr << "dft = " << dft << endln;
  // opserr << "dfc = " << dfc << endln;
  // opserr << "dmt = " << dmt << endln;
  // opserr << "dmc = " << dmc << endln;
  // opserr << "dst = " << dst << endln;
  // opserr << "dsc = " << dsc << endln << endln;



  dft += ddft * Ncycles;
  if (dft >= 1.0) dft = 0.999;

  dfc += ddfc * Ncycles;
  if (dfc >= 1.0) dfc = 0.999;

  dmt += ddmt * Ncycles;
  if (dmt >= 1.0) dmt = 0.999;

  dmc += ddmc * Ncycles;
  if (dmc >= 1.0) dmc = 0.999;

  dst += ddst * Ncycles;
  if (dst >= 1.0) dst = 0.999;

  dsc += ddsc * Ncycles;
  if (dsc >= 1.0) dsc = 0.999;

  // dft += 0.05/2;
  // dfc += 0.05/2;
  // dmt += 0.05/2;
  // dmc += 0.05/2;
  // dst += 0.05/2;
  // dsc += 0.05/2;


  // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  // !      Combine the derivatives computed for tensile/compressive stresses     !
  // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  // D11 = fiber = df
  // D22 = matrix = dm
  // D12 = shear = ds

  // Get damaged stresses
  Vector sigma(3);
  sigma = this->getStress();


  // D11
  if (sigma(0) >= 0.0) {
    D11 = dft + dfc;
  } else {
    D11 = dft + b * dfc;
    // D11 = b * dft + dfc; // According to Paepegem 2004 it should be this
  }

  // D22
  if (sigma(1) >= 0.0) {
    D22 = dmt + dmc;
  } else {
    D22 = dmt + b * dmc;
    // D22 = b * dmt + dmc; // According to Paepegem 2004 it should be this
  }

  // D12
  D12 = dst + dsc;



  // Check bounds
  if (D11 >= 1.0) D11 = 0.999;
  if (D22 >= 1.0) D22 = 0.999;
  if (D12 >= 1.0) D12 = 0.999;

  // if (D11 >= 0.999) opserr << "D11 = " << D11 << endln;
  // if (D22 >= 0.999) opserr << "D22 = " << D22 << endln;
  // if (D12 >= 0.999) opserr << "D12 = " << D12 << endln;




  // opserr << "D11 = " << D11 << endln;
  // opserr << "D22 = " << D22 << endln;
  // opserr << "D12 = " << D12 << endln << endln;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //     Update the plastic strains      !
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  strain1_p += dstrain1_p * Ncycles;
  strain2_p += dstrain2_p * Ncycles;
}

#define PARAM_ADVANCESTATE 1
#define PARAM_RESETSTRESSES 2
#define PARAM_CALCULATENJUMP 3

int VonPapaDamage :: setParameter(const char **argv, int argc, Parameter &param)
{
  // opserr << "VonPapaDamage :: setParameter called" << endln;
  // opserr << "argv[0] = " << argv[0] << endln;
  if (argc < 1)
    return -1;

  // int theMaterialTag;
  // theMaterialTag = atoi(argv[1]);
  // if (theMaterialTag == this->getTag()) { // I think this shouldn't be here
  if (strcmp(argv[0], "advanceDamageState") == 0) {    // enforce elastic/elastoplastic response
    // opserr << " VonPapaDamage :: PARAM_ADVANCESTATE  called" << endln;
    return param.addObject(PARAM_ADVANCESTATE, this);
  }
  if (strcmp(argv[0], "calculateNJUMP") == 0) {    // enforce elastic/elastoplastic response
    // opserr << " VonPapaDamage :: PARAM_ADVANCESTATE  called" << endln;
    return param.addObject(PARAM_CALCULATENJUMP, this);
  }
  if (strcmp(argv[0], "resetMaxStress") == 0) {    // enforce elastic/elastoplastic response
    return param.addObject(PARAM_RESETSTRESSES, this);
  }
  // }


  return -1;
}

int VonPapaDamage :: updateParameter(int parameterID, Information &info)
{
  // desde python  entregar Ncycles...
  // llamar advanceDamageState si es lo que se pide....

  // opserr << " VonPapaDamage :: updateParameter called" << endln;

  // called updateMaterialStage in tcl file
  if (parameterID == PARAM_ADVANCESTATE) {
    int Ncycles = (int) info.theDouble;
    // opserr << "Ncycles = " << Ncycles << endln;
    this->advanceDamageState(Ncycles);
  }
  if (parameterID == PARAM_RESETSTRESSES) {
    this->resetMaxStress();
  }
  if (parameterID == PARAM_CALCULATENJUMP) {
    int maxdamage = info.theDouble;
    // this->computeNJUMP(maxdamage);
    this->getNJUMP(maxdamage);
  }

  return -1;
}

double VonPapaDamage::proot_quadraticequ(double a, double b, double c)
{
  double proot = 0;

  if ( abs(a) > 1e-50 && pow(b, 2) - 4.0 * a * c >= 0)
  {
    proot = (-b + sqrt( pow(b, 2) - 4.0 * a * c)) / (2.0 * a); // Positive root
  } else {
    proot = -c / (b + 1e-10); // Positive root
  }

  return proot;
}


// void VonPapaDamage::computeNJUMP(double maxdamage)
// {
//   if(NJUMPVEC == 0)
//   {
//     NJUMPVEC = new int [NVonPapaMaterials]; //usar un ID
//   }

//   // calcular NJUMP
//   NJUMPVEC[i_current_material_point] = NJUMP;
//   i_current_material_point ++;
// }




#define RESPONSE_STRESS 1
#define RESPONSE_STRAIN 2
#define RESPONSE_NJUMP 3
#define RESPONSE_DAMAGESTATE 4


Response*
VonPapaDamage::setResponse (const char **argv, int argc, OPS_Stream &output)
{

  Response * theResponse = 0;

  output.tag("NdMaterialOutput");
  output.attr("matType", this->getClassType());
  output.attr("matTag", this->getTag());

  // opserr << "VonPapaDamage::setResponse called with arv[0] " << argv[0] << endln;
  if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0){
    output.tag("ResponseType","sigma11");
    output.tag("ResponseType","sigma22");
    output.tag("ResponseType","sigma12");
    theResponse = new MaterialResponse(this, RESPONSE_STRESS, this->getStress());
  }
  else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0){
    output.tag("ResponseType","e11");
    output.tag("ResponseType","e22");
    output.tag("ResponseType","e12");
    theResponse = new MaterialResponse(this, RESPONSE_STRAIN, this->getStrain());
  }  
  else if (strcmp(argv[0], "damagestate") == 0 || strcmp(argv[0], "DamageState") == 0){
    output.tag("ResponseType","dft");
    output.tag("ResponseType","dfc");
    output.tag("ResponseType","dmt");
    output.tag("ResponseType","dmc");
    output.tag("ResponseType","dst");
    output.tag("ResponseType","dsc");
    output.tag("ResponseType","D11");
    output.tag("ResponseType","D22");
    output.tag("ResponseType","D12");
    theResponse = new MaterialResponse(this, RESPONSE_DAMAGESTATE, this->getDamageState());
  }
  else if (strcmp(argv[0], "NJUMP") == 0) {
    double maxDamage = strtod(argv[1], nullptr);
    // opserr << "maxDamage = " << maxDamage << endln;
    theResponse = new MaterialResponse(this, RESPONSE_NJUMP, this -> getNJUMP(maxDamage) );
  }
  
  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int
VonPapaDamage::getResponse(int responseID, Information &matInfo)
{

  // opserr << "VonPapaDamage::getResponse called" << endln;

  switch (responseID) {
  case -1:
    return -1;
  case RESPONSE_STRESS:
    if (matInfo.theVector != 0)
      // *(matInfo.theVector) = getStress();
      return matInfo.setVector(this->getStress());
    // return 0;
  case RESPONSE_STRAIN:
    if (matInfo.theVector != 0)
      // *(matInfo.theVector) = getStrain();
      return matInfo.setVector(this->getStrain());
    // return 0;
  case RESPONSE_DAMAGESTATE:
    if (matInfo.theVector != 0)
      // *(matInfo.theVector) = getStrain();
      return matInfo.setVector(this->getDamageState());
    // return 0;
  case RESPONSE_NJUMP:
    if (matInfo.theVector != 0)
      // *(matInfo.theID) = getNJUMP();
      *(matInfo.theID) = NJUMP_local;
      // return matInfo.setID(NJUMP_local);
    // opserr << "NJUMP_local = "<< NJUMP_local << endln;
      return 0;

  default:
    return -1;
  }
}



// int VonPapaDamage :: activateParameter(int parameterID)
// {
//   return -1;
// }

// int VonPapaDamage :: setVariable(const char *variable, Information &)
// {
//   return -1;
// }

// int VonPapaDamage :: getVariable(const char *variable, Information &)
// {
//   return -1;
// }

const Vector& VonPapaDamage :: getDamageState() const
{
    static Vector damage_state(9);
    damage_state(0) = dft;
    damage_state(1) = dfc;
    damage_state(2) = dmt;
    damage_state(3) = dmc;
    damage_state(4) = dst;
    damage_state(5) = dsc;
    damage_state(6) = D11;
    damage_state(7) = D22;
    damage_state(8) = D12;

    return damage_state;
}
