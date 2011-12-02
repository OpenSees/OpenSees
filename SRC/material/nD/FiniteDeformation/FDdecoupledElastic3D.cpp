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
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================

#include <FDdecoupledElastic3D.h>

//-----------------------------------------------------------------------------------------------------------------------------------------------
FDdecoupledElastic3D::FDdecoupledElastic3D(int tag,
                                           int classTag,
                                           WEnergy *wEnergy_in,
                                           double rho_in= 0.0)
:FiniteDeformationElastic3D(tag, classTag, rho_in)
{
    if ( wEnergy_in )
    {
       W = wEnergy_in->newObj();
    }
    else
    {
      opserr << "FDdecoupledElastic3D:: FDdecoupledElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

FDdecoupledElastic3D::FDdecoupledElastic3D(int tag,
                                           WEnergy *wEnergy_in,
                                           double rho_in = 0.0)
:FiniteDeformationElastic3D(tag, ND_TAG_FDdecoupledElastic3D, rho_in)
{
    if ( wEnergy_in)
    {
       W = wEnergy_in->newObj();
    }
    else
    {
      opserr << "FDdecoupledElastic3D:: FDdecoupledElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

FDdecoupledElastic3D::FDdecoupledElastic3D(int tag,
                                           WEnergy *wEnergy_in)
:FiniteDeformationElastic3D(tag, ND_TAG_FDdecoupledElastic3D, 0.0)
{
    if ( wEnergy_in )
    {
       W = wEnergy_in->newObj();
    }
    else
    {
     opserr << "FDdecoupledElastic3D:: FDdecoupledElastic3D failed to construct the W Energy\n";
      exit(-1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------------
FDdecoupledElastic3D::FDdecoupledElastic3D( )
:FiniteDeformationElastic3D(0, 0, 0.0)
{
    W = 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------------
FDdecoupledElastic3D::FDdecoupledElastic3D(FDdecoupledElastic3D &fde3d )
{    
    F = fde3d.F;
    C = fde3d.C;
    Cinv = fde3d.Cinv;
    J = fde3d.J;
    lambda1 = fde3d.lambda1;
    lambda2 = fde3d.lambda2;
    lambda3 = fde3d.lambda3;
    lambda_wave1 = fde3d.lambda_wave1;
    lambda_wave2 = fde3d.lambda_wave2;
    lambda_wave3 = fde3d.lambda_wave3;

    Stiffness = fde3d.Stiffness;
    thisGreenStrain = fde3d.thisGreenStrain;
    thisPK2Stress = fde3d.thisPK2Stress;    
}

//------------------------------------------------------------------------------------------------------------------------------------------------
FDdecoupledElastic3D::~FDdecoupledElastic3D()
{
   if (W)
     delete W;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------
double FDdecoupledElastic3D::getRho(void)
{
   return rho;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
WEnergy *FDdecoupledElastic3D::getWEnergy(void)
{
   return W;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::setTrialF(const straintensor &f)
{
   FromForC = 0;
   F = f;
   C = F("ki")*F("kj");   C.null_indices();
   return this->ComputeTrials();
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::setTrialFIncr(const straintensor &df)
{
   return this->setTrialF(this->getF() + df);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::setTrialC(const straintensor &c)
{
   FromForC = 1;
   C = c;
   return this->ComputeTrials();
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::setTrialCIncr(const straintensor &dc)
{
   return this->setTrialC(this->getC() + dc);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor FDdecoupledElastic3D::getF(void)
{
   return F;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor FDdecoupledElastic3D::getC(void)
{
   return C;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
const double FDdecoupledElastic3D::getJ(void)
{
   return J;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FDdecoupledElastic3D::getlambda(void)
{
  Vector lambda(3);

  lambda(0) = lambda1;
  lambda(1) = lambda2;
  lambda(2) = lambda3;

  return lambda;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FDdecoupledElastic3D::getlambda_wave(void)
{
  Vector lambda_wave(3);
  lambda_wave(0) = lambda_wave1;
  lambda_wave(1) = lambda_wave2;
  lambda_wave(2) = lambda_wave3;

  return lambda_wave;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
const Vector FDdecoupledElastic3D::wa(void)
{
  Vector Wa(3);
  Vector lambda_wave(3);
  lambda_wave = this->getlambda_wave();
  Vector disowOverlambda = W->disowOdlambda(lambda_wave);
  double temp = disowOverlambda(0) * lambda_wave(0) +
                disowOverlambda(1) * lambda_wave(1) +
                disowOverlambda(2) * lambda_wave(2) ;
  temp = temp * (-0.3333333333333333333333333333);
  Wa(0) = temp + disowOverlambda(0) * lambda_wave(0);
  Wa(1) = temp + disowOverlambda(1) * lambda_wave(1);
  Wa(2) = temp + disowOverlambda(2) * lambda_wave(2);
  return Wa;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FDdecoupledElastic3D::Yab(void)
{
        Tensor Y(2, def_dim_2, 0.0);
        Tensor I_ij("I", 2, def_dim_2);
        Vector lambda_wave(3);
        lambda_wave = this->getlambda_wave();
        Tensor  d2 = W->d2isowOdlambda1dlambda2(lambda_wave);
        Vector  d1 = W->disowOdlambda(lambda_wave);
        Vector  d11 = W->d2isowOdlambda2(lambda_wave);
	d2.val(1,1) = d11(0);
        d2.val(2,2) = d11(1);
        d2.val(3,3) = d11(2);
        Vector tempi(3);
        double tempd = d1(0)*lambda_wave(0) + d1(1)*lambda_wave(1) + d1(2)*lambda_wave(2) ;
        double tempcd = 0.0;
        for (int i=0; i<3; i++)
        {
          tempi(i) = 0.0;
          for (int j=0; j<3; j++)
          {
              tempi(i) += d2.cval(i+1,j+1) * lambda_wave(i) * lambda_wave(j);
              tempcd   += d2.cval(i+1,j+1) * lambda_wave(i) * lambda_wave(j);
          }
        }
        for(int a=1; a<=3; a++)
        {
          for(int b=1; b<=3; b++)
          {
              Y.val(a,b) = d1(a-1)*I_ij.cval(a,b)*lambda_wave(b-1) + d2.cval(a,b)*lambda_wave(a-1)*lambda_wave(b-1) -
                           (  tempi(a-1) + tempi(b-1) + d1(a-1)*lambda_wave(a-1) + d1(b-1)*lambda_wave(b-1) ) / 3.0 +
                           ( tempcd + tempd ) / 9.0;
          }
        }
        return Y;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FDdecoupledElastic3D::FDisoStiffness(void)
{
  Tensor I_ij("I", 2, def_dim_2);
  Tensor I_ijkl( 4, def_dim_4, 0.0 );
  I_ijkl = I_ij("ij") * I_ij("kl");
  I_ij.null_indices();
  Tensor I_ikjl( 4, def_dim_4, 0.0 );
  I_ikjl = I_ijkl.transpose0110();
  Tensor I_iljk( 4, def_dim_4, 0.0 );
  I_iljk = I_ijkl.transpose0111();
  Tensor I4s = (I_ikjl+I_iljk)*0.5;
  Tensor  tempI = I4s - I_ijkl;

  Tensor CinvCinv = Cinv("ij") * Cinv("kl");
  CinvCinv.null_indices(); Cinv.null_indices(); 

  Tensor ICinv = ( CinvCinv.transpose0110() + CinvCinv.transpose0111() ) * (0.5);

  Tensor CinvCinv_ICinv = CinvCinv - ICinv;

  double I1 = lambda1*lambda1 + lambda2*lambda2 + lambda3*lambda3;

  Vector Wa = this->wa();
  Tensor yab = this->Yab();

  Tensor L_iso(2,def_dim_2,0.0);

  if(caseIndex == 0)
  {
    double d1 = (lambda1+lambda2)*(lambda1+lambda3)*(lambda1-lambda2)*(lambda1-lambda3);
    double d2 = (lambda2+lambda3)*(lambda2+lambda1)*(lambda2-lambda3)*(lambda2-lambda1);
    double d3 = (lambda3+lambda1)*(lambda3+lambda2)*(lambda3-lambda1)*(lambda3-lambda2);

    Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv*(J*J/(lambda1*lambda1)) ) * (1.0/d1);
    Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv*(J*J/(lambda2*lambda2)) ) * (1.0/d2);
    Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv*(J*J/(lambda3*lambda3)) ) * (1.0/d3);

    double d1p = 4.0 *lambda1*lambda1*lambda1*lambda1 - I1*lambda1*lambda1 - J*J /(lambda1*lambda1);
    double d2p = 4.0 *lambda2*lambda2*lambda2*lambda2 - I1*lambda2*lambda2 - J*J /(lambda2*lambda2);
    double d3p = 4.0 *lambda3*lambda3*lambda3*lambda3 - I1*lambda3*lambda3 - J*J /(lambda3*lambda3);

    Tensor Cm1M1M1Cm1 = Cinv("ij")*M1("kl") + M1("ij")*Cinv("kl");
    Cinv.null_indices(); M1.null_indices(); Cm1M1M1Cm1.null_indices();

    Tensor Cm1M2M2Cm1 = Cinv("ij")*M2("kl") + M2("ij")*Cinv("kl");
    Cinv.null_indices(); M2.null_indices(); Cm1M2M2Cm1.null_indices();

    Tensor Cm1M3M3Cm1 = Cinv("ij")*M3("kl") + M3("ij")*Cinv("kl");
    Cinv.null_indices(); M3.null_indices(); Cm1M3M3Cm1.null_indices();


    Tensor dM1M1d = I_ij("ij")*M1("kl") + M1("ij")*I_ij("kl");
    I_ij.null_indices(); M1.null_indices(); dM1M1d.null_indices();

    Tensor dM2M2d = I_ij("ij")*M2("kl") + M2("ij")*I_ij("kl");
    I_ij.null_indices(); M2.null_indices(); dM2M2d.null_indices();

    Tensor dM3M3d = I_ij("ij")*M3("kl") + M3("ij")*I_ij("kl");
    I_ij.null_indices(); M3.null_indices(); dM3M3d.null_indices();

    Tensor M1M1 = M1("ij") * M1("kl");
    M1.null_indices(); M1M1.null_indices();
    Tensor M2M2 = M2("ij") * M2("kl");
    M2.null_indices(); M2M2.null_indices();
    Tensor M3M3 = M3("ij") * M3("kl");
    M3.null_indices(); M3M3.null_indices();

    Tensor calM1 = ( tempI + (CinvCinv_ICinv -Cm1M1M1Cm1)*(J*J/(lambda1*lambda1)) + dM1M1d*(lambda1*lambda1) - M1M1*d1p ) *(1.0/d1);
    Tensor calM2 = ( tempI + (CinvCinv_ICinv -Cm1M2M2Cm1)*(J*J/(lambda2*lambda2)) + dM2M2d*(lambda2*lambda2) - M2M2*d2p ) *(1.0/d2);
    Tensor calM3 = ( tempI + (CinvCinv_ICinv -Cm1M3M3Cm1)*(J*J/(lambda3*lambda3)) + dM3M3d*(lambda3*lambda3) - M3M3*d3p ) *(1.0/d3);

    Tensor L_iso_1 = ( calM1*Wa(0) + calM2*Wa(1) + calM3*Wa(2) ) * 2.0;
    Tensor L_iso_2 =  M1("ij") * M1("kl") * yab.cval(1,1)  + M1("ij") * M2("kl") * yab.cval(1,2)  + M1("ij") * M3("kl") * yab.cval(1,3)  +
                      M2("ij") * M1("kl") * yab.cval(2,1)  + M2("ij") * M2("kl") * yab.cval(2,2)  + M2("ij") * M3("kl") * yab.cval(2,3)  +
                      M3("ij") * M1("kl") * yab.cval(3,1)  + M3("ij") * M2("kl") * yab.cval(3,2)  + M3("ij") * M3("kl") * yab.cval(3,3);
    L_iso = L_iso_1 + L_iso_2 ;
  }

  if(caseIndex == 11)
  {
    double d1 = (lambda1+lambda2)*(lambda1+lambda3)*(lambda1-lambda2)*(lambda1-lambda3);
    Tensor M1 = (I_ij - Cinv * (lambda2*lambda2)) * (1.0/(lambda1+lambda2)/(lambda1-lambda2));
    Tensor Mr = Cinv - M1;
    double d1p = 4.0 *lambda1*lambda1*lambda1*lambda1 - I1*lambda1*lambda1 - J*J /(lambda1*lambda1);
    Tensor Cm1M1M1Cm1 = Cinv("ij")*M1("kl") + M1("ij")*Cinv("kl");
    Cinv.null_indices(); M1.null_indices(); Cm1M1M1Cm1.null_indices();
    Tensor dM1M1d = I_ij("ij")*M1("kl") + M1("ij")*I_ij("kl");
    I_ij.null_indices(); M1.null_indices(); dM1M1d.null_indices();
    Tensor M1M1 = M1("ij") * M1("kl");
    M1.null_indices(); M1M1.null_indices();
    Tensor calM1 = ( tempI + (CinvCinv_ICinv -Cm1M1M1Cm1)*(J*J/(lambda1*lambda1)) + dM1M1d*(lambda1*lambda1) - M1M1*d1p ) *(1.0/d1);
    Tensor calMr = (ICinv + calM1) * (-1.0);
    Tensor L_iso_1 = ( calM1*Wa(0) + calMr*Wa(2) ) * 2.0;
    Tensor L_iso_2 =  M1("ij") * M1("kl") * yab.cval(1,1)  + M1("ij") * Mr("kl") * yab.cval(1,3)  +
                      Mr("ij") * M1("kl") * yab.cval(3,1)  + Mr("ij") * Mr("kl") * yab.cval(3,3);
    L_iso = L_iso_1 + L_iso_2 ;
  }

  if(caseIndex == 13)
  {
    double d3 = (lambda3+lambda1)*(lambda3+lambda2)*(lambda3-lambda1)*(lambda3-lambda2);
    Tensor M3 = (I_ij - Cinv * (lambda2*lambda2)) * (1.0/(lambda3+lambda2)/(lambda3-lambda2));
    Tensor Mr = Cinv - M3;
    double d3p = 4.0 *lambda3*lambda3*lambda3*lambda3 - I1*lambda3*lambda3 - J*J /(lambda3*lambda3);
    Tensor Cm1M3M3Cm1 = Cinv("ij")*M3("kl") + M3("ij")*Cinv("kl");
    Cinv.null_indices(); M3.null_indices(); Cm1M3M3Cm1.null_indices();
    Tensor dM3M3d = I_ij("ij")*M3("kl") + M3("ij")*I_ij("kl");
    I_ij.null_indices(); M3.null_indices(); dM3M3d.null_indices();
    Tensor M3M3 = M3("ij") * M3("kl");
    M3.null_indices(); M3M3.null_indices();
    Tensor calM3 = ( tempI + (CinvCinv_ICinv -Cm1M3M3Cm1)*(J*J/(lambda3*lambda3)) + dM3M3d*(lambda3*lambda3) - M3M3*d3p ) *(1.0/d3);
    Tensor calMr = (ICinv + calM3) * (-1.0);
    Tensor L_iso_1 = ( calM3*Wa(2) + calMr*Wa(0) ) * 2.0;
    Tensor L_iso_2 =  M3("ij") * M3("kl") * yab.cval(3,3)  + M3("ij") * Mr("kl") * yab.cval(3,1)  +
                      Mr("ij") * M3("kl") * yab.cval(1,3)  + Mr("ij") * Mr("kl") * yab.cval(1,1);
    L_iso = L_iso_1 + L_iso_2 ;
  }

    if(caseIndex == 2)
  {
    Vector lambda_wave(3);
    lambda_wave = this->getlambda_wave();
    Vector  d11 = W->d2isowOdlambda2(lambda_wave);
    Vector  d1 = W->disowOdlambda(lambda_wave);
    double G2linear = d11(1)*lambda_wave2*lambda_wave2 + d1(1)*lambda_wave2;

    L_iso = ( ICinv - CinvCinv * (1.0/3.0) ) * G2linear;
  }

    return L_iso;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor FDdecoupledElastic3D::FDvolStiffness(void)
{
   Tensor CinvCinv = Cinv("ij")*Cinv("kl") ;
   Cinv.null_indices(); CinvCinv.null_indices();
   Tensor ICinv = ( CinvCinv.transpose0110() + CinvCinv.transpose0111() ) * (0.5);
   double dWdJ = W->dvolwOdJ(J);
   double d2WdJ2 = W->d2volwOdJ2(J);
   double wj = d2WdJ2*J*J + J*dWdJ;

   Tensor L_vol = CinvCinv*wj - ICinv *2.0*J*dWdJ  ;

   return L_vol;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor& FDdecoupledElastic3D::getTangentTensor(void)
{
    return Stiffness;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor
&FDdecoupledElastic3D::getInitialTangentTensor(void)
{
    //tensor I2("I", 2, def_dim_2);
    //tensor I_ijkl = I2("ij")*I2("kl");
    //tensor I_ikjl = I_ijkl.transpose0110();
    //tensor I_iljk = I_ijkl.transpose0111();
    //tensor I4s = (I_ikjl+I_iljk)*0.5;
    //static tensor L0;
    //L0 = I_ijkl*( E*nu / ( (1.0+nu)*(1.0 - 2.0*nu) ) ) + I4s*( E / (1.0 + nu) );

    //return L0;
    return this->getTangentTensor();
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor FDdecoupledElastic3D::getStrainTensor(void)
{
   return thisGreenStrain;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FDdecoupledElastic3D::getStressTensor(void)
{
   return thisPK2Stress;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FDdecoupledElastic3D::getPK1StressTensor(void)
{
   stresstensor thisSPKStress;
   stresstensor thisFPKStress;

   if ( FromForC == 0 ) {
    thisSPKStress = this->getStressTensor();
    thisFPKStress = thisSPKStress("ij") * (F.transpose11())("jk") ;
   }

   if ( FromForC == 1 ) {
    opserr << "FDdecoupledElastic3D: unknown deformation gradient - cannot compute PK1 stress" << "\n";
    exit (-1);
   }

    return thisFPKStress;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor FDdecoupledElastic3D::getCauchyStressTensor(void)
{
   stresstensor thisSPKStress;
   stresstensor thisCauchyStress;

   if ( FromForC == 0 ) {
    thisSPKStress = this->getStressTensor();
    thisCauchyStress = F("ij") * thisSPKStress("jk") * (F.transpose11())("kl") * (1.0/J);
   }

   if ( FromForC == 1 ) {
    opserr << "FDdecoupledElastic3D: unknown deformation gradient - cannot compute Cauchy stress" << "\n";
    exit (-1);
   }

    return thisCauchyStress;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::commitState (void)
{
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::revertToLastCommit (void)
{
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::revertToStart (void)
{
   Tensor F0("I", 2, def_dim_2);
   F = F0;
   C = F0;
   Cinv = F0;

   Tensor ss_zero(2,def_dim_2,0.0);
   thisPK2Stress = ss_zero;
   thisGreenStrain = ss_zero;
   Stiffness = getInitialTangentTensor();

   J = 1.0;
   lambda1 = 1.0;
   lambda2 = 1.0;
   lambda3 = 1.0;
   lambda_wave1 = 1.0;
   lambda_wave2 = 1.0;
   lambda_wave3 = 1.0;

   caseIndex = 0;

   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * FDdecoupledElastic3D::getCopy (void)
{
    FDdecoupledElastic3D   *theCopy =
    new FDdecoupledElastic3D (this->getTag(), this->getWEnergy(), this->getRho());

    theCopy->F = F;
    theCopy->C = C;
    theCopy->Cinv = Cinv;
    theCopy->J = J;
    theCopy->lambda1 = lambda1;
    theCopy->lambda2 = lambda2;
    theCopy->lambda3 = lambda3;
    theCopy->lambda_wave1 = lambda_wave1;
    theCopy->lambda_wave2 = lambda_wave2;
    theCopy->lambda_wave3 = lambda_wave3;

    theCopy->Stiffness = Stiffness;
    theCopy->thisGreenStrain = thisGreenStrain;
    theCopy->thisPK2Stress = thisPK2Stress;

    return theCopy;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * FDdecoupledElastic3D::getCopy (const char *type)
{

  opserr << "FDdecoupledElastic3D::getCopy(const char *) - not yet implemented\n";

    return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const char* FDdecoupledElastic3D::getType (void) const
{
   return "ThreeDimentionalFD";
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::getOrder (void) const
{
   return 6;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::sendSelf (int commitTag, Channel &theChannel)
{
   int res = 0;

   static Vector data(2);

   data(0) = this->getTag();
   data(1) = rho;

   res += theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0)
    {
      opserr << "FDdecoupledElastic3D::sendSelf -- could not send Vector\n";
      return res;
    }

   return res;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::recvSelf (int commitTag,
                                          Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
   int res = 0;

   static Vector data(2);

   res += theChannel.recvVector(this->getDbTag(), commitTag, data);
   if (res < 0)
    {
      opserr << "FDdecoupledElastic3D::recvSelf -- could not recv Vector\n";
      return res;
    }

   this->setTag((int)data(0));
   rho = data(1);

   return res;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void FDdecoupledElastic3D::Print (OPS_Stream &s, int flag)
{
   s << "Finite Deformation Elastic 3D model" << endln;
   s << "\trho: " << rho << endln;
   return;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//int FDdecoupledElastic3D::setParameter(char **argv, int argc, Information &info)
//{
//   return -1;
//}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//int FDdecoupledElastic3D::updateParameter(int parameterID, Information &info)
//{
//   return -1;
//}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::ComputeTrials()
{
   // Cinv:
   Cinv = C.inverse();
   Cinv.symmetrize11();

   // J:
   J = sqrt(C.determinant());

   // lambda:
   tensor eigtensor = C.eigenvalues();
   lambda1 = sqrt(eigtensor.cval(1));
   lambda2 = sqrt(eigtensor.cval(2));
   lambda3 = sqrt(eigtensor.cval(3));

   // lambda_wave
   double JJJ = pow(J, -0.33333333333333333333333333333);
   lambda_wave1 = lambda1 *JJJ;
   lambda_wave2 = lambda2 *JJJ;
   lambda_wave3 = lambda3 *JJJ;

   // caseIndex, note lambda1 >= lambda2 >= lambda3 implied by C.eigenvalues()
   double diff12 = fabs(lambda1-lambda2);
   double diff23 = fabs(lambda2-lambda3);
   double perturbation = pow( d_macheps(), (0.4) );
   if ( diff12 >= perturbation && diff23 >= perturbation )
	caseIndex = 0;
   else if (diff12 >= perturbation && diff23 < perturbation )
	caseIndex = 11;
   else if (diff12 < perturbation && diff23 >= perturbation )
	caseIndex = 13;
   else if (diff12 < perturbation &&  diff23 < perturbation )
	caseIndex = 2;
   else   {opserr << "FDdecoupledElastic3D::getCaseIndex -- unknown case! \n";
	exit (-1);}

   Tensor I_ij("I", 2, def_dim_2);

   Tensor isoPK2Stress(2, def_dim_2, 0.0);

   Vector Wa = this->wa();

   double I1 = lambda1*lambda1+lambda2*lambda2+lambda3*lambda3;

   if (caseIndex == 0)
   {
     double d1 = (lambda1+lambda2)*(lambda1+lambda3)*(lambda1-lambda2)*(lambda1-lambda3);
     double d2 = (lambda2+lambda3)*(lambda2+lambda1)*(lambda2-lambda3)*(lambda2-lambda1);
     double d3 = (lambda3+lambda1)*(lambda3+lambda2)*(lambda3-lambda1)*(lambda3-lambda2);

     Tensor M1 = ( C - I_ij*(I1-lambda1*lambda1) + Cinv *(J*J/(lambda1*lambda1)) ) * (1.0/d1);
     Tensor M2 = ( C - I_ij*(I1-lambda2*lambda2) + Cinv *(J*J/(lambda2*lambda2)) ) * (1.0/d2);
     Tensor M3 = ( C - I_ij*(I1-lambda3*lambda3) + Cinv *(J*J/(lambda3*lambda3)) ) * (1.0/d3);

     isoPK2Stress = M1*Wa(0) + M2*Wa(1) + M3*Wa(2);
   }

   if (caseIndex == 11)
   {
     Tensor M1 = (I_ij - Cinv * (lambda2*lambda2)) * (1.0/(lambda1+lambda2)/(lambda1-lambda2));
     Tensor Mr = Cinv - M1;
     isoPK2Stress = Mr*Wa(2) + M1*Wa(0);
   }

   if (caseIndex == 13)
   {
     Tensor M3 = (I_ij - Cinv * (lambda2*lambda2)) * (1.0/(lambda3+lambda2)/(lambda3-lambda2));
     Tensor Mr = Cinv - M3;
     isoPK2Stress = Mr*Wa(0) + M3*Wa(2);
   }

   if (caseIndex == 2)
   {

   }

   double dWdJ = W->dvolwOdJ(J);
   Tensor volPK2Stress = Cinv * J * dWdJ;

   thisPK2Stress = volPK2Stress + isoPK2Stress; // This is PK2Stress

   thisGreenStrain = (C - I_ij) * 0.5; // This is Green Strain

   Tensor L_iso = this->FDisoStiffness();
   Tensor L_vol = this->FDvolStiffness();
   Stiffness = L_iso + L_vol; // This is Langrangian Tangent Stiffness

   return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------
int FDdecoupledElastic3D::getCaseIndex()
{
   return caseIndex;
}


