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

// the traditional neo-Hookean hyperelasticity:
// w = 0.5*lambda*(lnJ)^2 - G*(lnJ) + 0.5*G*(trace(C)-3)

#include <NeoHookeanCompressible3D.h>

//-----------------------------------------------------------------------------------------------------------------------------------------------
NeoHookeanCompressible3D::NeoHookeanCompressible3D(int tag,
                                                   int classTag,
                                                   double K_in,
                                                   double G_in,
                                                   double rho_in = 0.0)
:FiniteDeformationElastic3D(tag, classTag, rho_in), K(K_in), G(G_in)
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------------
NeoHookeanCompressible3D::NeoHookeanCompressible3D(int tag,
                                               double K_in,
                                               double G_in,
                                                   double rho_in = 0.0)
:FiniteDeformationElastic3D(tag, ND_TAG_NeoHookeanCompressible3D, rho_in), K(K_in), G(G_in)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------------
NeoHookeanCompressible3D::NeoHookeanCompressible3D( )
:FiniteDeformationElastic3D(0, 0, 0.0), K(0.0), G(0.0)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------------
NeoHookeanCompressible3D::~NeoHookeanCompressible3D()
{

}

//-------------------------------------------------------------------------------------------------------------------------------------------------
double NeoHookeanCompressible3D::getRho(void)
{
   return rho;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::setTrialF(const straintensor &f)
{
   FromForC = 0;
   F = f;
   C = F("ki")*F("kj");   C.null_indices();
   return this->ComputeTrials();
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::setTrialFIncr(const straintensor &df)
{
   return this->setTrialF(this->getF() + df);
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::setTrialC(const straintensor &c)
{
   FromForC = 1;
   C = c;
   return this->ComputeTrials();
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::setTrialCIncr(const straintensor &dc)
{
   return this->setTrialC(this->getC() + dc);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor NeoHookeanCompressible3D::getF(void)
{
   return F;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor NeoHookeanCompressible3D::getC(void)
{
   return C;
}

////------------------------------------------------------------------------------------------------------------------------------------------------------
//const double NeoHookeanCompressible3D::getJ(void)
//{
//   return J;
//}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor& NeoHookeanCompressible3D::getTangentTensor(void)
{
    return Stiffness;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor
&NeoHookeanCompressible3D::getInitialTangentTensor(void)
{
    return this->getTangentTensor();
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor NeoHookeanCompressible3D::getStrainTensor(void)
{
   return thisGreenStrain;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor NeoHookeanCompressible3D::getStressTensor(void)
{
   return thisPK2Stress;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor NeoHookeanCompressible3D::getPK1StressTensor(void)
{
   stresstensor thisSPKStress;
   stresstensor thisFPKStress;

   if ( FromForC == 0 ) {
    thisSPKStress = this->getStressTensor();
    thisFPKStress = thisSPKStress("ij") * (F.transpose11())("jk") ;
   }

   if ( FromForC == 1 ) {
    opserr << "NeoHookeanCompressible3D: unknown deformation gradient - cannot compute PK1 stress" << "\n";
    exit (-1);
   }

    return thisFPKStress;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor NeoHookeanCompressible3D::getCauchyStressTensor(void)
{
   stresstensor thisSPKStress;
   stresstensor thisCauchyStress;

   if ( FromForC == 0 ) {
    thisSPKStress = this->getStressTensor();
    thisCauchyStress = F("ij") * thisSPKStress("jk") * (F.transpose11())("kl") * (1.0/J);
   }

   if ( FromForC == 1 ) {
    opserr << "NeoHookeanCompressible3D: unknown deformation gradient - cannot compute Cauchy stress" << "\n";
    exit (-1);
   }

    return thisCauchyStress;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::commitState (void)
{
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::revertToLastCommit (void)
{
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::revertToStart (void)
{
   tensor F0("I", 2, def_dim_2);
   F = F0;
   C = F0;
   Cinv = F0;
   J = 1.0;

   tensor ss_zero(2,def_dim_2,0.0);
   thisPK2Stress = ss_zero;
   thisGreenStrain = ss_zero;
   
   Stiffness = getInitialTangentTensor();

   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * NeoHookeanCompressible3D::getCopy (void)
{
    NeoHookeanCompressible3D   *theCopy =
    new NeoHookeanCompressible3D (this->getTag(), K, G, rho);

    theCopy->F = F;
    theCopy->C = C;
    theCopy->Cinv = Cinv;
    theCopy->J = J;

    theCopy->Stiffness = Stiffness;
    theCopy->thisGreenStrain = thisGreenStrain;
    theCopy->thisPK2Stress = thisPK2Stress;

    return theCopy;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NDMaterial * NeoHookeanCompressible3D::getCopy (const char *type)
{
   opserr << "NeoHookeanCompressible3D::getCopy(const char *) - not yet implemented\n";
   return 0;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const char* NeoHookeanCompressible3D::getType (void) const
{
   return "ThreeDimentionalFD";
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::getOrder (void) const
{
   return 6;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::sendSelf (int commitTag, Channel &theChannel)
{
   int res = 0;
   // not yet implemented
   return res;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::recvSelf (int commitTag,
                                          Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
   int res = 0;
   // not yet implemented
   return res;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void NeoHookeanCompressible3D::Print (OPS_Stream &s, int flag)
{
   s << "Finite Deformation Elastic 3D model" << "\n";
   s << "\trho: " << rho << "\n";
   s << "\tK: " << K << "\n";
   s << "\tG: " << G << "\n";
   return;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//int NeoHookeanCompressible3D::setParameter(char **argv, int argc, Information &info)
//{
//   return -1;
//}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//int NeoHookeanCompressible3D::updateParameter(int parameterID, Information &info)
//{
//   return -1;
//}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int NeoHookeanCompressible3D::ComputeTrials()
{   
   tensor tensorI2("I", 2, def_dim_2);
   tensor tsr1;
   tensor tsr2;

   // Cinv:
   Cinv = C.inverse();
   Cinv.symmetrize11();

   // J:
   J = sqrt(C.determinant());

   // lame constants:
   double lambda = K - 2.0*G/3.0;
   double mu = G - lambda*log(J);

   // Pk2Stress:
   thisPK2Stress = (tensorI2-Cinv)*G + Cinv*lambda*log(J);
   
   // Green Strain:
   thisGreenStrain = (C - tensorI2) * 0.5; 
   
   // Langrangian Tangent Stiffness:
   tsr1 = Cinv("ij")*Cinv("kl");
     tsr1.null_indices();
   tsr2 = tsr1.transpose0110() + tsr1.transpose0111();
   Stiffness = tsr1*lambda + tsr2*mu;

   return 0;
}

