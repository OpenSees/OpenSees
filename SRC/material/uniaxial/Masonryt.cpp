// 
// Written by: Gonzalo Torrisi, UNCuyo, Argentina, 2015  
// 
// Description: This file contains the class implementation for Masonryt Material, a material model based 
// on Crisafulli's Stress-Strain Law originally developed to model cyclic axial behavior of Masonryt.Masonryt 
// constitutive law is defined in terms of Axial Force - Axial Displacement, considering varying Area for the  
// diagonal strut according to the rate of the Axial deformation. It is thus suitable for elements where Force-
// Deformation relation is considered (e.g. zerolength element) 

#include <stdlib.h> 
#include "Masonryt.h" 
#include <OPS_Globals.h> 
#include <float.h> 
#include <Vector.h> 
#include <Channel.h> 
#include <Information.h> 
#include <algorithm> 
#include <math.h> 
#include "classTags.h" 
#include "elementAPI.h" 
 
/******************************************************************************* 
* GETTING STARTED : Create a new Uniaxial Material class named "Masonryt"  
*******************************************************************************/ 
void * OPS_ADD_RUNTIME_VPV(OPS_Masonryt)
{ 

  // Pointer to a uniaxial material that will be returned 
  //UniaxialMaterial *theMaterial = 0; 
  // 
  // parse the input line for the material parameters 
  // 
  int    iData[1]; 
  double dData[21]; 
  int numData; 
  numData = 1; 
  if (OPS_GetIntInput(&numData, iData) != 0) { 
    opserr << "WARNING invalid uniaxialMaterial Masonryt tag" << endln; 
    return 0; 
  } 
  numData = 21; 
  if (OPS_GetDoubleInput(&numData, dData) != 0) { 
    opserr << "WARNING invalid Masonryt Material Parameters\n"; 
    return 0;   
  } 
  //  
  // create a new material 
  // 
  UniaxialMaterial* mat = new Masonryt(iData[0], dData[0], dData[1],dData[2], dData[3], dData[4], dData[5],dData[6], 
    dData[7],  dData[8],  dData[9],  dData[10], dData[11], dData[12], dData[13], dData[14], dData[15],dData[16], 
    dData[17], dData[18], dData[19], int(dData[20]));        
 
  if (mat == 0) { 
    opserr << "WARNING could not create uniaxialMaterial of type Masonryt\n"; 
    return 0; 
  } 
  // return the material 
  return mat; 
} 
/******************************************************************************* 
 * PART I:  Set material Constructors & Deconstructor & Initialize variables                                           
 *******************************************************************************/ 
Masonryt::Masonryt(int tag, double _Fm,     double _Ft,  
                     double _Um,    double _Uult,   double _Ucl,  
                   double _Emo,   double _Length, double _Area1,  
                   double _Area2, double _D1,     double _D2, 
                   double _Ach,   double _Are,    double _Ba,  
                   double _Bch,   double _Gun,    double _Gplu,   
                   double _Gplr,  double _Exp1,   double _Exp2,   
                   int _IENV): 
  UniaxialMaterial(tag, MAT_TAG_Masonryt), Fm(_Fm), Ft(_Ft),  
                   Um(_Um), Uult(_Uult), Ucl(_Ucl),  
                   Emo(_Emo), Length(_Length), Area1(_Area1), 
                   Area2(_Area2), D1(_D1), D2(_D2), 
                   Ach(_Ach), Are(_Are), Ba(_Ba),  
                   Bch(_Bch), Gun(_Gun), Gplu(_Gplu),  
                   Gplr(_Gplr), Exp1(_Exp1), Exp2(_Exp2),  
                   IENV(_IENV) 
{ 
  // Initialize variables 
    double  A1 = Emo*Um/Fm; 
    double  A2 = 1 - A1*Um/Uult; 
    //double  Ko = Area1*Emo/Length; 
	double Ko = Area1*Emo;
    double  DeltaU, Kfactor; 
    DeltaU  = 0.0; 
    Kfactor = 1; 
//  opserr << " Area1 "<<Area1<<" Length "<<Length<<" ko "<<Ko<<" Emo "<<Emo<<endln;
  // Initialize history variables 
    K      = Ko; 
    Et  = Emo; 
    Uun   = 0.0; 
    Sun  = 0.0; 
    Eun  = Gun*Emo; 
    Ure  = 0.0; 
    Sre  = 0.0; 
    Ere  = Eun; 
    Uch  = 0.0; 
    Sch  = 0.0; 
    Ech  = 0.0; 
    U1  = 0.0; 
    S1  = 0.0; 
    E1  = 0.0; 
    U2  = 0.0; 
    S2  = 0.0; 
    E2  = 0.0; 
    UunInt = 0.0; 
    UreInt = 0.0; 
    Upl    = 0.0; 
    FtRed  = Ft; 
    RuleNo = 3; 
    IVIR   = 1; 
    InnerCycleNo = 1; 
    Area   = Area1; 
	double Uma = Ft / Emo;
 
    cK      = Ko; 
    cEt     = Emo; 
    cUun   = 0.0; 
    cSun    = 0.0; 
    cEun    = Gun*Emo; 
    cUre    = 0.0; 
    cSre    = 0.0; 
    cEre    = Eun; 
    cUch    = 0.0; 
    cSch    = 0.0; 
    cEch    = 0.0; 
    cU1     = 0.0; 
    cS1     = 0.0; 
    cE1     = 0.0; 
    cU2     = 0.0; 
    cS2     = 0.0; 
    cE2     = 0.0; 
    cUunInt = 0.0; 
    cUreInt = 0.0; 
    cUpl    = 0.0; 
    cFtRed  = Ft; 
    cRuleNo = 3; 
    cIVIR   = 1; 
    cInnerCycleNo = 1;  
    cArea   = Area1; 


  this->revertToStart(); 
  this->revertToLastCommit();   
} 
 
Masonryt::Masonryt(void): 
UniaxialMaterial(0, MAT_TAG_Masonryt), Fm(0),  Ft(0), Um(0),  Uult(0),   Ucl(0), Emo(0),   Length(0), Area1(0), Area2(0),  D1(0), 
D2(0), Ach(0),  Are(0),  Ba(0), Bch(0),  Gun(0),    Gplu(0),  Gplr(0),  Exp1(0),   Exp2(0), IENV(0) 
{ 
 } 
Masonryt::~Masonryt(void) 
{ 
  // Does nothing 
} 
/**************************************************************************** 
 *  PART II : Material State determination --> PROGRAM CORE                                                          * 
 ****************************************************************************/ 
int 
Masonryt::setTrialStrain(double strain, double strainRate) 
{ 
  double  A1 = Emo*Um/Fm; 
  double  A2 = 1 - A1*Um/Uult; 
  //double  Ko = Area1*Emo/Length; 
  double Ko = Area1*Emo;
  double  DeltaU, Kfactor; 
  
  // retrieve history variables 
  Uun    = cUun; 
  Sun    = cSun; 
  Eun    = cEun; 
  Ure    = cUre; 
  Sre    = cSre; 
  Ere    = cEre; 
  Uch    = cUch; 
  Sch    = cSch; 
  Ech    = cEch; 
  U1     = cU1; 
  S1     = cS1; 
  E1     = cE1; 
  U2     = cU2; 
  S2     = cS2;  
  E2     = cE2; 
  UunInt = cUunInt; 
  UreInt = cUreInt; 
  Upl    = cUpl; 
  FtRed  = cFtRed; 
  RuleNo = cRuleNo; 
  IVIR   = cIVIR; 
  InnerCycleNo = cInnerCycleNo; 
  double Uma = Ft / Emo;
   
  // calculate current deformation etc 
  D      = strain; 
  //U      = D/Length; 
  U = D;
  DeltaU = U - cU; 
  if (fabs(DeltaU) <= DBL_EPSILON)    // ignore trivial strain change 
  { 
    S  = cS; 
    Et = cEt; 
  } 
  else 
  { 
// ------------------------------------------------------------------------- 
// Find case & Calculate Stress - Tangent:  
//-------------------------------------------------------------------------- 
this-> Stress_Tangent(U, DeltaU, cU, cS, cEt,  
             Um, Fm, Emo, Ft, Uult, Ucl, Ach, Are, 
           Ba, Bch, Gun, Gplu, Gplr, Exp1, Exp2,  
       U1, S1, E1, U2, S2, E2, S, Et,  
       FtRed, Upl, UunInt, UreInt, Uun, Sun, 
           Eun, Ure, Sre, Ere, Uch, Sch, Ech, 
           RuleNo, InnerCycleNo, IVIR ,UMAXIMA, INDIC, Uma); 
  } 
// CALCULATE FORCE & STIFFNESS (according to the level of the axial deformation) 
if ((Area1 == Area2) || (cArea == Area2)) 
{ 
    Area = Area2; 
} 
else 
{ 
    if (D > D1) 
    { 
    Area = Area1; 
    } 
    else if (D < D2) 
    { 
    Area = Area2; 
    } 
    else  
    { 
    Area = Area1-(Area1-Area2)*(D1-D)/(D1-D2); 
    } 
} 
cArea = Area; 
Kfactor = Et*Area/Emo/Area1;    // stiffness coefficient factor 
K = Kfactor*Ko;                               // current stiffness 
F = S*Area;                                      // current force 
return 0; 
 
} // end of Part II 
 
/**************************************************************************** 
 *  PART III : Get  Strain - Stress - Tangent                                                                                                
 ****************************************************************************/ 
 double 
 Masonryt::getStrain(void) 
 { 
   return D; 
 } 
  
 double 
 Masonryt::getStress(void) 
 { 
   return F; 
 } 
  
  double 
 Masonryt::getTangent(void) 
 { 
   return K; 
 } 
 
  double 
Masonryt::getInitialTangent(void) 
{ 
  //return Ko; 
  return Emo;
} 
/****************************************************************************** 
 *  PART IV : 1. Set up commitState(), revertToLastCommit() & revertToStart()  
 *            2. "getCopy" method -> to provide a clone of the material to a   
 *            calling   object (e.g. Element, Fiber, other Material)                                        
 ******************************************************************************/ 
 int 
 Masonryt::commitState(void) 
 { 
   cD      = D; 
   cF      = F; 
   cK      = K; 
   cU      = U; 
   cS      = S; 
   cEt     = Et; 
   cUun    = Uun; 
   cSun    = Sun; 
   cEun    = Eun; 
   cUre    = Ure; 
   cSre    = Sre; 
   cEre    = Ere;    
   cUch    = Uch; 
   cSch    = Sch; 
   cEch    = Ech; 
   cU1     = U1; 
   cS1     = S1; 
   cE1     = E1; 
   cU2     = U2; 
   cS2     = S2; 
   cE2     = E2; 
   cUunInt = UunInt; 
   cUreInt = UreInt; 
   cUpl    = Upl; 
   cFtRed  = FtRed; 
   cRuleNo = RuleNo; 
   cIVIR   = IVIR; 
   cInnerCycleNo = InnerCycleNo; 
   cArea   = Area; 
   cINDIC = INDIC;
   cUMAXIMA = UMAXIMA;
return 0; 
 } 
 int 
 Masonryt::revertToLastCommit(void) 
 { 
   D      = cD; 
   F      = cF; 
   K      = cK; 
   U      = cU; 
   S      = cS; 
   Et     = cEt; 
   Uun    = cUun; 
   Sun    = cSun; 
   Eun    = cEun; 
   Ure    = cUre; 
   Sre    = cSre; 
   Ere    = cEre; 
   Uch    = cUch; 
   Sch    = cSch; 
   Ech    = cEch; 
   U1     = cU1; 
   S1     = cS1; 
   E1     = cE1; 
   U2     = cU2; 
   S2     = cS2; 
   E2     = cE2; 
   UunInt = cUunInt; 
   UreInt = cUreInt; 
   Upl    = cUpl; 
   FtRed  = cFtRed; 
   RuleNo = cRuleNo; 
   IVIR   = cIVIR; 
   InnerCycleNo = cInnerCycleNo; 
   Area   = cArea; 
   UMAXIMA = cUMAXIMA;
   INDIC = cINDIC;

return 0;    
 } 
int 
 Masonryt::revertToStart(void) 
 {   
   cD        = 0.0; 
   cF        = 0.0; 
   cK        = Ko; 
   cU        = 0.0;  
   cS        = 0.0; 
   cEt       = Emo; 
   cUun      = 0.0; 
   cSun      = 0.0; 
   cEun      = Gun*Emo; 
   cUre      = 0.0; 
   cSre      = 0.0; 
   cEre      = Eun; 
   cUch     = 0.0; 
   cSch     = 0.0; 
   cEch     = 0.0; 
   cU1      = 0.0; 
   S1       = 0.0; 
   cE1       = 0.0; 
   cU2      = 0.0; 
    cS2       = 0.0; 
    cE2       = 0.0; 
   cUunInt   = 0.0; 
   cUreInt   = 0.0; 
   cUpl      = 0.0; 
   cFtRed    = Ft; 
   cRuleNo   = 3; 
   cIVIR     = 1; 
   cInnerCycleNo = 1; 
   cArea     = Area1; 
   cINDIC = 0;
   cUMAXIMA = 0;

 
   D        = 0.0; 
   F        = 0.0; 
   K        = Ko; 
  U        = 0.0; 
   S        = 0.0; 
   Et       = Emo; 
   Uun      = 0.0; 
   Sun      = 0.0; 
   Eun      = Gun*Emo; 
   Ure      = 0.0; 
   Sre      = 0.0; 
   Ere      = Eun; 
   Uch      = 0.0; 
   Sch      = 0.0; 
   Ech      = 0.0; 
   U1       = 0.0; 
   S1       = 0.0; 
   E1       = 0.0; 
   U2       = 0.0; 
   S2       = 0.0; 
   E2       = 0.0; 
   UunInt   = 0.0; 
   UreInt   = 0.0; 
   Upl      = 0.0; 
   FtRed    = Ft; 
   RuleNo   = 3; 
   IVIR     = 1; 
   InnerCycleNo = 1; 
   Area     = Area1; 
   Uma = Ft / Emo;
   UMAXIMA = 0.0;
return 0;       
 }  
UniaxialMaterial* 
Masonryt::getCopy(void) 
{ 
  Masonryt *theCopy = new Masonryt(this->getTag(), Fm, Ft, Um, Uult, Ucl,  
                                     Emo, Length, Area1, Area2, D1, D2, Ach, Are,  
                                     Ba, Bch, Gun, Gplu, Gplr, Exp1, Exp2, IENV); 
   
   theCopy->cD      = cD; 
   theCopy->cF      = cF; 
   theCopy->cK      = cK; 
   theCopy->cU      = cU; 
   theCopy->cS      = cS; 
   theCopy->cEt     = cEt; 
   theCopy->cUun    = cUun; 
   theCopy->cSun    = cSun; 
   theCopy->cEun    = cEun; 
   theCopy->cUre    = cUre; 
   theCopy->cSre    = cSre; 
   theCopy->cEre    = cEre;    
   theCopy->cUch    = cUch; 
   theCopy->cSch    = cSch; 
   theCopy->cEch    = cEch; 
   theCopy->cU1     = cU1; 
   theCopy->cS1     = cS1; 
   theCopy->cE1     = cE1; 
   theCopy->cU2     = cU2; 
   theCopy->cS2     = cS2; 
   theCopy->cE2     = cE2; 
   theCopy->cUunInt = cUunInt; 
   theCopy->cUreInt = cUreInt; 
   theCopy->cUpl    = cUpl; 
   theCopy->cFtRed  = cFtRed; 
   theCopy->cRuleNo = cRuleNo; 
   theCopy->cIVIR   = cIVIR; 
   theCopy->cInnerCycleNo = cInnerCycleNo; 
   theCopy->cArea   = cArea; 
   theCopy->cUMAXIMA = cUMAXIMA;
   theCopy->cINDIC = cINDIC;
 
  return theCopy; 
} 
/**************************************************************************** 
 *  PART IV : 1. "sendSelf()", "recvSelf()" methods for parallel processing                                         
 *            & database program(std::min)g  -> inherited from MovableObject                                    
 *            2. "Print()" method -> inherited by TaggedObject                                                                   
 ****************************************************************************/ 
 //1a. 
 int 
 Masonryt::sendSelf(int commitTag, Channel &theChannel) 
 { 
   static Vector data(53); 
   data(0)   = this->getTag(); 
   data(1)   = Fm; 
   data(2)   = Ft; 
   data(3)   = Um; 
   data(4)   = Uult; 
   data(5)   = Ucl; 
   data(6)   = Emo; 
   data(7)   = Length; 
   data(8)   = Area1; 
   data(9)   = Area2; 
   data(10)  = D1; 
   data(11)  = D2; 
   data(12)  = Ach; 
   data(13)  = Are; 
   data(14)  = Ba; 
   data(15)  = Bch; 
   data(16)  = Gun; 
   data(17)  = Gplu; 
   data(18)  = Gplr; 
   data(19)  = Exp1; 
   data(20)  = Exp2; 
   data(21)  = IENV; 
   data(22)  = cD; 
   data(23)  = cF; 
   data(24)  = cK; 
   data(25)  = cU; 
   data(26)  = cS; 
   data(27)  = cEt; 
   data(28)  = cUun; 
   data(29)  = cSun; 
   data(30)  = cEun; 
   data(31)  = cUre; 
   data(32)  = cSre; 
   data(33)  = cEre; 
   data(34)  = cUch; 
   data(35)  = cSch; 
   data(36)  = cEch; 
   data(37)  = cU1; 
   data(38)  = cS1; 
   data(39)  = cE1; 
   data(40)  = cU2; 
   data(41)  = cS2; 
   data(42)  = cE2; 
   data(43)  = cUunInt; 
   data(44)  = cUreInt; 
   data(45)  = cUpl; 
   data(46)  = cFtRed; 
   data(47)  = cRuleNo; 
   data(48)  = cIVIR; 
   data(49)  = cInnerCycleNo; 
   data(50)  = cArea; 
   data(51) = cUMAXIMA;
   data(52) = cINDIC;

   
   if (theChannel.sendVector(this->getDbTag(), commitTag, data)<0) 
   { 
     opserr << "Masonryt::sendSelf() - failed to sendSelf\n"; 
     return -1; 
     }  
return 0;    
 } 
 //1b. 
 int 
 Masonryt::recvSelf(int commitTag, Channel &theChannel,  
                         FEM_ObjectBroker &theBroker) 
 { 
   static Vector data(53); 
   if (theChannel.recvVector(this->getDbTag(), commitTag, data)<0) 
   { 
     opserr << "Masonryt::recvSelf() - failed to recvSelf\n"; 
     return -1; 
   } 
    
   this->setTag(int(data(0))); 
   Fm      = data(1); 
   Ft      = data(2); 
   Um      = data(3); 
   Uult    = data(4); 
   Ucl     = data(5); 
   Emo     = data(6); 
   Length  = data(7); 
 Area1   = data(8); 
 Area2   = data(9); 
 D1      = data(10); 
 D2      = data(11);    
 Ach     = data(12); 
   Are     = data(13); 
   Ba      = data(14); 
   Bch     = data(15); 
   Gun     = data(16); 
   Gplu    = data(17); 
   Gplr    = data(18); 
   Exp1    = data(19); 
   Exp2    = data(20); 
   IENV    = int(data(21)); 
   cD      = data(22); 
   cF      = data(23); 
   cK      = data(24);   
   cU      = data(25); 
   cS      = data(26); 
   cEt     = data(27); 
   cUun    = data(28); 
   cSun    = data(29); 
   cEun    = data(30); 
   cUre    = data(31); 
   cSre    = data(32); 
   cEre    = data(33); 
   cUch    = data(34); 
   cSch    = data(35); 
   cEch    = data(36); 
   cU1     = data(37); 
   cS1     = data(38); 
   cE1     = data(39); 
   cU2     = data(40); 
   cS2     = data(41); 
   cE2     = data(42); 
   cUunInt = data(43); 
   cUreInt = data(44); 
   cUpl    = data(45); 
   cFtRed  = data(46); 
   cRuleNo = int(data(47)); 
   cIVIR   = int(data(48)); 
   cInnerCycleNo = int(data(49)); 
   cArea   = data(50); 
   cUMAXIMA = data(51);
   cINDIC = int(data(52));
 
   D  = cD; 
   F  = cF; 
   K  = cK; 
   U  = cU; 
   S  = cS; 
   Et = cEt; 
return 0; 
 } 
  //2. 
  void 
  Masonryt::Print(OPS_Stream &s, int flag) 
  { 
    s << "Masonryt:: tag: " << this->getTag() << endln; 
//    s << "Masonryt::(deformation, force, tangent) " << D << " " << F << " " << K << endln; 
	s<< "Compressive stress  f'm= "<<Fm<<endln;
	s<< "Tensile stress      ft = "<<Ft<<endln;
	s<< "Def. at max. Stress eo = "<<Um<<endln;
	s<< "Ultimate deform.    eu = "<<Uult<<endln;
	s<< "Closing deform.     ecl= "<<Ucl<<endln;
	s<< "Elasticity modulus Emo = "<<Emo<<endln;
	s<< "% of Area           A2 = "<<Area2/Area1<<endln;
	s<< "Initial def. Area   e1 = "<<D1<<endln;
	s<< "Final def. Area     e2 = "<<D2<<endln;
	s<< "Ach                    = "<<Ach<<endln;
	s<< "Are                    = "<<Are<<endln;
	s<< "Ba                     = "<<Ba<<endln;
	s<< "Bch                    = "<<Bch<<endln;
	s<< "Gun                    = "<<Gun<<endln;
	s<< "Gplu                   = "<<Gplu<<endln;
	s<< "Gplr                   = "<<Gplr<<endln;
	s<< "Exp1                   = "<<Exp1<<endln;
	s<< "Exp2                   = "<<Exp2<<endln;
	s<< "Ienv                   = "<<IENV<<endln;
  } 
   
/*************************************************************************** 
*  PART VI : Define functions used in the main core of the program                                                
*          1. Comp_Envlp     : Calculates Stress-Tangent on the                                                              
*                              Compressive envelope, Sargin's equation                                                             
*          2. Unload_Reload  : Calculates Stress-Tangent using the                                                       
*                              Proposed equations (7.7) & (7.8)                                                                            
*          3. Plastic_Strain : Calculates Plastic Strain & Reduced                                                            
*                              tensile strength from eq.(7.19) & (7.35)                                                                
*          4. Stess_Tangent  : Calculates stress & tangent modulus at the                                             
*                              current step                                                                                                                 
****************************************************************************/ 
//1. 
void 
Masonryt::Comp_Envlp (double Uenv, double Um, double Fm, double Emo, double Uult, int IENV, double 
&Senv, double &Eenv ) 
{ 
// Calculate Parameters A1, A2 
 
double A1 = Emo*Um/Fm; 
double A2 = 1 - A1*Um/Uult; 
 
if ((Uenv > Um && IENV==1) || (Uenv > Uult && IENV==2)) 
{ 
    // from Sargin's equations (7.1a) - (7.2a) 
 
    double XX = Uenv/Um; 
     Senv  = Fm*( A1*XX + (A2-1)*pow(XX, 2) )/( 1 + (A1-2)*XX + A2*pow(XX, 2) ); 
     
    double EtNom = (Fm/Um) * (A1 + 2*(A2-1)*XX + (2 - A1 - 2*A2)*pow(XX, 2) ); 
    double EtDen = pow(1 + (A1-2)*XX + A2*pow(XX, 2), 2); 
     
     Eenv = (std::max)(EtNom/EtDen, 0.0); 
}     
else if (IENV==1) 
{ 
    // Parabola for the descending branch (Eq. 7.1b &  7.2b) 
 
     Senv = (std::min)( Fm*( 1- pow((Uenv-Um)/(Uult-Um),2) ) , 0.0); 
     Eenv = - 2*Fm*(Uenv - Um)/pow((Uult - Um),2); 
} 
else  
{ 
     Senv = 0.0; 
     Eenv = 0.0; 
}   
return; 
} //End of (1) Comp_Envlp 
//--------------------------------------------------------------------------------------------------------------------------
//2. 
void 
Masonryt::Unload_Reload( double U, double U1, double  U2, double S1, double S2, double E1, double E2,
double &S, double &Et ) 
{ 
  double Esec = (S2 - S1)/(U2 - U1);  
  double B1 = E1/Esec; 
    double B3 = 2 - E2/Esec*(1+B1);    // Parameters from eq.(7.9) 
    double B2 = B1 - B3; 
    double XX = (U - U1)/(U2-U1);      // Relative strain 
 
    S  = S1 + (S2 - S1)*(B1*XX + pow(XX, 2))/(1 + B2*XX + B3*pow(XX, 2)) ; 
    Et = Esec*((B1 + 2*XX)*(1 + B2*XX + B3*pow(XX, 2)) - (B1*XX + pow(XX, 2))*(B2 + 2*B3*XX))/pow(1 +    
B2*XX + B3*pow(XX, 2), 2); 
return; 
} //End of (2) Unload_Reload 
//------------------------------------------------------------------------------------------------------------------------------------------ 
//3. 
void 
Masonryt::Plastic_Strain( double Uun, double Sun, double Um, double Fm, double Emo, double Ft, double 
Ba, double &Upl, double &FtRed  ) 
{ 
	//GST: se llama AXSS25
   Upl = Uun-(Uun - Ba*fabs(Fm)/Emo)*Sun/(Sun - Ba*fabs(Fm));  
 
     if (Upl > Um && Upl <= 0.0 && FtRed!=0) 
     { 
     FtRed = (std::max)(Ft*(1-Upl/Um), 0.0); 
     } 
     else 
     { 
     FtRed = 0.0; 
     } 
return; 
} //End of (3) Plastic_Strain  
void
Masonryt::TRACCION(double Um, double &Upl, double Ft, double Emo, double &Et, double &S, double Uma, double U, double Ucl, double &UMAXIMA, int &INDIC)
{

	double UDEN = Upl;
	double UUTTT = 20 * Ft / Emo;
	double FE = 5.0;
	double UMAXT = Ft / Emo;
	UUTTT2 = (std::max)(Ucl, FE*UMAXT);


	if (UMAXIMA < UMAXT)
	{
		UMAXIMA = UMAXT;
	}

	 UMAXIMAT = UMAXIMA;
	if (( UUTTT2 - UMAXT) == 0.0)
	{
		UUTTT2 = 0.95*UUTTT2;
	}

	double FM2 = Ft*(UUTTT2 - UMAXIMAT) / (UUTTT2 - UMAXT);

	if (U > 0.0 && U < UMAXT && INDIC == 0)
	{
		Et = Emo;
		S = Et*U;
		INDIC = 1;
	}
	else if (U > 0.0 && U < UUTTT2 && U >= UMAXT)
	{
		double FMT3 = Ft;
		if (UMAXIMAT > UMAXT)
		{
			FMT3 = FM2;
		}
		if (UMAXIMAT == UDEN)
		{
			UMAXIMAT = 1.05*UMAXIMAT;
		}
		Et = FM2 / (UMAXIMAT - UDEN);

		double UU3 = UDEN;
		if (fabs(UDEN) > fabs(Um))
		{
			double UU3 = Um;
		}
		Et = FM2 / (UMAXIMAT - UU3);

		S = Et*(U - UU3);
		double aku = 800 * (U - UMAXT);
		double ST2 = Ft / (1 +sqrt(aku));
		FMT3 = ST2;

		if (S > FMT3)
		{
			S = FMT3;
		}
		INDIC = 1;
		UMAXIMA = U;
		if (UMAXIMAT > UMAXIMA)
		{
			UMAXIMA = UMAXIMAT;
		}
	}
	else if (U > 0.0 && U > UUTTT2)
	{
		S = 1.0e-20;
		Et = 1.0e-20;
	}
	else if (U<0.0 && U >Um)
	{
		  	UMAXT = Ft / Emo;
			UDEN = Upl;
			if (fabs(UDEN) > fabs(Um))
			{
				UDEN = Um;
			}
			double EMAX = FM2 / (UMAXIMAT - UDEN);
			S = EMAX*(U - UDEN);
			Et = EMAX;
			INDIC = 1;

	}
	else if (U < 0.0 && U <= Um)
	{
		S = 1.0e-20;
		Et = 1.0e-20;
	}
	else if (U > 0.0 && U < UMAXT && INDIC == 1)
	{
				UDEN = Upl;
				double EMAX = FM2 / (UMAXIMAT - UDEN);
				S = EMAX*(U - UDEN);
				Et = EMAX;
				INDIC = 1;
	}
	else
	{
		S = 1.0e-20;
		Et = 1.0e-20;
	}
	

	return;
 }










//------------------------------------------------------------------------------------------------------------------------------------------ 
//4. 
void 
Masonryt::Stress_Tangent(double U, double DeltaU, double cU, double cS, double cEt,  
    double Um, double Fm, double Emo, double Ft, double Uult, double Ucl, double Ach, double Are, 
    double Ba, double Bch, double Gun, double Gplu, double Gplr, double Exp1, double Exp2,  
    double &U1, double &S1, double &E1, double &U2, double &S2, double &E2,double &S, double &Et,  
    double &FtRed, double &Upl, double &UunInt, double &UreInt, double &Uun, double &Sun, 
    double &Eun, double &Ure, double &Sre, double &Ere, double &Uch, double &Sch, double &Ech, 
    int &RuleNo, int &InnerCycleNo, int &IVIR, double &UMAXIMA, int &INDIC, double Uma ) 
{ 
 double  Ub, Esec, Eplu, Eplr, Ecl, Senv, Eenv; 
 
// ----------------------------------------------------------------------------------------------------------------------------------------- 
// Find case & Calculate Stress - Tangent:  
//----------------------------------------------------------------------------------------------------------------------------------------- 
// RULE 1 : Skeleton/Envelope curve in compression (Sargin's equations) 
//----------------------------------------------------------------------------------------------------------------------------------------- 
if (RuleNo == 1) 
{ 
    if (DeltaU < 0.0)    // loading in compression  
    { 
//.................................................................................................................... 
//  CASE 1-3: ultimate strain has been reached => go to zero stress state   
//.................................................................................................................... 
        if (U <= Uult) 
        {       
            RuleNo = 3; 
           // Upl = 1000*Um; 
			Upl = U;
           // S = 1.0e-20; 
           // Et = 1.0e-20; 
          return; 
        } 
//.......................................................................... 
//  CASE 1-1: travels along the envelope curve    
//..........................................................................             
        else if  (U > Uult) 
        {    
            this-> Comp_Envlp( U, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
            S  = Senv; 
            Et = Eenv; 
        return; 
        } 
    } 
//............................................................................. 
//  CASE 1-2: unloading from Rule 1 using Rule 2    
//............................................................................. 
    else if (DeltaU >= 0.0) 
    {   
        RuleNo = 2; 
        InnerCycleNo = 1; 
        Uun = cU; 
        UunInt = Uun; 
        Sun = cS; 
        Eun = Gun*Emo; 
        Ure = Uun;  
        // Calculate plastic strain & reduced tensile strength 
        this-> Plastic_Strain(Uun, Sun, Um, Fm, Emo, Ft, Ba, Upl, FtRed); 
        // Calculate the parameters needed to implement Rule 2 
        U1 = cU; 
        S1 = cS; 
        E1 = Eun;   // when unloading from Rule 1, E1=Eun=Gun*Emo; 
        U2 = (std::max)( U1 - 1.5*S1/E1, Upl) ;  
        S2 = 0.0; 
        Esec = (S2 - S1)/(U2 - U1);    
        Eplu = Gplu*Emo/pow((1.0 + Uun/Um),Exp1);   // Eq.(7.21) 
        E2 = (std::min)(Eplu, 0.5*Esec);     
        // Calculate Stress & Tangent modulus using proposed equations (7.7) & (7.8) 
        this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2, S, Et ); 
       return; 
     } 
 }// end of RuleNo==1 
//------------------------------------------------------------------------------------------------------------------------------------------ 
// RULE 2 : Unloading from rule 1, 4 or 5   
//------------------------------------------------------------------------------------------------------------------------------------------ 
else if (RuleNo == 2) 
{ 
     if (U >= U2) // unloading enters tensile-stress state or zero stress state 
     { 
//.............................................................................. 
//  CASE 2-6: unloading from Rule 2 using Rule 6    
//.............................................................................. 
       if (FtRed > 0.0) 
       { 
       //    RuleNo = 6;  //GST: R=3
       //    U1 = Upl; 
		   //GST: call traccion tengo yo
      //     Et = Emo*FtRed/Ft; 
      //     S = Et*(U-U1); 
     //      if (S >=FtRed) 
     //       {  
      //          S = FtRed; 
      //          FtRed = 0; 
  //}  
		   RuleNo = 3;
		   Upl = U;
		   this->TRACCION(Um,Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
         return; 
     } 
//............................................................................. 
//  CASE 2-3: unloading from Rule 2 using Rule 3    
//............................................................................. 
       else if (FtRed == 0.0) 
      { 
          RuleNo = 3; 
		  Upl = U;
		  this->TRACCION(Um,Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
 //          S = 1.0e-20; 
 //          Et = 1.0e-20; 
return; 
       } 
    } 
    else if (U < U2) 
    {  
//...................................................................................... 
//  CASE 2-2: Typical unloading - travelling along Rule 2  
//...................................................................................... 
        if (DeltaU > 0.0) 
        {  
           // Calculate Stress & Tangent modulus using proposed equations (7.7) & (7.8) 
            this -> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et ); 
return; 
        }       
//            
//    Cases for RELOADING from RULE 2  before reaching apl (inner cycle)   
// 
        else if (DeltaU <= 0.0) 
        {  
        // Calculate new Reloading point: 
 
            if (UunInt < Upl) 
            { 
              UreInt = (std::min)(cU, Upl);  // strain at the inner loop at which reloading occurs  
              Ure = Ure + Are*(UunInt - UreInt)/ pow(InnerCycleNo, Exp2); // From Eq.(7.32) 
              InnerCycleNo = InnerCycleNo + 1; 
            } 
        // Reloading point at the envelope curve => Calculate stress-tangent: 
        this-> Comp_Envlp( Ure, Um, Fm, Emo, Uult, IENV, Senv, Eenv );  
         Sre = Senv; 
         Ere = Eenv; 
 
        // Calculate new Change point: 
        Ub  = Upl + Ach*(Uun - Sun/Eun - Upl);  // From Eq.(7.24) 
        Ech = Sun/(Uun - Ub);                   // From Eq.(7.25) 
        Sch = Bch*Sre/ pow(InnerCycleNo, 0.4);  // From Eq.(7.33) 
        Sch = (std::min)( Sch, 0.5*Sre);        // limit Eq.(7.34) 
        Uch = Ub + Sch/Ech;                     // From Eq.(7.23) 
 
//................................................................................................................................. 
//  CASE 2-3: Reloading from Rule 2 to zero-stress after Uult has been reached 
//................................................................................................................................. 
            if (Sre == 0.0) 
            {     
                //S = 1.0e-20; 
                //Et = 1.0e-20; 
                RuleNo = 3; 
               //Upl = 1000*Um; 
				Upl = U;
return; 
            } 
            else 
            { 
//............................................................................. 
//  CASE 2-4: Reloading from Rule 2 using Rule 4 
//............................................................................. 
                if (cS >= 0.9*Sch && cU >= 0.9*Uch)      // Case (A) for reloading from Rule 2, page205 
                { 
                    RuleNo = 4; 
                    U1 = Uch; 
                    S1 = Sch; 
                    E1 = Ech; 
                    S2 = cS; 
                    U2 = cU; 
                    E2 = (std::min)(1.2*cEt, 0.9*(S2-S1)/(U2-U1)); 
                } 
//............................................................................. 
//  CASE 2-5: Reloading from Rule 2 using Rule 5 
//............................................................................. 
                else                                    // Case (B) for reloading from Rule 2, page205 
                { 
                    RuleNo = 5; 
                    U1 = cU; 
                    S1 = cS; 
                    U2 = Ure; 
                    S2 = Sre; 
                    E1 = (std::min)( (std::max)(2.0*cEt, 1.2*(S2-S1)/(U2-U1)), Gun*Emo); 
                    E2 = Ere; 
                    if (E2<=0) 
                    { 
                       E2 = 0.5*(S2-S1)/(U2-U1); 
                    } 
                } 
            // Calculate stress-tangent for reloading cases : 
            this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et ); 
return; 
        } 
     } 
   } 
} // end of RuleNo==2 
//----------------------------------------------------------------------------------------------------------------------------- ------------ 
// RULE 3 : Zero stress state  
//----------------------------------------------------------------------------------------------------------------------------------------- 
else if (RuleNo == 3) 
{ 
       if (DeltaU >= 0.0)      // unloading from zero stress in tension 
       { 
//....................................................................................... 
//  CASE 3-6: Unloading - Passing from Rule 3 to Rule 6  
//....................................................................................... 
          // if (FtRed > 0.0) 
          // { 
			   this->TRACCION(Um,Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
			   // GST: call traccion
             //  RuleNo = 6; 
             //  U1 = 0.0; 
             //  Et = Emo; 
             //  Et = Et*FtRed/Ft; 
             //  S = Et*(U-U1); 
                
             //  if (S >= FtRed) 
              // { 
              //      S = FtRed; 
              //    FtRed = 0.0; 
              //    RuleNo = 3; 
             // } 
         return; 
         } 
//....................................................................................... 
//  CASE 3-3un: Unloading - travelling along Rule 3 
//....................................................................................... 
         //  else //if (FtRed == 0.0) 
          // { 
          //     S = 1.0e-20; 
          //     Et = 1.0e-20; 
//return; 
  //         } 
       //} 
       else if (DeltaU < 0.0) // reloading from zero stress in compression 
       { 
           if ((U <= Upl) || (U <= Ucl)) 
           { 
               if (IVIR == 1)               
               {        
//.................................................................................................... 
//  CASE 3-1: Reloading from Rule 3 to Rule 1 (for the 1st time) 
//.................................................................................................... 
                   if (U <= 0)  
                   { 
                       RuleNo = 1; 
                       IVIR = 0; 
                       this-> Comp_Envlp( U,Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
                       S = Senv; 
                       Et = Eenv; 
return; 
                   } 
//................................................................................................................................ 
//  CASE 3-4i: Reloading from Rule 3 using Rule 4 (when Sargin's not used yet) 
//................................................................................................................................                    
                   else if (U > 0)  // goes to compression Aa<0 without previous compressive-cycle (IVIR=1) 
                   {                        // => Rule 4 with special considerations (see page 210 Crisafulli phD) 
                       RuleNo = 4; 
                       U1 = -Ft/Emo; 
                       // Calculate stress & tangent for strain U1 at the Envelope curve 
                       this-> Comp_Envlp( U1, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
                       S1 = Senv; 
                       E1 = Eenv; 
 
                       // Define point 2 (see page 210 Crisafulli phD) 
                       U2 = U; 
                       S2 = 0.0; 
                       E2 = 0.01*Emo;  
 
                       // Calculate stress-strain at current step according to the proposed equations (7.7) & (7.8) 
                       Esec = (S2-S1)/(U2-U1);   
                       this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2, S, Et ); 
return; 
                    } 
                } 
               else //if (IVIR == 0) 
               { 
               // Calculate New Reloading point:    
 
                   if (UunInt < Upl) 
                   { 
                     UreInt = (std::min)(cU, Upl);  
                     Ure = Ure + Are*(UunInt - UreInt)/ pow(InnerCycleNo, Exp2); // From Eq.(7.32) 
                     InnerCycleNo = InnerCycleNo + 1; 
                   } 
                   this-> Comp_Envlp( Ure, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
                   Sre = Senv; 
                   Ere = Eenv; 
 
                  // Calculate New Change point: 
                   Ub  = Upl + Ach*(Uun - Sun/Eun - Upl);  // From Eq.(7.24) 
                   Ech = Sun/(Uun - Ub);                   // From Eq.(7.25) 
                   Sch = Bch*Sre/pow(InnerCycleNo, 0.4);   // From Eq.(7.33) 
                   Sch = (std::min)( Sch, 0.5*Sre);        // limit Eq.(7.34) 
                   Uch = Ub + Sch/Ech;                     // From Eq.(7.23) 
//........................................................................................................................................ 
//  CASE 3-4ii: Reloading from Rule 3 using Rule 4 (when Sargin's already used before) 
//........................................................................................................................................                       
                        if (Sre < 0.0) 
                        { 
                            RuleNo = 4; 
 
                            // Calculate parameters for reloading from Rule 3 
                            U1 = Uch; 
                            S1 = Sch; 
                            E1 = Ech; 
                            S2 = (std::min)(cS, 0.0); 
                            U2 = cU; 
                            Eplu = Gun*Emo/pow(1+Uun/Um, Exp1); // Eq.(7.21)- Tangent modulus corresponding to the plastic strain of the unloading curve  
                            Eplr = Gplr*Eplu;                   // Eq.(7.26)- Reloading modulus from Rule 3  
                            Esec = (S2 - S1)/(U2 - U1); 
 
                            if (Upl > Ucl) 
                            { 
                              E2 = Eplr; 
                            } 
                            else  
                            { 
                              Ecl = 0.15*Eplr; // & now consider cases of Eq.(7.39) 
                               if (Ucl>=U2 && U2>0.0) 
                               { 
                                 E2 = Ecl; 
                               } 
                               else 
                               { 
                                 E2 = Ecl + fabs(U2/Upl)*(Eplr - Ecl); 
                               } 
                             } 
                            double E2i = (std::min)(0.2*E1 , 0.9*Esec); 
                            E2 = (std::min)(E2i, E2);   // Limiting conditions in order to obtain a smooth transition curve 
                             
                           // Calculate stress-strain at current step according to proposed equations (7.7) & (7.8) 
                           this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et ); 
return; 
                       } 
//.................................................................................... 
//  CASE 3-3re: Reloading from Rule 3 using Rule 3 
//....................................................................................                         
                        else  
                        { 
                            RuleNo = 3;    // continues in compression under zero-stress 
                        //    Upl = 1000*Um; 
                        //    S = 1.0e-20; 
                        //    Et = 1.0e-20; 
							Upl = U;
							this->TRACCION(Um, Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
return; 
                        } 
                } 
           } 
//.................................................................................... 
//  CASE 3-3zs: From Rule 3 to zero-stress condition 
//.................................................................................... 
           else //if (U > Ucl || U > Upl) // zero stress condition 
           { 
            //  S = 1.0e-20; 
            //  Et = 1.0e-20; 
			   this->TRACCION(Um, Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
return; 
           } 
      } 
 } //end of RuleNo==3 
//----------------------------------------------------------------------------------------------------------------------------- -------------- 
// RULE 4 : First part of the reloading curve  
//----------------------------------------------------------------------------------------------------------------------------- --------------    
else if (RuleNo == 4) 
{          
    if (DeltaU < 0.0) // Reloading in compression 
    {                                                   
       if (U <= U1) 
       { 
//.................................................................................... 
//  CASE 4-5: Reloading from Rule 4 using Rule 5 
//.................................................................................... 
            if (IVIR == 0)     
            {  
                RuleNo = 5;    
                U2 = Ure; 
                S2 = Sre; 
                E2 = Ere; 
  if (E2<=0) 
  { 
                    E2 = 0.5*(S2-S1)/(E2-E1); 
  } 
                // Calculate stress-strain at current step according to proposed equations (7.7) & (7.8) 
                this->Unload_Reload( U, U1, U2, S1, S2, E1, E2, S, Et ); 
return; 
            }        
//................................................................................................................... 
//  CASE 4-1: Reloading from Rule 4 using Rule 1 (1st time Sargin is used) 
//.................................................................................................................. 
            else if (IVIR == 1)  // enters Rule(1) for the 1st time  
            { 
                RuleNo = 1; 
                IVIR = 0; 
                this-> Comp_Envlp( U, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
                S = Senv; 
                Et = Eenv; 
return; 
            } 
        } 
//.................................................................................... 
//  CASE 4-4: Travelling along Rule 4  
//.................................................................................... 
        else if (U > U1)  // travels along Rule (4) 
    { 
           this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et ); 
return; 
        } 
   } 
//............................................................................ 
//  CASE 4-2: Unloading from Rule 4 using Rule 2 
//............................................................................ 
    else if (DeltaU >= 0.0) // Unloading from Rule (4) using Rule (2) 
    { 
        RuleNo = 2; 
        UunInt = cU; 
        E1 = (std::min)(2*cEt, Gun*Emo); 
 
        // Calculate parameters in order to implement Rule 2 
        U1 = cU; 
        S1 = cS; 
        U2 = (std::max)( U1 - 1.5*S1/E1, Upl);   
        S2 = 0.0; 
        Upl = U2; 
        Esec = (S2 - S1)/(U2 - U1); 
        Eplu = Gplu*Emo/pow(1.0 + Uun/Um , Exp1); // Eq.(7.21) 
        E2 = (std::min)(Eplu, 0.5*Esec); 
 
        // Calculate stress-tangent at the current step according to proposed equations (7.7) & (7.8) 
        this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2, S, Et ); 
return; 
    } 
 } // end of RuleNo==4 
//----------------------------------------------------------------------------------------------------------------------------- ----------- 
// RULE 5 : Second part of the reloading curve  
//----------------------------------------------------------------------------------------------------------------------------- -----------  
else if (RuleNo == 5) 
{ 
//.......................................................................................................... 
//  CASE 5-1: Reloading from Rule 5 to the envelope curve (Rule 1) 
//.......................................................................................................... 
    if (U <= Ure)  
    { 
        RuleNo = 1; 
        this-> Comp_Envlp( U, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
        S = Senv; 
        Et = Eenv; 
return; 
    }        
    else if (U > Ure) 
    { 
//............................................................................ 
//  CASE 5-2: Unloading from Rule 5 using Rule 2 
//............................................................................         
        if (DeltaU > 0.0) // Unloading from Rule (5) using Rule (2) 
        { 
            RuleNo = 2; 
            UunInt = cU; 
            if (Uun > UunInt) 
            { 
                Uun = cU; 
                Ure = Uun; 
                Sun = cS; 
                InnerCycleNo = 1; 
 
                // Calculate new plastic strain & reduced tensile strength 
                this-> Plastic_Strain(Uun, Sun, Um, Fm, Emo, Ft, Ba, Upl, FtRed); 
            } 
 
            Eun = Gun*Emo; 
            E1 = Eun; 
 
            // Define points 1 & 2 for Rule 2 
            U1 = cU; 
            S1 = cS; 
            U2 = (std::max)(U1-1.5*S1/E1, Upl); 
            S2 = 0.0; 
            Esec = (S2 - S1)/(U2 - U1); 
            Eplu = Gplu*Emo/pow(1.0 + Uun/Um, Exp1); // Eq.(7.21) 
            E2 = (std::min)(Eplu, 0.5*Esec); 
        } 
//.......................................................................... 
//  CASE 5-5: Unloading travelling along Rule 5  
//..........................................................................              
   //elseif DeltaU < 0.0,  travels along Rule (5) with points 1,2 as already defined, end 
    
  // Calculate stress-tangent at the current step according to proposed equations (7.7) & (7.8) 
     this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et ); 
return; 
    } 
} // end of RuleNo==5 
//------------------------------------------------------------------------------------------------------------------------------------------ 
// RULE 6 : Tensile behavior (elastic)  
//-------------------------------------------------------------------------------------------------------------- ----------------------------   
else if (RuleNo==6) 
{ 
 //   Et = Emo*FtRed/Ft; 
 //   S  = Et*(U-U1);  
//.......................................................................... 
//  CASE 6-3: Reloading from Rule 6 to Rule 3  
//..........................................................................     
  //  if (S >= FtRed) 
  //  { 
  //      RuleNo = 3; 
  //      S = FtRed; 
  //      FtRed = 0; 
	this->TRACCION(Um,Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
return; 
    } 
if (S > 0.0)
{
	RuleNo = 3;
}
//.......................................................................... 
//  CASE 6-1: Reloading from Rule 6 to Rule 1  
//..........................................................................  
    else if (S < 0.0 && IVIR == 1 && U <= 0.0) 
    { 
        RuleNo = 1; 
        IVIR = 0; 
        this-> Comp_Envlp( U, Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
        S = Senv; 
        Et = Eenv;   
return; 
    } 
    else if (S < 0.0) 
    { 
        // Define new reloading point 
 
        if (UunInt < Upl) 
        { 
          UreInt = (std::min)(cU, Upl); // strain at the inner loop at which reloading occurs  
          Ure = Ure + Are*(UunInt - UreInt)/pow(InnerCycleNo, Exp2); // From Eq.(7.32) 
          InnerCycleNo = InnerCycleNo + 1; 
        } 
        // Calculate reloading strain & tangent 
 
        this-> Comp_Envlp( Ure,Um, Fm, Emo, Uult, IENV, Senv, Eenv ); 
        Sre = Senv; 
        Ere = Eenv; 
 
        // Define new Change point  
 
        Ub  = Upl + Ach*(Uun - Sun/Eun - Upl);         // From Eq.(7.24) 
        Ech = Sun/(Uun - Ub);                          // From Eq.(7.25) 
        Sch = Bch*Sre/pow(InnerCycleNo, 0.4);          // From Eq.(7.33) 
        Sch = (std::min)( Sch, 0.5*Sre);               // limit Eq.(7.34) 
        Uch = Ub + Sch/Ech;                            // From Eq.(7.23) 
//.......................................................................... 
//  CASE 6-4: Reloading from Rule 6 using Rule 4  
//..........................................................................         
        if (Sre < 0.0)  
        { 
            RuleNo = 4; 
            // Define points 1 & 2 (SS27) 
            U1 = Uch; 
            S1 = Sch; 
            E1 = Ech; 
            S2 = (std::min)(cS, 0.0); 
            U2 = cU; 
            Eplu = Gun*Emo/pow(1+Uun/Um, Exp1);   // Eq.(7.21) 
            Eplr = Gplr*Eplu;      // Eq.(7.26) 
            Esec = (S2 - S1)/(U2 - U1); 
 
            if (Upl > Ucl) 
            { 
              E2 = Eplr; 
            } 
            else  
            { 
              Ecl = 0.15*Eplr; // & now consider cases of Eq.(7.39) 
              if (Ucl>=U2 && U2>0.0) 
              { 
                E2 = Ecl; 
              } else { 
                E2 = Ecl + fabs(U2/Upl)*(Eplr - Ecl); 
              } 
            } 
           double E2i = (std::min)(0.2*E1, 0.9*Esec); 
           E2 = (std::min)(E2i, E2); // Limiting conditions in order to obtain a smooth transition curve  
  
        // Calculate stress-tangent at the current step according to proposed equations (7.7) & (7.8) 
 
        Esec = (S2-S1)/(U2-U1); 
        this-> Unload_Reload( U, U1, U2, S1, S2, E1, E2,  S, Et );  
return; 
       }   
//.......................................................................... 
//  CASE 6-3: Reloading from Rule 6 using Rule 3  
//..........................................................................             
        else //if (Sre >= 0.0)  // accumulates strains under zero-stress 
        { 
            RuleNo = 3; 
           // Upl = 1000*Um; 
           // S = 1.0e-20; 
           // Et = 1.0e-20; 
			Upl = U;
			this->TRACCION(Um,Upl, Ft, Emo, Et, S, Uma, U, Ucl, UMAXIMA, INDIC);
return; 
        } 
return; 
    } 
 } // end of RuleNo==6 
 

