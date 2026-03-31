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

/* 
** Description: This file contains the class implementation for a       
** Constitutive (FEM) model for FRP and Steel-Confined Concrete for 
** Circular Concrete Sections  **Created in May 2011**                 
** written and developed by:                                        
** Konstantinos G. Megalooikonomou (C)                              
*/ 

#include <FRPConfinedConcrete.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <cmath>
#include <cstdio>
#include <float.h>
#include <iostream>
#include <algorithm>
#include <cctype>

static double fpc,Ec, Ec1, Ec2, R, A, Rcore, Acore, Acover, beta1, beta2, Ash, rs, eyh;
static const double pi = 3.1415926;
static double min(double a, double b);
static int numFRPConfinedConcrete = 0;

// NOTE: units should b in Newton(N) and MegaPascal(MPa)

void *
OPS_FRPConfinedConcrete(void)
{
  if (numFRPConfinedConcrete == 0) {
    numFRPConfinedConcrete++;
    opserr << "FRPConfinedConcrete uniaxial material - Developed by Konstantinos G. Megalooikonomou University of Roma Tre Copyright 2009";
  }

  opserr << "Due to known issues and unreliable results, this material has been" << endln;
  opserr << "temporarily removed from the compiled versions of OpenSees (Tcl and Py)" << endln;
  opserr << "The material source code remains available. Compile at your own risk." << endln;
  return 0;

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[18];  // size of arg list
  int numData;

  if (OPS_GetNumRemainingInputArgs() != 19) {
	  opserr << "WARNING invalid #args: uniaxialMaterial FRPConfinedConcrete $tag $fpc1 $fpc2 $epsc0";
	  opserr << " $D $c $Ej $Sj $tj $eju $S $fyl $fyh $dlong $dtrans $Es $v0 $k $useBuck\n";
	  return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial FRPConfinedConcrete tag" << endln;
    return 0;
  }


  numData = 18;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
	  opserr << "WARNING invalid Material Properties: fpc1: Concrete Core Compressive Strength \n";
	  opserr << "fpc2: Concrete Cover Compressive Strength \n";
	  opserr << "epsc0: Strain Corresponding to Unconfined Concrete Strength \n";
	  opserr << "D = Diameter of the Circular Section \n";
	  opserr << "c = concrete cover \n";
	  opserr << "Ej = Elastic Modulus of the Jacket \n";
	  opserr << "Sj = Clear Spacing of the FRP strips - zero if it's continuous \n";
	  opserr << "tj = Thickness of the FRP Jacket\n";
	  opserr << "eju = Rupture strain of the Jacket\n";
	  opserr << "S = Spacing of the stirrups\n";	  
	  opserr << "fyl = Yielding Strength of longitudinal steel bars\n";
	  opserr << "fyh = Yielding Strength of the hoops\n";
	  opserr << "dlong = Diameter of the longitudinal bars\n";
	  opserr << "dtrans = diameter of the stirrups\n";
	  opserr << "Es = Steel's Elastic modulus\n";
	  opserr << "vo = Poisson's coefficient for concrete\n";
      opserr << "k = reduction factor (0.5-0.8) for the rupture strain of the FRP\n";
	  opserr << "useBuck = FRP Jacket Failure Criterion due to Buckling of Longitudinal Compressive Steel Bars (0 = not include it, 1= to include it)\n";
    return 0;	
  }

    theMaterial = new FRPConfinedConcrete(iData[0],dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],
		dData[7],dData[8],dData[9],dData[10],dData[11],dData[12],dData[13],dData[14],dData[15],dData[16],dData[17]);    

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type FRPConfinedConcrete\n";
    return 0;
  }

  return theMaterial;
}

/*
Input:
D = Diameter of the Circular Section, c = concrete cover, Ej = Elastic Modulus of the Jacket, Sj = Clear Spacing of the FRP strips - zero if it's continuous, tj = Thickness of the FRP Jacket, eju = Rupture strain of the Jacket,
S = Spacing of the stirrups, fyl = Yielding Strength of longitudinal steel bars, fyh = Yielding Strength of the hoops, dlong = Diameter of the longitudinal bars, dtrans = diameter of the stirrups
Es = Steel's Elastic modulus, vo = Poisson's coefficient for concrete, k = reduction factor (0.5-0.8) for the rupture strain of the FRP, useBuck = FRP Jacket Failure Criterion due to Buckling of Longitudinal Compressive Steel Bars (0 = not include it, 1= to include it).
*/

 FRPConfinedConcrete::FRPConfinedConcrete(int tag, 
	 double fpc1_,
	 double fpc2_, 
	 double epsc0_, 
	 double D_, 
	 double c_, 
	 double Ej_, 
	 double Sj_, 
	 double tj_, 
	 double eju_, 
	 double S_, 
	 double fyl_,
	 double fyh_, 
	 double dlong_, 
	 double dtrans_, 
	 double Es_, 
	 double v0_, 
	 double k_,
	 double useBuck_)
 :UniaxialMaterial(tag, MAT_TAG_FRPConfinedConcrete),CminStrain(0.0), CendStrain(0.0),Cstrain(0.0), Cstress(0.0), CaLatstress(0.0) ,
   CbLatstress(0.00001),CLatStrain(0.0) ,CConvFlag(false) ,CConfRat(1.0) ,CConfStrain(epsc0),CLBuck(0.0)
{
  fpc1 = fpc1_;
  fpc2 = fpc2_;
  epsc0 = epsc0_;
  D=D_;
  c=c_;
  Ej=Ej_;
  Sj=Sj_;
  tj = tj_;
  eju = eju_;
  S=S_;
  fyl = fyl_;
  fyh = fyh_;
  dlong = dlong_;
  dtrans = dtrans_;
  Es = Es_;
  v0=v0_;
  k = k_;
  useBuck = useBuck_;

   double Ec0;
   R=D/2.;
   A=pi*pow(R,2);
   //Regions
   //Core Region
   Rcore=R-c;
   Acore=pi*pow(Rcore,2);
   //Cover Region
   Acover=A-Acore;

   //Concrete
   fpc  = (Acore/A)*fpc1 + (Acover/A)*fpc2;
   beta1= 5700.0/(sqrt(fpc1))-500;
   beta2= 5700.0/ (sqrt(fpc2))-500;

   //Steel
   Ash = pi*pow(dtrans,2)/4;
   rs  = (4*Ash)/(S*2*Rcore);
   eyh = fyh/Es;
  
   // Initial tangent
   Ec1 = 5700*sqrt(fpc1);
   Ec2 = 5700*sqrt(fpc2);
   Ec  = (Acore/A)*Ec1 + (Acover/A)*Ec2;
   Ec0 = Ec;
   Ctangent = Ec0;
   CunloadSlope = Ec0;
   Ttangent = Ec0;
  
  // Set trial values
  this->revertToLastCommit();

   // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////

  //bucklong criterion always initialized to false
  buckCrInit = false;
}

FRPConfinedConcrete::FRPConfinedConcrete():UniaxialMaterial(0, MAT_TAG_FRPConfinedConcrete),
 fpc1(0.0),fpc2(0.0), epsc0(0.0),
 CminStrain(0.0), CunloadSlope(0.0), CendStrain(0.0),
 Cstrain(0.0), Cstress(0.0),CaLatstress(0.0) ,CbLatstress(0.00001),CLatStrain(0.0) ,CConvFlag(true) ,CConfRat(1.0),CConfStrain(epsc0),CLBuck(0.0)
{
  // Set trial values
  this->revertToLastCommit();

  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

FRPConfinedConcrete::~FRPConfinedConcrete()
{
  if (SHVs != 0)
	  delete SHVs;
}

double FRPConfinedConcrete::getInitialTangent( ) {return Ec;}

int FRPConfinedConcrete::setTrialStrain (double strain, double strainRate)
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tstrain = Cstrain;
   TaLatstress = CaLatstress;
   TbLatstress = CbLatstress;
   TConvFlag = CConvFlag;
   TConfRat = CConfRat;
   TConfStrain = CConfStrain ;
   TLatStrain = CLatStrain ;
   TLBuck = CLBuck;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    return 0;
  }

  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain < Cstrain) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }

  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }

  //  Tstrain = strain;  
  return 0;
}
int 
FRPConfinedConcrete::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
  // Reset trial history variables to last committed state

  TminStrain = CminStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Tstrain = Cstrain;
  TaLatstress = CaLatstress;
  TbLatstress = CbLatstress;
  TConvFlag = CConvFlag;
  TConfRat = CConfRat;
  TConfStrain = CConfStrain ;
  TLatStrain = CLatStrain ;
  TLBuck = CLBuck;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON) {
    stress = Tstress;
    tangent = Ttangent;
    return 0;
  }

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    stress = 0;
    tangent = 0;
    return 0;
  }

  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain <= Cstrain) {

    TminStrain = CminStrain;
    TendStrain = CendStrain;

    reload ();

    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  //  Tstrain = strain;  
  stress = Tstress;
  tangent =  Ttangent;

  Tstrain = strain;

  return 0;
}

void FRPConfinedConcrete::determineTrialState (double dStrain)
{  
  TminStrain = CminStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*dStrain;
  
  // Material goes further into compression
  if (Tstrain <= Cstrain) {
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
}

void FRPConfinedConcrete::reload ()
{
  if (Tstrain <= TminStrain) {
    
    TminStrain = Tstrain;
    
    // Determine point on envelope
    envelope ();
    unload ();
  }
  else if (Tstrain <= TendStrain) {
    Ttangent = TunloadSlope;
    Tstress = Ttangent*(Tstrain-TendStrain);
  }
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
}
void FRPConfinedConcrete::envelope ( )
{ 
  //  double fla =0.0, flb =0.0;  Delete this line
  double arrayLatA[6],arrayLatB[6],arrayLatC[6];
  double tol = 0.000001 ;//  Add this line
  
  TConvFlag = false;
  double number = 1.0 ;
  bool changedStrain = false;
  
  if (Tstrain < 0.0) {
    Tstrain = -Tstrain ;
    changedStrain = true;
  }
  
  while(!TConvFlag) {
    number = number + 1 ;
    
    flat(TaLatstress,arrayLatA);
    double ya = arrayLatA[0];
    flat(TbLatstress,arrayLatB);
    double yb = arrayLatB[0];
    
    if (yb*ya >0) {
      TbLatstress = 0.1*number; //Change this line
      continue;
    }
    
    double dx = ya*(TaLatstress - TbLatstress)/(ya - yb);
    double flc  = TaLatstress - dx;
    flat(flc, arrayLatC);
    double yc = arrayLatC[0];
    double sigc = arrayLatC[1];
    double flj = arrayLatC[2];
    double fcc = arrayLatC[3];
    double et_cover = arrayLatC[4];
    double el_cover = arrayLatC[5];
    if (yc == 0 || fabs(yc) <= tol) { 
      Tstress=-sigc;
      TConfRat = fcc/fpc;
      TConfStrain = (5*(TConfRat-1)+1)*epsc0;
      TaLatstress = flj;
      TLatStrain  = et_cover;
      double dStrain = Tstrain - Cstrain;
      if (changedStrain == true)
	dStrain = -Tstrain - Cstrain;
      else
	dStrain = Tstrain - Cstrain;

      Ttangent= (Tstress-Cstress)/dStrain;
      TConvFlag = true ;
      return;
    }
    if  (yb*yc > 0) { TbLatstress = TaLatstress; //Change this line
      yb = ya;
      TaLatstress = flc; //Change this line
      ya = yc;
    }
    else{
      TaLatstress = flc; //Change this line
      ya = yc;
    }

    int Maxiter = 100;   //Change this line
    int iter = 0 ;
    
    while (fabs(yc) > tol && iter<=Maxiter) {
      
      dx = ya*(TaLatstress - TbLatstress)/(ya - yb);
      flc  = TaLatstress - dx; //Change this line
      flat(flc, arrayLatC);
      yc = arrayLatC[0];
      sigc = arrayLatC[1];
      flj = arrayLatC[2];
      fcc = arrayLatC[3];
      et_cover = arrayLatC[4];
      el_cover = arrayLatC[5];
      if  (yb*yc > 0) {
	TbLatstress = TaLatstress; //Change this line
	yb = ya;
	TaLatstress = flc;//  Change this line
	ya = yc;
      } else {
	TaLatstress = flc;//  Change this line
	ya = yc;
      }
    }
    if (iter<=Maxiter) {

      Tstress=-sigc;
      TConfRat = fcc/fpc;
	TConfStrain = (5*(TConfRat-1)+1)*epsc0;
	TaLatstress = flj;
	TLatStrain = et_cover;
	double dStrain = 0;
	if (changedStrain == true)
	  dStrain = -Tstrain - Cstrain;
	else
	  dStrain = Tstrain - Cstrain;
	//      Ttangent= fabs((Tstress-Cstress)/dStrain); Delete this line
	Ttangent= (Tstress-Cstress)/dStrain;
	TConvFlag = true ;

    } else {
      TConvFlag = false;
    }
  }
  double et_cover = arrayLatC[4];
  if (et_cover >= k*eju){
    opserr << "FRP Rupture" ; 
  }
  //2017 adds. 

	if(useBuck>0.5){ //0 false, 1 true
		double eyl = fyl/Es;
		double el_cover = arrayLatC[5];
		double vcover = el_cover / Tstrain;
   
		if ( (Tstrain>=eyl) && (vcover >= 0.5 ) ){
			if (!buckCrInit){
				opserr << "Initiation of Buckling of Long.Bar under Compression";
				double Ib = (pi*(pow(dlong,4)))/64;
				double EIred = 0.5*Es*Ib*sqrt(fyl/400);
				double Abar = (pi*(pow(dlong,2)))/4;
				//Ash = pi*dtrans^2/4; Ash Has already been calculated in the constructor function
				double Dcore = D - 2*c;
				int mBuck = 1; //buckling mode for critical load
				double Pcr = fyl*Abar;
				bool convFlag2 = false; //secondary convergence flag
				double n;
				bool convFlag = myRegulaFalsi(Pcr,EIred,Es,Ash,Dcore,S,mBuck,n,convFlag2);
	  
				double LBuck;
				if (!convFlag){
					LBuck = S;
				}
				else{
					if (n<=1){
						LBuck = S;
					}
					else{// if (n>1)
						n = floor(n + 0.5);//round(n)
						LBuck = (n+1)*S;
					}
				}
				TLBuck = LBuck;
				buckCrInit = true;// never going inside again, LBuck computed only once 
			}	//end of lbuck calculation for once

			double esl = Tstrain;     
			double theta = (6.9/(pow((TLBuck/dlong),2)))-0.05;
			double wa = (esl*dlong)/((0.035*cos(theta)+theta)/(cos(theta)-0.035*theta)); 
			double wb = (esl/((0.07*cos(theta)+theta)/(cos(theta)-0.07*theta))+0.035)*dlong;
			double w = std::max(wa,wb);
			double ecbuck;

			if (Sj<0.000001){ //practically zero
				ecbuck =( 2*pi*fabs(((w-c))))/(pi*D);
			}
			else{ //nonzero
				ecbuck = (2*pi*fabs((w-c/2)))/(pi*D);
			}

			if( (ecbuck > 0) && (ecbuck >= eju)){
				opserr << "FRP Rupture due to Buckling of Long.Bar under compression";
			}
		}//end of (Tstrain>=eyl) && (vcover >= 0.5 )
	}//end of useBuck
  
  if (changedStrain == true)
    Tstrain = -Tstrain;
}

bool FRPConfinedConcrete::myRegulaFalsi(double Pcr, double EIred, double Es, double Ash, double Dcore, double S, int mBuck, double& xRes, bool& returnFlag){
	
	int insIter = 0; int insMaxIter = 1000;
	double xA = 0; 
	double xB = 10;
	double xNew;

	double fxA =  PCriticalSolve(xA,Pcr,EIred,Es,Ash,Dcore,S,mBuck);
	double fxAhat =  PCriticalSolve(xB,Pcr,EIred,Es,Ash,Dcore,S,mBuck);

	xNew = xA - (fxA)*(xA - xB)/(fxA - fxAhat);
	double fxNew = PCriticalSolve(xNew,Pcr,EIred,Es,Ash,Dcore,S,mBuck);
	
	while (( fabs(fxNew)> 1E-6) && (insIter<=insMaxIter)){
		insIter++;
		if  (fxAhat*fxNew>0){
			xB = xNew;
			fxAhat = fxNew;
		}
		else{
			xA = xNew;
			fxA = fxNew;
		}
		
		xNew = xA - (fxA)*(xA - xB)/(fxA - fxAhat);
		fxNew = PCriticalSolve(xNew,Pcr,EIred,Es,Ash,Dcore,S,mBuck);

		if ((fabs(xA-xB)<1E-12) && (fxA*fxAhat<0)){
			returnFlag = true;	
			break;
		}
	}


	xRes = xNew;
	if(insIter<=insMaxIter){
		return true;
	}
	else{
		return false;
	}
}

double  
FRPConfinedConcrete::PCriticalSolve(double n, double Pcr, double EIred, double Es, double Ash, double Dcore, double S, int mBuck){//has to become zero
	double f = Pcr-((pow(pi,2)*EIred)/pow(((n+1)*S),2))*(pow(mBuck,2.0) + (1/(pow(mBuck,2.0)))*((((n*Es*Ash)/((n+1)*S*pi*Dcore))*pow(((n+1)*S),4))/(pow(pi,4)*EIred)));
	return f;
}

void 
FRPConfinedConcrete::flat (double flcover_n, double arrayLat[6] )
{

  double v, els, fls, ksl, fls_n, flcore, fcc_core, ecc_core, x_core, Esec_core, r_core, fc_core, fcc_cover, ecc_cover,
	   x_cover, Esec_cover, r_cover, fc_cover, sig, fcc, el_core, el_cover, et_cover, rj, ke, flj, y ;
  //Braga, Gigliotti & Laterza Model

  v=v0*(1+0.2*(Tstrain/epsc0)-pow((Tstrain/epsc0),2.0)+1.55*pow((Tstrain/epsc0),3.0));
  els=v*Tstrain;

  if (els < eyh) {
    fls=(Ec1*Es*Ash*v*Tstrain)/(Rcore*Ec1*S+Es*Ash*(1-v)*(v*Tstrain+1));
  }  else {
    fls = 0.5*rs*fyh;
  }


  ksl= (45*pow((dlong/S),3))/(45*pow((dlong/S),3)+(dtrans/dlong)*(dtrans/(pi*Rcore/2)));
  fls_n=ksl*fls;


  //Spoelstra&Monti iteration
  flcore=flcover_n+fls_n;
  //Vertical Stresses
  fcc_core=fpc1*(2.254*sqrt(1+7.94*(flcore/fpc1))-2*(flcore/fpc1)-1.254); 
  ecc_core=epsc0*(1+5*(fcc_core/fpc1-1));
  x_core = Tstrain /ecc_core;
  Esec_core = fcc_core/ecc_core ;

  r_core  = Ec1/(Ec1-Esec_core);
  fc_core = (fcc_core*x_core*r_core)/(r_core-1+pow(x_core,r_core));
  fcc_cover=fpc2*(2.254*sqrt(1+7.94*(flcover_n/fpc2))-2*(flcover_n/fpc2)-1.254);
  ecc_cover=epsc0*(1+5*(fcc_cover/fpc2-1));
  x_cover = Tstrain /ecc_cover;
  Esec_cover = fcc_cover/ecc_cover ;
  r_cover  = Ec2/(Ec2-Esec_cover);
  fc_cover = (fcc_cover*x_cover*r_cover)/(r_cover-1+pow(x_cover,r_cover));
  sig=(Acore/A)*fc_core + (Acover/A)*fc_cover;
  fcc=(Acore/A)*fcc_core + (Acover/A)*fcc_cover;
  el_core=(Ec1*Tstrain-fc_core)/(2*beta1*fc_core);
  el_cover=(Ec2*Tstrain-fc_cover)/(2*beta2*fc_cover);
  et_cover=(Rcore*(1+el_core)+c*(1+el_cover))/(Rcore+c)-1;
  rj=(4*tj)/D;
  ke=pow((1-(Sj/(2*D))),2); // 14th fib Bulletin
  flj=(0.5*ke*rj*Ej*et_cover);
  y=flj-flcover_n;

  arrayLat[0] = y;
  arrayLat[1] = sig;
  arrayLat[2] = flj;
  arrayLat[3] = fcc;
  arrayLat[4] = et_cover;
  arrayLat[5] = el_cover;
}

double FRPConfinedConcrete::ComputeTendStrain()
{
  
  double tempStrain = TminStrain;
  
  double eta = tempStrain/TConfStrain;
  
  double ratio = 0.707*(eta-2.0) + 0.834;
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  TendStrain = ratio*TConfStrain;
 
  return TendStrain;
}

void FRPConfinedConcrete::unload ( )
{
double  Ec = getInitialTangent ();
  ComputeTendStrain();
  
  double temp1 = TminStrain - TendStrain;

  double fcc = TConfRat*fpc;

  double dstrain = Tstrain - Cstrain;
            
  double x = (Tstrain - dstrain)/(-TConfStrain);
  
  double Esec = fcc/TConfStrain;
		   
  double r = Ec/(Ec-Esec);
		   
  double sig = -(fcc*x*r)/(r-1+pow(x,r));
  
  double Ec0 = sig/(TminStrain-TendStrain);
  
  double temp2 = Tstress/Ec0;
  
  if (temp1 > -DBL_EPSILON) {   // temp1 should always be negative
    TunloadSlope = Ec0;
  }
  else if (temp1 <= temp2) {
    TendStrain = TminStrain - temp1;
    TunloadSlope = Tstress/temp1;
  }
  else {
    TendStrain = TminStrain - temp2;
    TunloadSlope = Ec0;
  }
}

double FRPConfinedConcrete::getStress ()
{
  return Tstress;
}

double FRPConfinedConcrete::getStrain ()
{
   return Tstrain;
}

double FRPConfinedConcrete::getTangent ()
{
  return Ttangent;
}

double FRPConfinedConcrete::getLatstress ()
{
   return TaLatstress;
}

double FRPConfinedConcrete::getLatStrain ()
{
   return TLatStrain;
}

int FRPConfinedConcrete::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CendStrain = TendStrain;
   CunloadSlope = TunloadSlope;
   CbLatstress = TbLatstress;
   CConvFlag = TConvFlag;
   CConfRat = TConfRat;
   CConfStrain = TConfStrain ;
   CLBuck = TLBuck;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;
   CLatStrain = TLatStrain;
   CaLatstress = TaLatstress;

   return 0;
}
int FRPConfinedConcrete::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   TbLatstress = CbLatstress;
   TConvFlag = CConvFlag;
   TConfRat = CConfRat;
   TConfStrain = CConfStrain ;
   TLBuck = CLBuck;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   TLatStrain = CLatStrain;
   TaLatstress = CaLatstress;

   return 0;
}

int FRPConfinedConcrete::revertToStart ()
{ 
  double Ec = 5700*sqrt(fpc);
  double Ec0 = Ec;

   // History variables

   CminStrain = 0.0;
   CendStrain = 0.0;
   CunloadSlope = Ec0;
   CbLatstress = 0.00001;
   CConvFlag = 0.0;
   CConfRat = 1.0;
   CConfStrain = 0.0 ;
   CLBuck = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;
   CLatStrain = 0.0;
   CaLatstress = 0.0;


   // Reset trial variables and state
   this->revertToLastCommit();

      // Quan April 2006---
   if (SHVs !=0) {SHVs->Zero();}
   parameterID=0;

   return 0;
}

UniaxialMaterial* FRPConfinedConcrete::getCopy()
{
   FRPConfinedConcrete* theCopy = new FRPConfinedConcrete(this->getTag(),
							  fpc1, fpc2, epsc0, D, c, Ej, Sj, tj, eju, S, fyl, fyh, dlong, dtrans, Es, v0, k, useBuck);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CunloadSlope = CunloadSlope;
   theCopy->CendStrain = CendStrain;
   theCopy->CbLatstress = CbLatstress;
   theCopy->CConvFlag = CConvFlag;
   theCopy->CConfRat = CConfRat;
   theCopy->CConfStrain = CConfStrain;
   theCopy->CLBuck = CLBuck;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;
   theCopy->TLatStrain = CLatStrain;
   theCopy->TaLatstress = CaLatstress;

   return theCopy;
}

int FRPConfinedConcrete::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(31);
   data(0) = this->getTag();

   // Material properties
   data(1)  = fpc1;
   data(2)  = fpc2;
   data(3)  = epsc0;
   data(4)  = D;
   data(5)  = c;
   data(6)  = Ej;
   data(7)  = Sj;
   data(8)  = tj;
   data(9)  = eju;
   data(10) = S;
   data(11) = fyl;
   data(12) = fyh;
   data(13) = dlong;
   data(14) = dtrans;
   data(15) = Es;
   data(16) = v0;
   data(17) = k;
   data(18) = useBuck;

   // History variables from last converged state
   data(19) = CminStrain;
   data(20) = CunloadSlope;
   data(21) = CendStrain;
   data(22) = CbLatstress;
   data(23) = CConfRat;
   data(24) = CConfStrain;
   data(25) = CLBuck;

   // State variables from last converged state
   data(26) = Cstrain;
   data(27) = Cstress;
   data(28) = Ctangent;
   data(29) = CLatStrain;
   data(30) = CaLatstress;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "FRPConfinedConcrete::sendSelf() - failed to send data\n";

   return res;
}

int FRPConfinedConcrete::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(31);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "FRPConfinedConcrete::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

   // Material properties
   fpc1  = data(1);
   fpc2  = data(2);
   epsc0 = data(3);
   D     = data(4);
   c     = data(5);
   Ej    = data(6);
   Sj    = data(7);
   tj    = data(8);
   eju   = data(9);
   S     = data(10);
   fyl   = data(11);
   fyh   = data(12);
   dlong = data(13);
   dtrans= data(14);
   Es    = data(15);
   v0    = data(16);
   k     = data(17);
   useBuck = data(18);

   // History variables from last converged state
   CminStrain   = data(19);
   CunloadSlope = data(20);
   CendStrain   = data(21);
   CbLatstress  = data(22);
   CConfRat     = data(23);
   CConfStrain  = data(24);
   CLBuck       = data(25);

   // State variables from last converged state
   Cstrain      = data(26);
   Cstress      = data(27);
   Ctangent     = data(28);
   CLatStrain   = data(29);
   CaLatstress  = data(30);

   // Set trial state variables
	  Tstrain = Cstrain;
      Tstress = Cstress;
     Ttangent = Ctangent;
   TLatStrain = CLatStrain;
   TaLatstress = CaLatstress;
   }

   return res;
}

void FRPConfinedConcrete::Print (OPS_Stream& s, int flag)
{
   s << "  FRPConfinedConcrete: Constitutive (FEM) Model for FRP and Tie - Confined Concrete for Circular Concrete Sections, tag: " << this->getTag() << endln;
   s << "  Compressive Strength of Concrete Core: "  << fpc1 << endln;
   s << "  Compressive Strength of Concrete Cover: " << fpc2 << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  Diameter of the Section: " << D << endln;
   s << "  Concrete Cover: " << c << endln;
   s << "  Elastic Modulus of the Jacket " << Ej << endln;
   s << "  Clear Spacing of FRP Strips (zero if continuous): " << Sj << endln;
   s << "  Thickness of the Jacket: " << tj << endln;
   s << "  Ultimate Strain of the Jacket: " << eju << endln;
   s << "  Spacing of the Stirrups: " << S << endln;
   s << "  Yielding Strength of Longitudinal Steel Bars: " << fyl << endln;
   s << "  Yielding Strength of Stirrups: " << fyh << endln;
   s << "  Diameter of Longitudinal Bars: " << dlong << endln;
   s << "  Diameter of Stirrups " << dtrans << endln;
   s << "  Poisson's Coefficient for Concrete" << v0 << endln;
   s << "  Elastic Modulus for Steel " << Es << endln;
   s << "  Reduction Factor for FRP Ultimate Strain (0.5-0.8) " << k << endln;
   s << "  FRP Jacket Failure Criterion due to Buckling of Longitudinal Compressive Steel Bars (0 = not include it, 1= to include it) " << useBuck << endln;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
FRPConfinedConcrete::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"fc1") == 0)          {// Compressive strength of Concrete Core
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"fc2") == 0)     {// Compressive Strength of Concrete Cover
    return param.addObject(2, this);
  }
  else if (strcmp(argv[0],"epsco") == 0)   {// Strain at compressive strength
    return param.addObject(3, this);
  }
  else if (strcmp(argv[0],"D") == 0)      {// Diameter of the Circular Section
    return param.addObject(4, this);
  }
  else if (strcmp(argv[0],"c") == 0)      {// Concrete Cover
    return param.addObject(5, this);
  }
  else if (strcmp(argv[0],"Ej") == 0)     {// Elastic Modulus for Jacket
    return param.addObject(6, this); 
  }
  else if (strcmp(argv[0],"Sj") == 0)     {// Jacket's Clear Spacing
    return param.addObject(7, this); 
  }
  else if (strcmp(argv[0],"tj") == 0)     {// Thickness of the Jacket
    return param.addObject(8, this);
  }
  else if (strcmp(argv[0],"eju") == 0)    {// Ultimate Strain of the Jacket
    return param.addObject(9, this);
  }
  else if (strcmp(argv[0],"S") == 0)      {// Spacing of the Stirrups
    return param.addObject(10, this);
  }
  else if (strcmp(argv[0],"fyl") == 0)    {// Yielding Strength of Longitudinal Steel Bars
    return param.addObject(11, this);
  }
  else if (strcmp(argv[0],"fyh") == 0)    {// Yielding Strength of Stirrups
    return param.addObject(12, this);
  }
  else if (strcmp(argv[0],"dlong") == 0)  {// Diameter of Longitudinal bar
    return param.addObject(13, this);
  }
  else if (strcmp(argv[0],"dtrans") == 0) {// Diameter of Stirrups
    return param.addObject(14, this);
  }
  else if (strcmp(argv[0],"Es") == 0)     {// Elastic modulus for Steel
    return param.addObject(15, this);
  }
  else if (strcmp(argv[0],"vo") == 0)     {// Poisson's Coefficient for Concrete
    return param.addObject(16, this);
  }
  else if (strcmp(argv[0],"k") == 0)      {// Reduction Factor For FRP Ultimate Strain
    return param.addObject(17, this);
  }
  else if (strcmp(argv[0],"useBuck") == 0)      {// FRP Jacket Failure Criterion due to Buckling of Longitudinal Compressive Steel Bars (0 = not include it, 1= to include it) 
    return param.addObject(18, this);
  }
  
  return -1;
}

int
FRPConfinedConcrete::updateParameter(int parameterID, Information &info)
{
	  
        switch (parameterID) {
        case 1:
                this->fpc1   = info.theDouble;
                break;
	    case 2:
                this->fpc2   = info.theDouble;
                break;
        case 3:
                this->epsc0  = info.theDouble;
                break;
        case 4:
                this->D      = info.theDouble;
                break;
        case 5:
                this->c      = info.theDouble;
                break;
        case 6:
                this->Ej     = info.theDouble;
                break;
		case 7:
                this->Sj     = info.theDouble;
                break;
		case 8:
                this->tj     = info.theDouble;
                break;
		case 9:
                this->eju    = info.theDouble;
                break;
		case 10:
                this->S      = info.theDouble;
                break;
		 case 11:
                this->fyl    = info.theDouble;
                break;
	    case 12:
                this->fyh    = info.theDouble;
                break;
		case 13:
                this->dlong  = info.theDouble;
                break;
		case 14:
                this->dtrans = info.theDouble;
                break;
		case 15:
                this->Es     = info.theDouble;
                break;
		case 16:
                this->v0 = info.theDouble;
                break;
		case 17:
                this->k = info.theDouble;
                break;
		case 18:
                this->useBuck = info.theDouble;
                break;
        default:
                break;
        }
	double Ec = 5700*sqrt(fpc);
        double Ec0 = Ec;
        Ctangent = Ec0;
        CunloadSlope = Ec0;
        Ttangent = Ec0;
        TunloadSlope = CunloadSlope;

        return 0;
}
int
FRPConfinedConcrete::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}

double
FRPConfinedConcrete::getStressSensitivity(int gradNumber, bool conditional)
{      
  // Initialize return value
  double TstressSensitivity = 0.0;
  double dktdh = 0.0;
  double TstrainSensitivity = 0.0;
  
  
  // Pick up sensitivity history variables
  double CminStrainSensitivity = 0.0;
  double CunloadSlopeSensitivity = 0.0;
  double CendStrainSensitivity = 0.0;
  double CstressSensitivity = 0.0;
  double CstrainSensitivity = 0.0;
  if (SHVs != 0) {
    CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
    CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
    CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
    CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
    CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
  }
  
  
  // Assign values to parameter derivatives (depending on what's random)
  double fpcSensitivity = 0.0;
  double epsc0Sensitivity = 0.0;
  
  if (parameterID == 1) {
    fpcSensitivity = 1.0;
  }
  else if (parameterID == 2) {
    epsc0Sensitivity = 1.0;
  }
  
  
  // Strain increment 
  double dStrain = Tstrain - Cstrain;
  
  // Evaluate stress sensitivity 
  if (dStrain < 0.0) {                                    // applying more compression to the material
    
    if (Tstrain < CminStrain) {                     // loading along the backbone curve
      
      if (Tstrain > epsc0) {                  //on the parabola
	
	TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
	  + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
		  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
	
	dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
	  * (1.0-Tstrain/epsc0)
	  - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
	  / (epsc0*epsc0);
      }    
    }
    else if (Tstrain < CendStrain) {        // reloading after an unloading that didn't go all the way to zero stress
      //cerr << "RELOADING AFTER AN UNLOADING THAT DIDN'T GO ALL THE WAY DOWN" << endl;
      TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
	+ CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);
      
      dktdh = CunloadSlopeSensitivity;
    }
    else {
      
      TstressSensitivity = 0.0;
      dktdh = 0.0;
      
    }
  }
  else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
    //cerr << "UNLOADING, BUT NOT ALL THE WAY DOWN" << endl;
    TstressSensitivity = CstressSensitivity 
      + CunloadSlopeSensitivity*dStrain
      + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);
    
    dktdh = CunloadSlopeSensitivity;
  }
  else {                                                                  // unloading all the way down to zero stress
    //cerr << "UNLOADING ALL THE WAY DOWN" << endl;
    
    TstressSensitivity = 0.0;
    dktdh = 0.0;
    
  }
  
  return TstressSensitivity;
}

int
FRPConfinedConcrete::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
        // Initialize unconditaional stress sensitivity
        double TstressSensitivity = 0.0;
        double dktdh = 0.0;


        // Assign values to parameter derivatives (depending on what's random)
        double fpcSensitivity = 0.0;
        double epsc0Sensitivity = 0.0;


        if (parameterID == 1) {
                fpcSensitivity = 1.0;
        }
        else if (parameterID == 2) {
                epsc0Sensitivity = 1.0;
        }

        // Pick up sensitivity history variables
        double CminStrainSensitivity = 0.0;
        double CunloadSlopeSensitivity = 0.0;
        double CendStrainSensitivity = 0.0;
        double CstressSensitivity = 0.0;
        double CstrainSensitivity = 0.0;
        
        if (SHVs == 0) {
                SHVs = new Matrix(5,numGrads);
                CunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
        }
        else {
                CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
                CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
                CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
                CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
                CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
        }


        // Strain increment 
        double dStrain = Tstrain - Cstrain;

        // Evaluate stress sensitivity 
        if (dStrain < 0.0) {                                    // applying more compression to the material

                if (Tstrain < CminStrain) {                     // loading along the backbone curve

                        if (Tstrain > epsc0) {                  //on the parabola
                                
                                TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
                                              + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
                                                  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
                                
                                dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
                                          * (1.0-Tstrain/epsc0)
                                          - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
                                          / (epsc0*epsc0);
                        }
				}
                else if (Tstrain < CendStrain) {        // reloading after an unloading that didn't go all the way to zero stress

                        TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
                                      + CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);

                        dktdh = CunloadSlopeSensitivity;
                }
                else {

                        TstressSensitivity = 0.0;
                        dktdh = 0.0;

                }
        }
        else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
        
                TstressSensitivity = CstressSensitivity 
                                       + CunloadSlopeSensitivity*dStrain
                                           + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);

                dktdh = CunloadSlopeSensitivity;
        }
        else {                                                                  // unloading all the way down to zero stress

                TstressSensitivity = 0.0;
                dktdh = 0.0;

        }

        // Commit some history variables
        (*SHVs)(3,(gradNumber-1)) = TstressSensitivity;
        (*SHVs)(4,(gradNumber-1)) = TstrainSensitivity;


        // Possibly update history variables for the three ordinary history variable derivatives
        double epsTemp, epsTempSensitivity;
        double eta, etaSensitivity;
        double ratio, ratioSensitivity;
        double temp1, temp1Sensitivity;
        double temp2, temp2Sensitivity;
        double TminStrainSensitivity;
        double TunloadSlopeSensitivity;
        double TendStrainSensitivity;

        if (dStrain<0.0 && Tstrain<CminStrain) {

                TminStrainSensitivity = TstrainSensitivity;

                epsTemp = Tstrain;

                epsTempSensitivity = TstrainSensitivity;

                eta = epsTemp/epsc0;

                etaSensitivity = (epsTempSensitivity*epsc0-epsTemp*epsc0Sensitivity) / (epsc0*epsc0);

                if (eta < 2.0) {

                        ratio = 0.145 * eta*eta + 0.13*eta;

                        ratioSensitivity = 0.29 * eta * etaSensitivity + 0.13 * etaSensitivity;

                }
                else {

                        ratio = 0.707*(eta-2.0) + 0.834;

                        ratioSensitivity = 0.707 * etaSensitivity;
                }

                temp1 = Tstrain - ratio * epsc0;

                temp1Sensitivity = TstrainSensitivity - ratioSensitivity * epsc0
                                                          - ratio * epsc0Sensitivity;

                temp2 = Tstress * epsc0 / (2.0*fpc); 
                
                temp2Sensitivity = (2.0*fpc*(TstressSensitivity*epsc0+Tstress*epsc0Sensitivity)
                        -2.0*Tstress*epsc0*fpcSensitivity) / (4.0*fpc*fpc);

                if (temp1 == 0.0) {

		  TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
                }
                else if (temp1 < temp2) {

		  TendStrainSensitivity = TstrainSensitivity - temp1Sensitivity;

		  TunloadSlopeSensitivity = (TstressSensitivity*temp1-Tstress*temp1Sensitivity) / (temp1*temp1);

                }
                else {

		  TendStrainSensitivity = TstrainSensitivity - temp2Sensitivity;

		  TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
                }
        }
        else {
                TminStrainSensitivity = CminStrainSensitivity;
                TunloadSlopeSensitivity = CunloadSlopeSensitivity;
                TendStrainSensitivity = CendStrainSensitivity;
        }



        (*SHVs)(0,(gradNumber-1)) = TminStrainSensitivity;
        (*SHVs)(1,(gradNumber-1)) = TunloadSlopeSensitivity;
        (*SHVs)(2,(gradNumber-1)) = TendStrainSensitivity;

        return 0;
}
int
FRPConfinedConcrete::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

