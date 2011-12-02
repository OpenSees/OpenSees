/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/SmoothConcrete02.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <SmoothConcrete02.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

SmoothConcrete02::SmoothConcrete02
(int tag, double FPC, double EPSC0, double FPCU, double EPSCU, double ETENSION,
 double p_gamma, double p_beta)
 :UniaxialMaterial(tag, MAT_TAG_SmoothConcrete01)
{
	fpc=FPC;
	epsc0=EPSC0;
	fpcu=FPCU;
	epscu=EPSCU;
	etension=ETENSION;
	gamma = p_gamma;
	beta = p_beta;
	epsilon = fabs(epsc0)*1.0e-6;

	// Make all concrete parameters negative
	if (fpc > 0.0) fpc = -fpc;
	if (epsc0 > 0.0)epsc0 = -epsc0;
	if (fpcu > 0.0)fpcu = -fpcu;
	if (epscu > 0.0)epscu = -epscu;
	if (etension == 0.0) etension=0.01*(-epsc0); 

	Cstrain=0.0;	
	Cstress=0.0;	
	CminStrain=0.0;	
	CendStrain=0.0;	
	Celastic=0.0;	
	Ctangent=2.0*fpc/epsc0;
	Cloading=true;	
	Cstate=0;		

	// Set trial values
	this->revertToLastCommit();

//	iopen=0;
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	DTstressDhTstrainFix=0.0;
	DeminDhUpdated=0.0;
	UnconditionalDerivative=0;
// AddingSensitivity:END //////////////////////////////////////

/*	path1=false;
	path2=false;
	path3=false;
	path4=false;
	path5=false;
	path6=false;
	path7=false;
	path8=false;
	path9=false;
	path10=false;
	path11=false;
	path12=false;
	path13=false;
	path14=false;
	path15=false;
	path16=false;
	path17=false;
	path18=false;
	path19=false;
	path20=false;
	path21=false;
	path22=false;
	path23=false;
	path24=false;
	path25=false;
	path26=false;
	path27=false;
	path28=false;
	path29=false;
	path30=false;
*/

}
SmoothConcrete02::SmoothConcrete02()
:UniaxialMaterial(0, MAT_TAG_SmoothConcrete01)
{
	fpc=0.0;
	epsc0=0.0;
	fpcu=0.0;
	epscu=0.0;
	etension=0.0;

	Cstrain=0.0;	
	Cstress=0.0;	
	CminStrain=0.0;	
	CendStrain=0.0;	
	Celastic=0.0;	
	Cloading=true;
	Ctangent=0.0;
	Cstate=0;	

	// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}
SmoothConcrete02::~SmoothConcrete02 ()
{
	if(SHVs!=0){
		delete SHVs;
		SHVs=0;
	}
}

UniaxialMaterial* SmoothConcrete02::getCopy ()
{
   SmoothConcrete02* theCopy = new SmoothConcrete02(this->getTag(),
                    fpc, epsc0, fpcu, epscu, etension, gamma, beta);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CendStrain = CendStrain;
   theCopy->Celastic = Celastic;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;
   theCopy->Cloading = Cloading;
   theCopy->Cstate = Cstate;

   return theCopy;
}

int SmoothConcrete02::revertToStart ()
{
	Cstrain=0.0;	
	Cstress=0.0;	
	CminStrain=0.0;	
	CendStrain=0.0;	
	Celastic=0.0;	
	Cloading=true;	
	Cstate=0;		
	Ctangent=2.0*fpc/epsc0;

	// Reset trial variables and state
	this->revertToLastCommit();
	return 0;
}
int SmoothConcrete02::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   Telastic = Celastic;
   Tloading = Cloading;
   Tstate = Cstate;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int SmoothConcrete02::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(14);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = fpcu;
   data(4) = epscu;
   data(5) = etension;

   // History variables from last converged state
   data(6) = CminStrain;
   data(7) = Celastic;
   data(8) = CendStrain;

   // State variables from last converged state
   data(9) = Cstrain;
   data(10) = Cstress;
   data(11) = Ctangent;
   if(Cloading) data(12) = 1.0;
   else data(12) = 0.0;
   data(13) = (double)Cstate;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "SmoothConcrete01::sendSelf() - failed to send data\n";

   return res;
}

int SmoothConcrete02::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
	int res = 0;
	static Vector data(14);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);

	if (res < 0) {
		opserr << "SmoothConcrete01::recvSelf() - failed to receive data\n";
		this->setTag(0);      
	}
	else {
		this->setTag(int(data(0)));

		// Material properties
		fpc = data(1);
		epsc0 = data(2);
		fpcu = data(3);
		epscu = data(4);
		etension = etension;

		// History variables from last converged state
		CminStrain = data(6);
		Celastic = data(7);
		CendStrain = data(8);

		// State variables from last converged state
		Cstrain = data(9);
		Cstress = data(10);
		Ctangent = data(11);
		if(data(12)==1.0) Cloading=true;
		else Cloading=false;
		Cstate=(int)data(13);
   }

   return res;
}

void SmoothConcrete02::Print (OPS_Stream& s, int flag)
{
   s << "SmoothConcrete01, tag: " << this->getTag() << endln;
   s << "  fpc: " << fpc << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  fpcu: " << fpcu << endln;
   s << "  epscu: " << epscu << endln;
}
void SmoothConcrete02::unload 
(const double& emin, const double& sigmin, double& eend, double& elastic)
{
	double etemp,ratio,de;
	if(emin<=epscu){
		etemp=epscu;
	}else{
		etemp=emin;
	}
	double gzai=etemp/epsc0;
	if(gzai>=2){
		ratio=0.707*(gzai-2.0)+0.84;
	}else {
		ratio=gzai*(0.145*gzai+0.13);
	}
	double de1=emin-ratio*epsc0;
	double de2=sigmin*epsc0/(2.0*fpc);
	if(de1<=de2) {
		de=de1;
	}else{
		de=de2;
	}
	eend=emin-de;
	elastic=sigmin/de;
}
void SmoothConcrete02::unloadwithSens1 
(const double& emin, const double& sigmin, const double& DsigminDemin,
 double& eend, double& elastic, double& DeendDemin, double& DelasticDemin)
{
	double etemp,gzai,DetempDemin,DgzaiDemin,ratio;
	double DratioDemin,de1,de2,de,DeDemin;
	if(emin<=epscu) {
		etemp=epscu;
		DetempDemin=0.0;
	}else {
		etemp=emin;
		DetempDemin=1.0;
	}
	gzai=etemp/epsc0;
	DgzaiDemin=DetempDemin/epsc0;
	if(gzai>=2.0) {
		ratio=0.707*(gzai-2.0)+0.84;
		DratioDemin=0.707*DgzaiDemin;
	}else{
		ratio=gzai*(0.145*gzai+0.13);
		DratioDemin=(0.29*gzai+0.13)*DgzaiDemin;
	}
	de1=emin-ratio*epsc0;
	de2=sigmin*epsc0/(2.0*fpc);
	if(de1<=de2){
		de=de1;
		DeDemin=1.0-DratioDemin*epsc0;
	}else{
		de=de2;
		DeDemin=DsigminDemin*epsc0/(2.0*fpc);
	}
	eend=emin-de;
	DeendDemin=1.0-DeDemin;
	elastic=sigmin/de;
	DelasticDemin=DsigminDemin/de-sigmin/(de*de)*DeDemin;
}

void SmoothConcrete02::unloadwithSens2 
(const double& emin, const double& sigmin, 
 const double& DsigminDemin, const double& DsigminDh,
 const double& DfpcDh, const double& DfpcuDh, 
 const double& Depsc0Dh, const double& DepscuDh,
 double& eend, double& elastic, 
 double& DeendDemin, double& DelasticDemin,
 double& DeendDh, double& DelasticDh)
{
	double etemp,gzai,DetempDemin,DgzaiDemin,ratio;
	double DetempDh,DgzaiDh,DratioDh;
	double DratioDemin,de1,de2,de,DeDh,DeDemin;
	if(emin<=epscu) {
		etemp=epscu;
		DetempDemin=0.0;
		DetempDh=DepscuDh;
	}else {
		etemp=emin;
		DetempDemin=1.0;
		DetempDh=0.0;
	}
	gzai=etemp/epsc0;
	DgzaiDemin=DetempDemin/epsc0;
	DgzaiDh=DetempDh/epsc0-etemp/(epsc0*epsc0)*Depsc0Dh;
	if(gzai>=2.0) {
		ratio=0.707*(gzai-2.0)+0.84;
		DratioDemin=0.707*DgzaiDemin;
		DratioDh=0.707*DgzaiDh;
	}else{
		ratio=0.145*gzai*gzai+0.13*gzai;
		DratioDemin=(0.29*gzai+0.13)*DgzaiDemin;
		DratioDh=(0.29*gzai+0.13)*DgzaiDh;
	}
	de1=emin-ratio*epsc0;
	de2=sigmin*epsc0/(2.0*fpc);
	if(de1<=de2){
		de=de1;
		DeDemin=1.0-DratioDemin*epsc0;
		DeDh=-DratioDh*epsc0-ratio*Depsc0Dh;
	}else{
		de=de2;
		DeDemin=DsigminDemin*epsc0/(2.0*fpc);
		DeDh=DsigminDh*epsc0/(2.0*fpc)
			+sigmin/(2.0*fpc)*(Depsc0Dh-epsc0/fpc*DfpcDh);
	}
	eend=emin-de;
	DeendDemin=1.0-DeDemin;
	DeendDh=-DeDh;
	elastic=sigmin/de;
	DelasticDemin=DsigminDemin/de-sigmin/(de*de)*DeDemin;
	DelasticDh=DsigminDh/de-sigmin/(de*de)*DeDh;
}

void SmoothConcrete02::solemin2
(const double& sigma, const double& eps, 
 double& emin, double& eend, double& elastic) 
{
// initial value
	double temin;
	double tsigmin,tanmin,DtanminDemin;
	double teend,telastic,DeendDemin,DelasticDemin;

	int iend=0;
	temin=emin;
	while(iend==0){
		envelopewithSens1(temin,tsigmin,tanmin,DtanminDemin); 
		unloadwithSens1(temin,tsigmin,tanmin, teend, telastic, 
			            DeendDemin, DelasticDemin);
		double fff=sigma-telastic*(eps-teend);
		if( fabs(fff/sigma)<1.0e-6) break;
		double DfDemin=-DelasticDemin*(eps-teend)+telastic*DeendDemin;
		temin-=fff/DfDemin;
	}
	emin=temin;
	eend=teend;
	elastic=telastic;
}


void SmoothConcrete02::solemin2Sens
(const double& emin, const double& eps, 
 const double& DsigmaDh, const double& DepsDh,
 const double& DfpcDh, const double& DfpcuDh, 
 const double& Depsc0Dh, const double& DepscuDh, 
 double& DeminDh, double& DeendDh, double& DelasticDh)
{
// initial value
	double sigmin,DsigminDemin,DtanminDemin;
	double DsigminDheminFix,DtanminDheminFix;
	double DeendDemin,DelasticDemin; 
	double DeendDheminFix,DelasticDheminFix; 
	double eend,elastic,rhs,lhs;

	envelopewithSens1(emin,sigmin,DsigminDemin,DtanminDemin); 
	envelopeSens(emin, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				 DsigminDheminFix, DtanminDheminFix);
	unloadwithSens2(emin, sigmin, DsigminDemin, DsigminDheminFix,
					DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					eend, elastic, 
					DeendDemin, DelasticDemin, 
					DeendDheminFix, DelasticDheminFix);
	
	rhs=DsigmaDh-DelasticDheminFix*(eps-eend)
				-elastic*(DepsDh-DeendDheminFix);
	lhs=DelasticDemin*(eps-eend)-elastic*DeendDemin;
	DeminDh=rhs/lhs;
	DeendDh=DeendDheminFix+DeendDemin*DeminDh;
	DelasticDh=DelasticDheminFix+DelasticDemin*DeminDh;
}
void SmoothConcrete02::solemin1
(const double& sigma, const double& eps, 
 double& emin, double& eend, double& elastic, 
 double& amin, double& bmin, double& cmin, double& dmin)
{
	double temin,tsigmin,tanmin,DtanminDemin;
	double teend,telastic,DeendDemin,DelasticDemin;
	double a,b,c,d;
	double deminend,DdeminendDemin,emin2,emin3;
	double Demin2Demin,Demin3Demin,sigmin2,Dsigmin2Demin;
	double sigmin3,tanmin3,Dtanmin3Demin3,Dsigmin3Demin,Dtanmin3Demin;
	double DaDemin,DbDemin,DcDemin,DdDemin,fff,DfDemin;

	int iend=0;
	temin=emin;
	while(iend==0){
		envelopewithSens1(temin,tsigmin,tanmin,DtanminDemin); 
		unloadwithSens1(temin,tsigmin,tanmin,
			            teend,telastic,DeendDemin,DelasticDemin);
		deminend=temin-teend;
		DdeminendDemin=1.0-DeendDemin;
		emin2=temin-gamma*deminend;
		emin3=temin+beta*deminend;
		Demin2Demin=1.0-gamma*DdeminendDemin;
		Demin3Demin=1.0+beta*DdeminendDemin;
		sigmin2=(1.0-gamma)*telastic*deminend;
		Dsigmin2Demin=(1.0-gamma)*(DelasticDemin*deminend+telastic*DdeminendDemin);
		envelopewithSens1(emin3,sigmin3,tanmin3,Dtanmin3Demin3); 
		Dsigmin3Demin=tanmin3*Demin3Demin;
		Dtanmin3Demin=Dtanmin3Demin3*Demin3Demin;
		polinomial(emin2,sigmin2,telastic,emin3,sigmin3,tanmin3,a, b, c, d); 
		polinomialwithSens(emin2,sigmin2,telastic,emin3,sigmin3,tanmin3,a,b,c,d,
			               Demin2Demin,Dsigmin2Demin,DelasticDemin,
						   Demin3Demin,Dsigmin3Demin,Dtanmin3Demin,
						   DaDemin,DbDemin,DcDemin,DdDemin);
		fff=sigma-a*eps*eps*eps-b*eps*eps-c*eps-d;
		if( fabs(fff/sigma)<1.0e-6) break;
		DfDemin=-DaDemin*eps*eps*eps-DbDemin*eps*eps-DcDemin*eps-DdDemin;
		temin-=fff/DfDemin;
	}
	emin=temin;
	eend=teend;
	elastic=telastic;
	amin=a;
	bmin=b;
	cmin=c;
	dmin=d;
}

void SmoothConcrete02::solemin1Sens
(const double& emin, const double eps, const double& DsigmaDh, const double& DepsDh,
 const double& DfpcDh, const double& DfpcuDh, 
 const double& Depsc0Dh, const double& DepscuDh, 
 double& DeminDh, double& DaDh, double& DbDh, double& DcDh, double& DdDh) 
{ 
	double sigmin,DsigminDemin,DtanminDemin;
	double DsigminDheminFix,DtanminDheminFix;
	double eend,elastic,DeendDemin,DelasticDemin;
	double DeendDheminFix,DelasticDheminFix;
	double deminend,DdeminendDemin,DdeminendDheminFix; 
	double emin2,Demin2Demin,Demin2DheminFix;
	double emin3,Demin3Demin,Demin3DheminFix;
	double sigmin2,Dsigmin2Demin,Dsigmin2DheminFix;
	double sigmin3,Dsigmin3Demin3,Dtanmin3Demin3;
	double Dsigmin3Dhemin3Fix,Dtanmin3Dhemin3Fix;
	double Dsigmin3Demin,Dsigmin3DheminFix,Dtanmin3Demin,Dtanmin3DheminFix;
	double a,b,c,d,DaDemin,DbDemin,DcDemin,DdDemin;
	double DaDheminFix,DbDheminFix,DcDheminFix,DdDheminFix;
	double rhs,lhs,eps3,eps2;

	/// emin: solution obtained already //
	envelopewithSens1(emin, sigmin, DsigminDemin, DtanminDemin);
	envelopeSens(emin, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				 DsigminDheminFix, DtanminDheminFix);
	unloadwithSens2(emin, sigmin, DsigminDemin, DsigminDheminFix,
					DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,					
					eend, elastic,DeendDemin, DelasticDemin, 
					DeendDheminFix, DelasticDheminFix);
	deminend=emin-eend;
	DdeminendDemin=1.0-DeendDemin;
	DdeminendDheminFix=-DeendDheminFix;
	emin2=emin-gamma*deminend;
	Demin2Demin=1.0-gamma*DdeminendDemin;
	Demin2DheminFix=-gamma*DdeminendDheminFix;
	emin3=emin+beta*deminend;
	Demin3Demin=1.0+beta*DdeminendDemin;
	Demin3DheminFix=beta*DdeminendDheminFix;
	sigmin2=(1.0-gamma)*elastic*deminend;
	Dsigmin2Demin=(1.0-gamma)*(DelasticDemin*deminend+elastic*DdeminendDemin);
	Dsigmin2DheminFix=(1.0-gamma)*
		(DelasticDheminFix*deminend+elastic*DdeminendDheminFix);

	envelopewithSens1(emin3, sigmin3, Dsigmin3Demin3, Dtanmin3Demin3);
	envelopeSens(emin3, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				 Dsigmin3Dhemin3Fix, Dtanmin3Dhemin3Fix);
	Dsigmin3Demin=Dsigmin3Demin3*Demin3Demin;
	Dsigmin3DheminFix=Dsigmin3Dhemin3Fix+Dsigmin3Demin3*Demin3DheminFix;
	Dtanmin3Demin=Dtanmin3Demin3*Demin3Demin;
	Dtanmin3DheminFix=Dtanmin3Dhemin3Fix+Dtanmin3Demin3*Demin3DheminFix;

	polinomial(emin2,sigmin2,elastic,emin3,sigmin3,Dsigmin3Demin3,
		       a, b, c, d); 
	polinomialwithSens
		(emin2,sigmin2,elastic,emin3,sigmin3,Dsigmin3Demin3, a, b, c, d,
		 Demin2Demin,Dsigmin2Demin,DelasticDemin,
		 Demin3Demin,Dsigmin3Demin,Dtanmin3Demin,
		 DaDemin, DbDemin, DcDemin, DdDemin);
	polinomialwithSens
		(emin2,sigmin2,elastic,emin3,sigmin3,Dsigmin3Demin3, a, b, c, d,
		 Demin2DheminFix,Dsigmin2DheminFix,DelasticDheminFix,
		 Demin3DheminFix,Dsigmin3DheminFix,Dtanmin3DheminFix,
		 DaDheminFix, DbDheminFix, DcDheminFix, DdDheminFix);

	eps2=eps*eps;
	eps3=eps2*eps;
	rhs=DsigmaDh-DaDheminFix*eps3-DbDheminFix*eps2-DcDheminFix*eps-DdDheminFix
				-(3.0*a*eps2+2.0*b*eps+c)*DepsDh;

	lhs=DaDemin*eps3+DbDemin*eps2+DcDemin*eps+DdDemin;

	DeminDh=rhs/lhs;
	DaDh=DaDheminFix+DaDemin*DeminDh;
	DbDh=DbDheminFix+DbDemin*DeminDh;
	DcDh=DcDheminFix+DcDemin*DeminDh;
	DdDh=DdDheminFix+DdDemin*DeminDh;
}

int SmoothConcrete02::setTrialStrain (double strain, double strainRate)
{
	bool upTstress=false;
	bool upTminStrain=false;
	bool upTendStrain=false;
	bool upTelastic=false;
	bool upTtangent=false;
	bool upTstate=false;
	// Set trial strain
	Tstrain = strain;

//	if(Tstrain <= -2.00e-3){ 
//		double abc=0.0;
//	}
	double delta_e =Tstrain - Cstrain;

	if(fabs(delta_e)<DBL_EPSILON){
		/// no change because incremental strain is too small ///
		this->revertToLastCommit();
		return 0;
	}	

	if(delta_e>0.0) {
		Tloading=false;
	}else{ 
		Tloading=true;
	}

	if(Cstate==0){
		if(Tloading){
			envelope(Tstrain, Tstress, Ttangent);	
			Tstate=1;
			TminStrain=Tstrain;  
			unload(TminStrain, Tstress, TendStrain, Telastic);
			upTstate=true;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}else{
			// Trial unloading //
			Tstress=0.0;
			Ttangent=0.0;
			upTstress=true;
			upTtangent=true;
			if(Tstrain<etension){
				Tstate=-1;
				upTstate=true;
				TminStrain=-Tstrain;
				TendStrain=Tstrain;
				Telastic=0.0;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}else{
				Tstate=-2;
				upTstate=true;
				TminStrain=-etension;
				TendStrain=etension;
				Telastic=0.0;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}
		}
	}else if(Cstate==1){
		if(Tloading){
			envelope(Tstrain, Tstress, Ttangent);
			Tstate=1;
			TminStrain=Tstrain;
			unload(TminStrain, Tstress, TendStrain, Telastic);
			upTstate=true;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}else{
			TminStrain=CminStrain;
			TendStrain=CendStrain;
			Telastic=Celastic;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
			double deminend=TminStrain-TendStrain;
			double emin2=TminStrain-gamma*deminend;
			double eend2=TendStrain+gamma*deminend;
			double eend3=TendStrain-beta*deminend;
			if(Tstrain<=eend2){
				Tstress=Telastic*(Tstrain-TendStrain);
				Ttangent=Telastic;
				if(Tstrain<emin2) Tstate=2;
				else Tstate=3;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}else if(Tstrain<eend3){
				double sigend2=gamma*Telastic*deminend;
				double a,b,c,d;
				polinomial(eend2, sigend2, Telastic, 
						   eend3, 0.0, 0.0, a,b,c,d);
				double Tstrain2=Tstrain*Tstrain;
				double Tstrain3=Tstrain2*Tstrain;
				Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
				Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
				Tstate=4;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}else{
				Tstress=0.0;
				Ttangent=0.0;
				Tstate=5;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}
		}
	}else if(Cstate==3||Cstate==4||Cstate==5){
		double deminend=CminStrain-CendStrain;
		double emin3=CminStrain+beta*deminend;
		if(Tstrain<=emin3){
			envelope(Tstrain, Tstress, Ttangent);
			Tstate=1;
			TminStrain=Tstrain;
			unload(TminStrain, Tstress, TendStrain, Telastic );
			upTstate=true;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}else{
			double eend2=CendStrain+gamma*deminend;
			double eend3=CendStrain-beta*deminend;
			double emin2=CminStrain-gamma*deminend;
			if(Tstrain<emin2){
				double sigmin2=(1.0-gamma)*Celastic*deminend;
				double sigmin3,elast3,a,b,c,d;
				envelope(emin3, sigmin3, elast3);
				polinomial(emin2, sigmin2, Celastic, 
						   emin3, sigmin3, elast3, a,b,c,d);
				double Tstrain2=Tstrain*Tstrain;
				double Tstrain3=Tstrain2*Tstrain;
				Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
				Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
				Tstate=2;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}else if(Tstrain<=eend2){
				Tstress=Celastic*(Tstrain-CendStrain);
				Ttangent=Celastic;
				Tstate=3;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}else if(Tstrain<eend3){
				double sigend2=gamma*Celastic*deminend;
				double a,b,c,d;
				polinomial(eend2, sigend2, Celastic, 
						   eend3, 0.0,0.0, a,b,c,d);
				double Tstrain2=Tstrain*Tstrain;
				double Tstrain3=Tstrain2*Tstrain;
				Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
				Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
				Tstate=4;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}else{
				Tstress=0.0;
				Ttangent=0.0;
				Tstate=5;
				upTstate=true;
				upTstress=true;
				upTtangent=true;
			}
			TminStrain=CminStrain;
			TendStrain=CendStrain;
			Telastic=Celastic;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}
	}else if(Cstate==2){
		if(Tloading){
			if(Cloading){
				double deminend=CminStrain-CendStrain;
				double emin3=CminStrain+beta*deminend;
				if(Tstrain<=emin3){
					envelope(Tstrain, Tstress, Ttangent);
					Tstate=1;
					TminStrain=Tstrain;
					unload(TminStrain, Tstress, TendStrain, Telastic );
					upTstate=true;
					upTstress=true;
					upTtangent=true;
					upTminStrain=true;
					upTendStrain=true;
					upTelastic=true;
				}else{
					double emin2=CminStrain-gamma*deminend;
					double sigmin3,elast3,sigmin2,a,b,c,d;
					envelope(emin3, sigmin3, elast3);
					sigmin2=(1.0-gamma)*Celastic*deminend;
					polinomial(emin2, sigmin2, Celastic, 
						       emin3, sigmin3,elast3, a,b,c,d);
					double Tstrain2=Tstrain*Tstrain;
					double Tstrain3=Tstrain2*Tstrain;
					Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
					Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
					Tstate=2;
					TminStrain=CminStrain;
					TendStrain=CendStrain;
					Telastic=Celastic;
					upTstate=true;
					upTstress=true;
					upTtangent=true;
					upTminStrain=true;
					upTendStrain=true;
					upTelastic=true;
				}
			}else{
				double a,b,c,d;
				TminStrain=CminStrain; 
				solemin1(Cstress, Cstrain, TminStrain, TendStrain, Telastic,
			             a,b,c,d);		 
				double deminend=TminStrain-TendStrain;
				double emin3=TminStrain+beta*deminend;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				if(Tstrain<=emin3){
					envelope(Tstrain, Tstress, Ttangent);
					Tstate=1;
					TminStrain=Tstrain;
					unload(TminStrain, Tstress, TendStrain, Telastic);
					upTstate=true;
					upTstress=true;
					upTtangent=true;
					upTminStrain=true;
					upTendStrain=true;
					upTelastic=true;
				}else{
					double Tstrain2=Tstrain*Tstrain;
					double Tstrain3=Tstrain2*Tstrain;
					Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
					Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
					Tstate=2;
					upTstate=true;
					upTstress=true;
					upTtangent=true;
				}
			}
		}else{
			if(Cloading){
				TminStrain=CminStrain; 
				solemin2(Cstress, Cstrain,TminStrain,TendStrain, Telastic); 
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				double deminend=TminStrain-TendStrain;
				double emin2=TminStrain-gamma*deminend;
				double eend2=TendStrain+gamma*deminend;
				double eend3=TendStrain-beta*deminend;
				if(Tstrain<=eend2){
					Tstress=Telastic*(Tstrain-TendStrain);
					Ttangent=Telastic;
					upTstress=true;
					upTtangent=true;
					if(Tstrain>=emin2){
						Tstate=3;
					}else{
						Tstate=2;
					}
					upTstate=true;
				}else if(Tstrain<eend3){
					double sigend2=gamma*Telastic*deminend;
					double a,b,c,d;
					polinomial(eend2, sigend2, Celastic, 
							   eend3, 0.0,0.0, a,b,c,d);
					double Tstrain2=Tstrain*Tstrain;
					double Tstrain3=Tstrain2*Tstrain;
					Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
					Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
					Tstate=4;
					upTstress=true;
					upTtangent=true;
					upTstate=true;
				}else{
					Tstress=0.0;
					Ttangent=0.0;
					Tstate=5;
					upTstress=true;
					upTtangent=true;
					upTstate=true;
				}
			}else{
				TminStrain=CminStrain; 
				TendStrain=CendStrain; 
				Telastic=Celastic; 
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				double deminend=TminStrain-TendStrain;
				double emin2=TminStrain-gamma*deminend;
				double eend2=TendStrain+gamma*deminend;
				double eend3=TendStrain-beta*deminend;
				if(Tstrain<=eend2){
					Tstress=Telastic*(Tstrain-TendStrain);
					Ttangent=Telastic;
					upTstress=true;
					upTtangent=true;
					if(Tstrain>=emin2){
						Tstate=3;
					}else{
						Tstate=2;
					}
					upTstate=true;
				}else if(Tstrain<eend3){
					double sigend2=gamma*Telastic*deminend;
					double a,b,c,d;
					polinomial(eend2, sigend2, Celastic, 
							   eend3, 0.0,0.0, a,b,c,d);
					double Tstrain2=Tstrain*Tstrain;
					double Tstrain3=Tstrain2*Tstrain;
					Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
					Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
					Tstate=4;
					upTstress=true;
					upTtangent=true;
					upTstate=true;
				}else{
					Tstress=0.0;
					Ttangent=0.0;
					Tstate=5;
					upTstress=true;
					upTtangent=true;
					upTstate=true;
				}
			}
		}
	}else if(Cstate==-1){
		if(Tloading){
			if(Tstrain<=CminStrain){
				envelope(Tstrain, Tstress, Ttangent);
				Tstate=1;
				TminStrain=Tstrain;
				unload(TminStrain, Tstress, TendStrain, Telastic);
				upTstate=true;
				upTstress=true;
				upTtangent=true;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}else{
				TminStrain=CminStrain;
				TendStrain=CendStrain;
				Telastic=0.0;
				Tstate=-3;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				upTstate=true;
				double sigb,tanb,a,b,c,d;
				envelope(CminStrain, sigb, tanb);
				polinomial(CendStrain, 0.0, 0.0, 
						   CminStrain, sigb, tanb, a,b,c,d);
				double Tstrain2=Tstrain*Tstrain;
				double Tstrain3=Tstrain2*Tstrain;
				Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
				Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
				upTstress=true;
				upTtangent=true;
			}
		}else{
			Tstress=0.0;
			Ttangent=0.0;
			upTstress=true;
			upTtangent=true;
			if(Tstrain<etension){
				Tstate=-1;
				TminStrain=-Tstrain;
				TendStrain=Tstrain;
				Telastic=0.0;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				upTstate=true;
			}else{
				Tstate=-2;
				TminStrain=-etension;
				TendStrain=etension;
				Telastic=0.0;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
				upTstate=true;
			}
		}
	}else if(Cstate==-2){
		if(Tloading){
			if(Tstrain<=CminStrain){
				envelope(Tstrain, Tstress, Ttangent);
				Tstate=1;
				TminStrain=Tstrain;
				unload(TminStrain, Tstress, TendStrain, Telastic);
				upTstate=true;
				upTstress=true;
				upTtangent=true;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}else{
				if(Tstrain<=CendStrain){
					Tstate=-3;
					upTstate=true;
					double sigb,tanb,a,b,c,d;
					envelope(CminStrain, sigb, tanb);
					polinomial(CendStrain, 0.0, 0.0, 
							  CminStrain, sigb, tanb, a,b,c,d);
					double Tstrain2=Tstrain*Tstrain;
					double Tstrain3=Tstrain2*Tstrain;
					Tstress=a*Tstrain3+b*Tstrain2+c*Tstrain+d;
					Ttangent=3.0*a*Tstrain2+2.0*b*Tstrain+c;
					upTstress=true;
					upTtangent=true;
				}else{
					Tstate=-2;
					Tstress=0.0;
					Ttangent=0.0;
					upTstate=true;
					upTstress=true;
					upTtangent=true;
				}
				TminStrain=CminStrain;
				TendStrain=CendStrain;
				Telastic=Celastic;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}
		}else{
			Tstate=-2;
			Tstress=0.0;
			Ttangent=0.0;
			Telastic=0.0;
			TminStrain=CminStrain;
			TendStrain=CendStrain;
			upTstate=true;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}
	}else if(Cstate==-3){
		if(Tstrain<=CminStrain){
			envelope(Tstrain, Tstress, Ttangent);
			Tstate=1;
			TminStrain=Tstrain;
			unload(TminStrain, Tstress, TendStrain, Telastic);
			upTstate=true;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}else if(Tstrain<CendStrain){
			Tstate=-3;
			upTstate=true;
			double sigb,tanb,a,b,c,d;
			envelope(CminStrain, sigb, tanb);
			polinomial(CendStrain, 0.0, 0.0, 
  					   CminStrain, sigb, tanb, a,b,c,d);
			Tstress=a*Tstrain*Tstrain*Tstrain+b*Tstrain*Tstrain+c*Tstrain+d;
			Ttangent=3.0*a*Tstrain*Tstrain+3.0*b*Tstrain*Tstrain+c;
			TendStrain=CendStrain;
			TminStrain=CminStrain;
			Telastic=Celastic;
			upTstress=true;
			upTtangent=true;
			upTminStrain=true;
			upTendStrain=true;
			upTelastic=true;
		}else{
			Tstress=0.0;
			Ttangent=0.0;
			upTstress=true;
			upTtangent=true;
			if(Tstrain<etension){
				Tstate=-1;
				TminStrain=-Tstrain;
				TendStrain=Tstrain;
				Telastic=0.0;
				upTstate=true;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}else{
				Tstate=-1;
				TminStrain=-etension;
				TendStrain=etension;
				Telastic=0.0;
				upTstate=true;
				upTminStrain=true;
				upTendStrain=true;
				upTelastic=true;
			}
		}
	}
	
//	if(!upTstate){opserr<<" upTstate is false\n"; exit(-1); }
//	if(!upTminStrain){opserr<<" upTminStrain is false\n"; exit(-1); }
//	if(!upTendStrain){opserr<<" upTendStrain is false\n"; exit(-1); }
//	if(!upTelastic){opserr<<" upTelastic is false\n"; exit(-1); }
//	if(!upTtangent){opserr<<" upTtangent is false\n"; exit(-1); }
//	if(!upTstress){opserr<<" upTstress is false\n"; exit(-1); }

	return 0;
}

int SmoothConcrete02::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CendStrain = TendStrain;
   Celastic = Telastic;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;
   Cloading = Tloading;
   Cstate = Tstate;

//	output2.setf(ios::right);
//	output2.setf(ios::scientific, ios::floatfield);
//	output2<<setw(15)<<setprecision(5)<<CminStrain;
//	output2<<setw(15)<<setprecision(5)<<Cstress;
//	output2<<setw(15)<<setprecision(5)<<Cstrain;
//	output2<<"\n";
//	output2.flush();
   return 0;
}

int
SmoothConcrete02::setParameter(const char **argv, int argc, Parameter &param)
{
	if (strcmp(argv[0],"fc") == 0) {// Compressive strength
		return param.addObject(1, this);
	}
	if (strcmp(argv[0],"epsco") == 0) {// Strain at compressive strength
		return param.addObject(2, this);
	}
	if (strcmp(argv[0],"fcu") == 0) {// Crushing strength
		return param.addObject(3, this);
	}
	if (strcmp(argv[0],"epscu") == 0) {// Strain at crushing strength
		return param.addObject(4, this);
	}
	else {
		opserr << "WARNING: Could not set parameter in SmoothConcrete01! " << endln;
		return -1;
	}
}
int
SmoothConcrete02::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->fpc = info.theDouble;
		break;
	case 2:
		this->epsc0 = info.theDouble;
		break;
	case 3:
		this->fpcu = info.theDouble;
		break;
	case 4:
		this->epscu = info.theDouble;
		break;
	default:
		break;
	}
        
	// Make all concrete parameters negative
	if (fpc > 0.0)
		fpc = -fpc;

	if (epsc0 > 0.0)
		epsc0 = -epsc0;

	if (fpcu > 0.0)
		fpcu = -fpcu;

	if (epscu > 0.0)
		epscu = -epscu;

	// Initial tangent
	Ctangent = 2*fpc/epsc0;
	Ttangent = 2*fpc/epsc0;
	return 0;
}


int
SmoothConcrete02::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return parameterID;
}



double
SmoothConcrete02::getStrainSensitivity(int gradNumber)
{
	return 0.0; 
}

double
SmoothConcrete02::getTangentSensitivity(int gradNumber)
{
	return 0.0; 
}

double
SmoothConcrete02::getDampTangentSensitivity(int gradNumber)
{
	return 0.0; 
}

double
SmoothConcrete02::getRhoSensitivity(int gradNumber)
{
	return 0.0; 
}

void SmoothConcrete02::polinomial 
(const double& eps1, const double& sig1, const double& tan1, 
 const double& eps2, const double& sig2, const double& tan2,
 double& a, double& b, double& c, double& d) 
{
/// COMPUTE COEFFICIENTS OF THE 3RD ORDER POLYNOMIAL ///
	double de=eps1-eps2;
	a=-2.0*(sig1-sig2)/(de*de*de)+(tan1+tan2)/(de*de);
	b=(tan1-tan2)/(2.0*de)-3.0*a*(eps1+eps2)/2.0;
	c=3.0*a*eps1*eps2-(tan1*eps2-tan2*eps1)/de;
	d=sig1-eps1*(a*eps1*eps1+b*eps1+c);
}
void SmoothConcrete02::polinomialwithSens
(const double& eps1, const double& sig1, const double& tan1, 
 const double& eps2, const double& sig2, const double& tan2,
 const double& a, const double& b, const double& c, const double& d, 
 const double& Deps1Dh, const double& Dsig1Dh, const double& Dtan1Dh,
 const double& Deps2Dh, const double& Dsig2Dh, const double& Dtan2Dh,
 double& DaDh, double& DbDh, double& DcDh, double& DdDh)
{
	double de=eps1-eps2;
	double DeDh=Deps1Dh-Deps2Dh;
	double de2=de*de;
	double de3=de2*de;
	double de4=de3*de;

	DaDh=-2.0*(Dsig1Dh-Dsig2Dh)/de3+6.0*(sig1-sig2)/de4*DeDh
		+(Dtan1Dh+Dtan2Dh)/de2-2.0*(tan1+tan2)/de3*DeDh;

	DbDh=(Dtan1Dh-Dtan2Dh)/(2.0*de)-(tan1-tan2)/(2.0*de2)*DeDh
		-3.0*(eps1+eps2)*DaDh/2.0-3.0*a/2.0*(Deps1Dh+Deps2Dh);

	DcDh=3.0*DaDh*eps1*eps2+3.0*a*Deps1Dh*eps2+3.0*a*eps1*Deps2Dh
		-(eps2*Dtan1Dh+tan1*Deps2Dh-tan2*Deps1Dh-eps1*Dtan2Dh)/de
		+(tan1*eps2-tan2*eps1)/de2*DeDh;

	DdDh=Dsig1Dh-eps1*(DaDh*eps1*eps1+DbDh*eps1+DcDh)
		-(3.0*a*eps1*eps1+2.0*b*eps1+c)*Deps1Dh;
}

double
SmoothConcrete02::getStressSensitivity(int gradNumber, bool conditional)
{
//	if(iopen==0){
//	conc02_output2.open("SmoothConcrete02PATH.txt", ios::out);
//	conc02_output2 << "\n";
//	conc02_output2.flush();
//	iopen=1;
//	}


	UnconditionalDerivative=gradNumber;

	double DeminDh,DstressDh,DstrainDh;

    double 	DfpcDh=0.0;
    double 	DfpcuDh=0.0;
    double 	Depsc0Dh=0.0;
    double 	DepscuDh=0.0;	
	switch (parameterID) {
		case 1: DfpcDh=1.0;	break;
		case 2: Depsc0Dh=1.0; break;
		case 3: DfpcuDh=1.0; break;
		case 4: DepscuDh=1.0; break;
		default: break;
	}

	if (SHVs != 0) {
		DeminDh   = (*SHVs)(0,(gradNumber-1));
		DstressDh = (*SHVs)(1,(gradNumber-1));
		DstrainDh = (*SHVs)(2,(gradNumber-1));
	}else{
		DeminDh = 0.0;
		DstressDh = 0.0;
		DstrainDh = 0.0;
	}

	int path=0;
	double DsigDh=0.0;
	double DtanDh=0.0;

	if(Cstate==0){
		if(Tloading){
			envelopeSens 
			(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 DsigDh, DtanDh);
			path=1;
		}else{
			DsigDh=0.0;
			path=2;
		}
	}else if(Cstate==1){
		if(Tloading){
			path=3;
			envelopeSens 
			(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 DsigDh, DtanDh);
		}else{
			double sigmin,DsigminDemin,DtanminDemin;
			double DsigminDheminFix,DtanminDheminFix;
			double eend,elastic,DeendDemin,DelasticDemin;
			double DeendDheminFix,DelasticDheminFix;
			double DeendDh,DelasticDh;
			envelopewithSens1 
			(CminStrain, sigmin, DsigminDemin, DtanminDemin);
			envelopeSens 
			(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 DsigminDheminFix, DtanminDheminFix);
			unloadwithSens2 
			(CminStrain, sigmin, DsigminDemin, DsigminDheminFix,
			 DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 eend, elastic, DeendDemin, DelasticDemin,
			 DeendDheminFix, DelasticDheminFix);
			 DeendDh=DeendDheminFix+DeendDemin*DeminDh;
			 DelasticDh=DelasticDheminFix+DelasticDemin*DeminDh; 
			if(Tstate==2||Tstate==3){
				path=4;
				DsigDh=DelasticDh*(Tstrain-CendStrain)-elastic*DeendDh;
			}else if(Tstate==4){
				path=5;
				double deminend = CminStrain-eend;
				double DdeminendDh = DeminDh-DeendDh;
				double eend2=eend+gamma*deminend;
				double Deend2Dh=DeendDh+gamma*DdeminendDh;
				double eend3=eend-beta*deminend;
				double Deend3Dh=DeendDh-beta*DdeminendDh;
				double sigend2=gamma*elastic*deminend;
				double Dsigend2Dh=gamma*
					(DelasticDh*deminend+elastic*DdeminendDh);
				double a,b,c,d,DaDh,DbDh,DcDh,DdDh;
				polinomial 
				(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d); 
				polinomialwithSens
				(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d,
				 Deend3Dh, 0.0, 0.0, Deend2Dh, Dsigend2Dh, DelasticDh,
				 DaDh, DbDh, DcDh, DdDh);
				DsigDh=DaDh*Tstrain*Tstrain*Tstrain+DbDh*Tstrain*Tstrain
					   +DcDh*Tstrain+DdDh; 
			}else{
				path=6;
				DsigDh=0.0;
			}
		}
	}else if(Cstate==3||Cstate==4||Cstate==5){
		if(Tstate==1){
			path=7;
			envelopeSens 
			(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 DsigDh, DtanDh);
		}else {
			double sigmin,DsigminDemin,DtanminDemin;
			double DsigminDheminFix,DtanminDheminFix;
			double eend,elastic,DeendDemin,DelasticDemin;
			double DeendDheminFix,DelasticDheminFix;
			double DeendDh,DelasticDh;
			envelopewithSens1 
			(CminStrain, sigmin, DsigminDemin, DtanminDemin);
			envelopeSens 
			(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 DsigminDheminFix, DtanminDheminFix);
			unloadwithSens2 
			(CminStrain, sigmin, DsigminDemin, DsigminDheminFix,
			 DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			 eend, elastic, DeendDemin, DelasticDemin, 
			 DeendDheminFix, DelasticDheminFix);
			DeendDh=DeendDheminFix+DeendDemin*DeminDh;
			DelasticDh=DelasticDheminFix+DelasticDemin*DeminDh; 
			double deminend = CminStrain-eend;
			double DdeminendDh = DeminDh-DeendDh;
			if(Tstate==2){
				path=8;
				double emin2=CminStrain-gamma*deminend;
				double Demin2Dh=DeminDh-gamma*DdeminendDh;
				double emin3=CminStrain+beta*deminend;
				double Demin3Dh=DeminDh+beta*DdeminendDh;
				double sigmin2=(1.0-gamma)*elastic*deminend;
				double Dsigmin2Dh=
				(1.0-gamma)*(DelasticDh*deminend+elastic*DdeminendDh);
				double sigmin3, tanmin3, Dtanmin3Demin3;
				double Dsigmin3Dhemin3Fix, Dtanmin3Dhemin3Fix;
				double Dsigmin3Dh, Dtanmin3Dh, a, b, c, d;
				double DaDh, DbDh, DcDh, DdDh;
				envelopewithSens1 
				(emin3, sigmin3, tanmin3, Dtanmin3Demin3);
				envelopeSens 
				(emin3, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			    Dsigmin3Dhemin3Fix, Dtanmin3Dhemin3Fix);
				Dsigmin3Dh=Dsigmin3Dhemin3Fix+tanmin3*Demin3Dh;
				Dtanmin3Dh=Dtanmin3Dhemin3Fix+Dtanmin3Demin3*Demin3Dh;
				polinomial 
				(emin2, sigmin2, elastic, emin3, sigmin3, tanmin3, a, b, c, d); 
				polinomialwithSens
				(emin2, sigmin2, elastic, emin3, sigmin3, tanmin3, a, b, c, d, 
				 Demin2Dh, Dsigmin2Dh, DelasticDh, Demin3Dh, Dsigmin3Dh, Dtanmin3Dh,
				 DaDh, DbDh, DcDh, DdDh);
				DsigDh=DaDh*Tstrain*Tstrain*Tstrain+DbDh*Tstrain*Tstrain
					  +DcDh*Tstrain+DdDh; 
			}else if(Tstate==3){
				path=9;
				DsigDh=DelasticDh*(Tstrain-CendStrain)-elastic*DeendDh;
			}else if(Tstate==4){
				double eend2=eend+gamma*deminend;
				double Deend2Dh=DeendDh+gamma*DdeminendDh;
				double eend3=eend-beta*deminend;
				double Deend3Dh=DeendDh-beta*DdeminendDh;
				double sigend2=gamma*elastic*deminend;
				double Dsigend2Dh=gamma*
					(DelasticDh*deminend+elastic*DdeminendDh);
				double a,b,c,d,DaDh,DbDh,DcDh,DdDh;
				polinomial 
				(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d); 
				polinomialwithSens
				(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d,
				 Deend3Dh, 0.0, 0.0, Deend2Dh, Dsigend2Dh, DelasticDh,
				 DaDh, DbDh, DcDh, DdDh);
				DsigDh=DaDh*Tstrain*Tstrain*Tstrain+DbDh*Tstrain*Tstrain
					  +DcDh*Tstrain+DdDh; 
				path=10;
			}else{
				path=11;
				DsigDh=0.0;
			}
		}
	}else if(Cstate==2){
		if(Tloading){
			if(Cloading){
				if(Tstate==1){
					path=12;
					envelopeSens 
					(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					DsigDh, DtanDh);
				}else{
					path=13;
					double sigmin,tanmin,DtanminDemin;
					double DsigminDheminFix,DtanminDheminFix;
					double eend,elastic,DeendDemin,DelasticDemin;
					double DeendDheminFix,DelasticDheminFix;
					double DeendDh,DelasticDh;
					envelopewithSens1 
					(CminStrain, sigmin, tanmin, DtanminDemin);
					envelopeSens 
					(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					DsigminDheminFix, DtanminDheminFix);
					unloadwithSens2 
					(CminStrain, sigmin, tanmin, DsigminDheminFix,
					DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					eend, elastic, DeendDemin, DelasticDemin, 
					DeendDheminFix, DelasticDheminFix);
					DeendDh=DeendDheminFix+DeendDemin*DeminDh;
					DelasticDh=DelasticDheminFix+DelasticDemin*DeminDh; 
					double deminend = CminStrain-eend;
					double DdeminendDh = DeminDh-DeendDh;
					double emin2=CminStrain-gamma*deminend;
					double Demin2Dh=DeminDh-gamma*DdeminendDh;
					double emin3=CminStrain+beta*deminend;
					double Demin3Dh=DeminDh+beta*DdeminendDh;
					double sigmin2=(1.0-gamma)*elastic*deminend;
					double Dsigmin2Dh=
					(1.0-gamma)*(DelasticDh*deminend+elastic*DdeminendDh);
					double sigmin3, tanmin3, Dtanmin3Demin3;
					double Dsigmin3Dhemin3Fix, Dtanmin3Dhemin3Fix;
					double Dsigmin3Dh, Dtanmin3Dh, a, b, c, d;
					double DaDh, DbDh, DcDh, DdDh;
					envelopewithSens1 
					(emin3, sigmin3, tanmin3, Dtanmin3Demin3);
					envelopeSens 
					(emin3, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				    Dsigmin3Dhemin3Fix, Dtanmin3Dhemin3Fix);
					Dsigmin3Dh=Dsigmin3Dhemin3Fix+tanmin3*Demin3Dh;
					Dtanmin3Dh=Dtanmin3Dhemin3Fix+Dtanmin3Demin3*Demin3Dh;
					polinomial 
					(emin2, sigmin2, elastic, emin3, sigmin3, tanmin3, a, b, c, d); 
					polinomialwithSens
					(emin2, sigmin2, elastic, emin3, sigmin3, tanmin3, a, b, c, d, 
					 Demin2Dh, Dsigmin2Dh, DelasticDh, Demin3Dh, Dsigmin3Dh, Dtanmin3Dh,
					 DaDh, DbDh, DcDh, DdDh);
					DsigDh=DaDh*Tstrain*Tstrain*Tstrain+DbDh*Tstrain*Tstrain
					      +DcDh*Tstrain+DdDh; 
				}
			}else{
				if(Tstate==1){
					path=14;
					envelopeSens 
					(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					 DsigDh, DtanDh);
				}else {
					//// DeminDh is updated here !! /////
					path=15;
					double DaDh,DbDh,DcDh,DdDh;
					solemin1Sens
					(TminStrain, Cstrain, DstressDh, DstrainDh,
					 DfpcDh,DfpcuDh,Depsc0Dh,DepscuDh, 
					 DeminDh, DaDh, DbDh,DcDh,DdDh);
					DeminDhUpdated=DeminDh;
				    DsigDh=DaDh*Tstrain*Tstrain*Tstrain
						  +DbDh*Tstrain*Tstrain
						  +DcDh*Tstrain+DdDh;
				}
			}
		}else{
			if(Cloading){
				if(Tstate==5){
					path=16;
					DsigDh=0.0;
				}else{
					//// DeminDh is updated here !! /////
					double DeendDh,DelasticDh;
					solemin2Sens
					(TminStrain, Cstrain, DstressDh, DstrainDh,
					 DfpcDh,DfpcuDh,Depsc0Dh,DepscuDh, 
					 DeminDh, DeendDh, DelasticDh);
					DeminDhUpdated=DeminDh;
					if(Tstate==2||Tstate==3){
						path=17;
						DsigDh=DelasticDh*(Tstrain-TendStrain)
							  -Telastic*DeendDh;
					}else{
						path=18;
						double deminend=TminStrain-TendStrain;
						double DdeminendDh=DeminDh-DeendDh;
						double eend2=TendStrain+gamma*deminend;
						double Deend2Dh=DeendDh+gamma*DdeminendDh;
						double eend3=TendStrain-beta*deminend;
						double Deend3Dh=DeendDh-beta*DdeminendDh;
						double sigend2=gamma*Telastic*deminend;
						double Dsigend2Dh=
						gamma*(DelasticDh*deminend+Telastic*DdeminendDh);
						double a,b,c,d,DaDh,DbDh,DcDh,DdDh;
						polinomial
						(eend2,sigend2,Telastic,eend3,0.0,0.0, a, b, c, d); 
						polinomialwithSens
						(eend2,sigend2,Telastic,eend3,0.0,0.0, a, b, c, d, 
						Deend2Dh,Dsigend2Dh,DelasticDh,
						Deend3Dh,0.0,0.0,DaDh, DbDh, DcDh, DdDh);
						DsigDh=DaDh*Tstrain*Tstrain*Tstrain
						  +DbDh*Tstrain*Tstrain
						  +DcDh*Tstrain+DdDh;
					}
				}
			}else{
				if(Tstate==5){
					DsigDh=0.0;
					path=19;
				}else{
					double sigmin,tanmin,dtanmindemin;
					double DsigminDheminFix,DtanminDheminFix;
					double eend,elastic,DeendDheminFix,DelasticDheminFix;
					double DeendDh,DelasticDh,DeendDemin,DelasticDemin;
					envelopewithSens1 
					(CminStrain, sigmin, tanmin, dtanmindemin);
					envelopeSens 
					(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					DsigminDheminFix, DtanminDheminFix);
					unloadwithSens2 
					(CminStrain, sigmin, tanmin, DsigminDheminFix,
					DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
					eend, elastic, DeendDemin, DelasticDemin, 
					DeendDheminFix, DelasticDheminFix);
					DeendDh=DeendDheminFix+DeendDemin*DeminDh;
					DelasticDh=DelasticDheminFix+DelasticDemin*DeminDh; 
					if(Tstate==2||Tstate==3){
						DsigDh=DelasticDh*(Tstrain-CendStrain)-elastic*DeendDh;
					    path=20;
					}else if(Tstate==4){
					    path=21;
						double deminend = CminStrain-eend;
						double DdeminendDh = DeminDh-DeendDh;
						double eend2=eend+gamma*deminend;
						double Deend2Dh=DeendDh+gamma*DdeminendDh;
						double eend3=eend-beta*deminend;
						double Deend3Dh=DeendDh-beta*DdeminendDh;
						double sigend2=gamma*elastic*deminend;
						double Dsigend2Dh=gamma*
						(DelasticDh*deminend+elastic*DdeminendDh);
						double a,b,c,d,DaDh,DbDh,DcDh,DdDh;
						polinomial 
						(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d); 
						polinomialwithSens
						(eend3, 0.0, 0.0, eend2, sigend2, elastic, a, b, c, d,
						Deend3Dh, 0.0, 0.0, Deend2Dh, Dsigend2Dh, DelasticDh,
						DaDh, DbDh, DcDh, DdDh);
						DsigDh=DaDh*Tstrain*Tstrain*Tstrain+DbDh*Tstrain*Tstrain
							+DcDh*Tstrain+DdDh; 
					}
				}
			}
		}
	}else if(Cstate==-1){
		if(Tloading){
			if(Tstate==1){
				envelopeSens 
				(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				 DsigDh, DtanDh);
				path=22;
			}else{
				double sigmin,tanmin,DtanminDemin;
				double DsigminDheminFix,DtanminDheminFix;
				double a,b,c,d,DaDh,DbDh,DcDh,DdDh,DsigminDh,DtanminDh;
				envelopewithSens1 
				(CminStrain, sigmin, tanmin, DtanminDemin);
				envelopeSens 
				(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				DsigminDheminFix, DtanminDheminFix);
				DsigminDh=DsigminDheminFix+tanmin*DeminDh;
				DtanminDh=DtanminDheminFix+DtanminDemin*DeminDh;
				polinomial
				(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d); 
				polinomialwithSens
				(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d, 
				 -DeminDh,   0.0,0.0,DeminDh   ,DsigminDh,DtanminDh, 
				 DaDh, DbDh, DcDh, DdDh);
				DsigDh=DaDh*Tstrain*Tstrain*Tstrain
						  +DbDh*Tstrain*Tstrain
						  +DcDh*Tstrain+DdDh;
				path=23;
			}
		}else{
			DsigDh=0.0;
			path=24;
		}
	}else if(Cstate==-2){
		if(Tloading){
			if(Tstate==1){
				path=25;
				envelopeSens 
				(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				 DsigDh, DtanDh);
			}else if(Tstate==-3){
				path=26;
				double sigmin,tanmin,DtanminDemin;
				double DsigminDheminFix,DtanminDheminFix;
				double a,b,c,d,DaDh,DbDh,DcDh,DdDh,DsigminDh,DtanminDh;
				envelopewithSens1 
				(CminStrain, sigmin, tanmin, DtanminDemin);
				envelopeSens 
				(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
				DsigminDheminFix, DtanminDheminFix);
				DsigminDh=DsigminDheminFix+tanmin*DeminDh;
				DtanminDh=DtanminDheminFix+DtanminDemin*DeminDh;
				polinomial
				(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d); 
				polinomialwithSens
				(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d, 
				 -DeminDh,   0.0,0.0,DeminDh   ,DsigminDh,DtanminDh, 
				 DaDh, DbDh, DcDh, DdDh);
				DsigDh=DaDh*Tstrain*Tstrain*Tstrain
						  +DbDh*Tstrain*Tstrain
						  +DcDh*Tstrain+DdDh;
			}
		}else{
			path=27;
			DsigDh=0.0;
		}
	}else if(Cstate==-3){
		if(Tstate==1){
			envelopeSens 
			(Tstrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			DsigDh, DtanDh);
			path=28;
		}else if(Tstate==-3){
			double sigmin,tanmin,DtanminDemin;
			double DsigminDheminFix,DtanminDheminFix;
			double a,b,c,d,DaDh,DbDh,DcDh,DdDh,DsigminDh,DtanminDh;
			envelopewithSens1 
			(CminStrain, sigmin, tanmin, DtanminDemin);
			envelopeSens 
			(CminStrain, DfpcDh, DfpcuDh, Depsc0Dh, DepscuDh,
			DsigminDheminFix, DtanminDheminFix);
			DsigminDh=DsigminDheminFix+tanmin*DeminDh;
			DtanminDh=DtanminDheminFix+DtanminDemin*DeminDh;
			polinomial
			(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d); 
			polinomialwithSens
			(-CminStrain,0.0,0.0,CminStrain,sigmin,tanmin,a,b,c,d, 
			 -DeminDh,   0.0,0.0,DeminDh   ,DsigminDh,DtanminDh, 
			 DaDh, DbDh, DcDh, DdDh);
			DsigDh=DaDh*Tstrain*Tstrain*Tstrain
					  +DbDh*Tstrain*Tstrain
					  +DcDh*Tstrain+DdDh;
			path=29;
		}else{
			path=30;
			DsigDh=0.0;
		}
	}
	DTstressDhTstrainFix=DsigDh;

//	opserr<< "path=          "<<path<<"\n";
	return DsigDh;

}

void SmoothConcrete02::envelope 
(const double& eps, double& sig, double& tan)
{
	if(eps>0.0){
		sig=0.0;
		tan=0.0;
	}else if(eps>=epsc0){
		double rrr=eps/epsc0;
		sig=fpc*rrr*(2.0-rrr);
		tan=2.0*fpc/epsc0*(1.0-rrr);
	}else if(eps>=epscu){
		double a0,b0,c0,d0;
		polinomial(epsc0, fpc, 0.0,  epscu, fpcu, 0.0,
				   a0,b0,c0,d0);
		double eps2=eps*eps;
		double eps3=eps2*eps;
		sig=a0*eps3+b0*eps2+c0*eps+d0;
		tan=3.0*a0*eps2+2.0*b0*eps+c0;
	}else{
		sig=fpcu;
		tan=0.0;
	}
}
void SmoothConcrete02::envelopewithSens1 
(const double& eps, double& sig, double& tan, double& DtanDeps)
{
	if(eps>0.0){
		sig=0.0;
		tan=0.0;
		DtanDeps=0.0;
	} else if(eps>=epsc0){
		double rrr=eps/epsc0;
		sig=fpc*rrr*(2.0-rrr);
		tan=2.0*fpc/epsc0*(1-rrr);
		DtanDeps=-2.0*fpc/(epsc0*epsc0);
	}else if(eps>=epscu){
		double a0,b0,c0,d0;
		double eps2=eps*eps;
		double eps3=eps2*eps;
		polinomial(epsc0, fpc, 0.0,  epscu, fpcu, 0.0,
				   a0,b0,c0,d0);
		sig=a0*eps3+b0*eps2+c0*eps+d0;
		tan=3.0*a0*eps2+2.0*b0*eps+c0;
		DtanDeps=6.0*a0*eps+2.0*b0;
	}else{
		sig=fpcu;
		tan=0.0;
		DtanDeps=0.0;
	}
}
void SmoothConcrete02::envelopeSens 
(const double& eps, const double& DfpcDh, const double& DfpcuDh, 
 const double& Depsc0Dh, const double& DepscuDh,
 double& DsigDh, double& DtanDh)
{
	if(eps>0.0){
		DsigDh=0.0;
		DtanDh=0.0;
	}else if(eps>=epsc0){
		double rrr=eps/epsc0;
		DsigDh=DfpcDh*rrr*(2.0-rrr)
               +2.0*fpc*rrr*Depsc0Dh/epsc0*(rrr-1.0);
		DtanDh=2.0/epsc0*
		(DfpcDh*(1.0-rrr)-Depsc0Dh*fpc/epsc0*(1.0-2.0*rrr));
	}else if(eps>=epscu){
		double a0, b0, c0, d0;
		double Da0Dh, Db0Dh, Dc0Dh, Dd0Dh;
		polinomial(epsc0, fpc, 0.0,  epscu, fpcu, 0.0,
				   a0,b0,c0,d0);
		polinomialwithSens
		(epsc0, fpc, 0.0, epscu, fpcu, 0.0, a0, b0, c0, d0, 
		 Depsc0Dh, DfpcDh, 0.0, DepscuDh, DfpcuDh, 0.0,
		 Da0Dh, Db0Dh, Dc0Dh, Dd0Dh);
		double eps2=eps*eps;
		double eps3=eps2*eps;
		DsigDh=Da0Dh*eps3+Db0Dh*eps2+Dc0Dh*eps+Dd0Dh;
		DtanDh=3.0*Da0Dh*eps2+2.0*Db0Dh*eps+Dc0Dh;
	}else{
		DsigDh=DfpcuDh;
		DtanDh=0.0;
	}
}
int
SmoothConcrete02::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	double DeminDh,DstressDh,DstrainDh;

	if (SHVs == 0) {
		SHVs = new Matrix(3,numGrads);
		DeminDh = 0.0;
		DstressDh = 0.0;
		DstrainDh = 0.0;
//		output.open("SmoothConcrete02.txt", ios::out);
//		output << "\n";
//		output.flush();
	}else{
		DeminDh   = (*SHVs)(0,(gradNumber-1));
		DstressDh = (*SHVs)(1,(gradNumber-1));
		DstrainDh = (*SHVs)(2,(gradNumber-1));
	}

	if(gradNumber!=UnconditionalDerivative){
		this->getStressSensitivity(gradNumber,false);
	}

	DstressDh=DTstressDhTstrainFix+Ttangent*TstrainSensitivity;
	DstrainDh=TstrainSensitivity;
////// update DeminDh //////
	if(Cstate==0){
		if(Tloading){
			DeminDh=TstrainSensitivity;  
		}else{
			if(Tstate==-1){
				DeminDh=TstrainSensitivity;  
			}else{
				DeminDh=0.0;
			}
		}
	}else if(Cstate==1){
		if(Tloading){
			DeminDh=TstrainSensitivity;  
		}
	}else if(Cstate==3||Cstate==4||Cstate==5){
		if(Tstate==1){
			DeminDh=TstrainSensitivity;  
		}
	}else if(Cstate==2){
		if(Tloading){
			if(Cloading){
				if(Tstate==1){
					DeminDh=TstrainSensitivity;  
				}
			}else{
				if(Tstate==1){
					DeminDh=TstrainSensitivity;  
				}else{
					DeminDh=DeminDhUpdated;
				}
			}
		}else{
			if(Cloading){
				DeminDh=DeminDhUpdated;
			}
		}
	}else if(Cstate==-1){
		if(Tloading){
			if(Tstate==1){
				DeminDh=TstrainSensitivity;  
			}
		}else{
			if(Tstate==1){
				DeminDh=TstrainSensitivity;  
			}
		}
	}else if(Cstate==-2){
		if(Tloading){
			if(Tstate==1){
				DeminDh=TstrainSensitivity;  
			}
		}
	}else if(Cstate==-3){
		if(Tstate==1){
			DeminDh=TstrainSensitivity;  
		}else if(Tstate==-1){
			DeminDh=-TstrainSensitivity;  
		}else if(Tstate==-2){
			DeminDh=0.0;
		}
	}

	if(abs(DstressDh)>1.0e5){
		double abc=0.0;
	}

	(*SHVs)(0,(gradNumber-1))=DeminDh;
	(*SHVs)(1,(gradNumber-1))=DstressDh;
	(*SHVs)(2,(gradNumber-1))=DstrainDh;

/*	if(gradNumber==numGrads){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		for(int i=0;i<gradNumber;i++){
			output<<setw(15)<<setprecision(5)<<(*SHVs)(0,i);
		}
		for(i=0;i<gradNumber;i++){
			output<<setw(15)<<setprecision(5)<<(*SHVs)(1,i);
		}
		for(i=0;i<gradNumber;i++){
			output<<setw(15)<<setprecision(5)<<(*SHVs)(2,i);
		}
		output<<"\n";
		output.flush();
		output<<"\n";
		output<<"\n";
		output<<"path1=  "<<path1<<"\n";
		output<<"path2=  "<<path2<<"\n";
		output<<"path3=  "<<path3<<"\n";
		output<<"path4=  "<<path4<<"\n";
		output<<"path5=  "<<path5<<"\n";
		output<<"path6=  "<<path6<<"\n";
		output<<"path7=  "<<path7<<"\n";
		output<<"path8=  "<<path8<<"\n";
		output<<"path9=  "<<path9<<"\n";
		output<<"path10= "<<path10<<"\n";
		output<<"path11=  "<<path11<<"\n";
		output<<"path12=  "<<path12<<"\n";
		output<<"path13=  "<<path13<<"\n";
		output<<"path14=  "<<path14<<"\n";
		output<<"path15=  "<<path15<<"\n";
		output<<"path16=  "<<path16<<"\n";
		output<<"path17=  "<<path17<<"\n";
		output<<"path18=  "<<path18<<"\n";
		output<<"path19=  "<<path19<<"\n";
		output<<"path20=  "<<path20<<"\n";
		output<<"path21=  "<<path21<<"\n";
		output<<"path22=  "<<path22<<"\n";
		output<<"path23=  "<<path23<<"\n";
		output<<"path24=  "<<path24<<"\n";
		output<<"path25=  "<<path25<<"\n";
		output<<"path26=  "<<path26<<"\n";
		output<<"path27=  "<<path27<<"\n";
		output<<"path28=  "<<path28<<"\n";
		output<<"path29=  "<<path29<<"\n";
		output<<"path30=  "<<path30<<"\n";
	}
*/
	return 0;
}



