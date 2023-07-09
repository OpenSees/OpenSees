// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete05.cpp,v $
// Created: 07/04
// Written: KO 
// Description: This file contains the class implementation for 
// Concrete05. 

#include <Concrete05.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>


Concrete05::Concrete05
(int tag, double FPCC, double EPCC, double EC, double RC, double XCRN, double FT, double ET, double RT, double XCRP)
  :UniaxialMaterial(tag, MAT_TAG_Concrete05),
   fpcc(FPCC), epcc(EPCC), Ec(EC), rc(RC), xcrn(XCRN), ft(FT), et(ET), rt(RT), xcrp(XCRP),// input
   Ceunn(0.0), Cfunn(0.0), Ceunp(0.0), Cfunp(0.0), Cer(0.0), Cfr(0.0), Cer0n(0.0), 
   Cfr0n(0.0), Cer0p(0.0), Cfr0p(0.0), Ce0(0.0), Cea(0.0), Ceb(0.0), Ced(0.0), Cinc(0.0), 
   Crule(0.0), Cstrain(0.0), Cstress(0.0), Ctangent(0.0) // history and state variables

// Come back to rule on previous line

{

// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

}

Concrete05::Concrete05():UniaxialMaterial(0, MAT_TAG_Concrete05),
   fpcc(0.0), epcc(0.0), Ec(0.0), rc(0.0), xcrn(0.0), ft(0.0), et(0.0), rt(0.0), xcrp(0.0),// input
   Ceunn(0.0), Cfunn(0.0), Ceunp(0.0), Cfunp(0.0), Cer(0.0), Cfr(0.0), Cer0n(0.0), 
   Cfr0n(0.0), Cer0p(0.0), Cfr0p(0.0), Ce0(0.0), Cea(0.0), Ceb(0.0), Ced(0.0), Cinc(0.0), 
   Crule(0.0), Cstrain(0.0), Cstress(0.0), Ctangent(0.0) // history and state variables

// Last four lines and the comma just before are not necessary

{

// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Concrete05::~Concrete05 ()
{
   // Does nothing
}

int Concrete05::setTrialStrain (double strain, double strainRate)

{
  // Set trial strain
  this->revertToLastCommit();
  Tstrain = strain;

  if (Cinc==0.0)	{			// monotonic, first data point
	  if (Tstrain < 0.0)   {	        // negative envelope,	  
		  fcEtnf(Tstrain);
		  Tinc=-1.0;
//		  opserr << "e" << Tstrain << "\n";
	  }	  
	  else if (Tstrain > 0.0)	{  
		  fcEtpf(Tstrain,Ce0);
		  Tinc=1.0;
	  }
	  else	{
		  Tstress=0.0;
		  Ttangent=Ec;
		  Trule=0.0;
		  Tinc=0.0;
	  }
	  
	  Teunn=0.0;	
	  Tfunn=0.0;
	  Teunp=0.0;
	  Tfunp=0.0;
	  Ter=0.0;
	  Tfr=0.0;
	  Ter0n=0.0;
	  Tfr0n=0.0;
	  Ter0p=0.0;
	  Tfr0p=0.0;
	  Te0=0.0;
	  Tea=0.0;
	  Teb=0.0;
	  Ted=0.0;
  }

  else	{	//cyclic

	  if (Tstrain > Cstrain)	{
		  Tinc=1.0;
	  }
	  else if (Tstrain < Cstrain)	{
		  Tinc=-1.0;
	  }
	  else	{
		  Tinc=Cinc;
	  }

	  Teunn=Ceunn;
	  Tfunn=Cfunn;
	  Teunp=Ceunp;
	  Tfunp=Cfunp;
	  Ter=Cer;
	  Tfr=Cfr;
	  Ter0n=Cer0n;
	  Tfr0n=Cfr0n;
	  Ter0p=Cer0p;
	  Tfr0p=Cfr0p;
	  Te0=Ce0;
	  Tea=Cea;
	  Teb=Ceb;
	  Ted=Ced;
	  Trule=Crule;
	  
	  esplnf(Teunn,Tfunn);
	  Eplnf(Teunn);
	  Esecnf(Teunn,Tfunn);
	  delenf(Teunn);
	  delfnf(Teunn,Tfunn);
	  fnewnf(Teunn,Tfunn);
	  Enewnf(Teunn,Tfunn);
	  esrenf(Teunn);
	  freErenf(Teunn);
	  fnewstnf(Tfunn,delfn,Teunn,Ter0n,espln);
	  Enewstnf(fnewstn,Tfr0n,Teunn,Ter0n);
//	  opserr << "esrestn" << esrestn << "\n";
	  esrestnf(Teunn,delen,Ter0n,espln);   
//	  opserr << "esrestn" << esrestn << "\n";
	  freErestnf(Teunn,Tfunn,Ter0n);

	  esplpf(Teunp,Tfunp,Te0,espln);
	  Eplpf(Te0,Teunp);
	  Esecpf(Te0,Teunp,Tfunp,espln);  
	  delepf(Teunp,Te0);
	  delfpf(Tfunp,Teunp,Te0); 
	  fnewpf(Tfunp,Teunp,Te0);
	  Enewpf(Teunp,Tfunp,Te0,espln);
	  esrepf(Teunp,Te0);
	  freErepf(Teunp,Te0);
	  fnewstpf(Tfunp,delfp,Teunp,Ter0p,esplp,Te0);
	  Enewstpf(fnewstp,Tfr0p,Teunp,Ter0p);
	  esrestpf(Teunp,delep,Ter0p,esplp);
	  freErestpf(Teunp,Tfunp,Ter0p,Te0,espln);

	  if (Cinc==-1.0)	{

		  if (Tstrain > Cstrain)	{	// Start reversal from negative direction to positive direction		
			  
			  if (Crule==1.0 || Crule==5.0 || Crule==7.0)	{	// Rules [3,9,8,2,6]
			  
				  Teunn=Cstrain;
				  Tfunn=Cstress;
				  e0eunpfunpf(Te0,Teunp,Tfunp,Teunn,Tfunn);
				  esplnf(Teunn,Tfunn);
			  	  Eplnf(Teunn);			
//				  Esecpf(Te0,Teunp,Tfunp,espln);	
			  	  fnewpf(Tfunp,Teunp,Te0);
			  	  Enewpf(Teunp,Tfunp,Te0,espln);
			  	  esrepf(Teunp,Te0);
			  	  freErepf(Teunp,Te0);
				  
				  if (Tstrain<=espln)	{			// Rule 3
					  r3f(Teunn,Tfunn,espln,Epln);
				  	  Trule=3.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
				  	  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain<=Teunp)	{		// Rule 9
				  	  r9f(espln,Epln,Teunp,fnewp,Enewp);
				  	  Trule=9.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain<=esrep)	{		// Rule 8
				  	  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
				      Trule=8.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else	{							// Rules 2 and 6
				  	  fcEtpf(Tstrain,Te0);
			  	  }
	
//			  	  Tstress=0.0;
//			  	  Ttangent=0.0;
//			  	  Teunn=0.0;	
//			  	  Tfunn=0.0;
//			  	  Teunp=0.0;
//			  	  Tfunp=0.0;
//			  	  Ter=0.0;
//			  	  Tfr=0.0;
//			  	  Ter0n=0.0;
//			  	  Tfr0n=0.0;
//			  	  Ter0p=0.0;
//			  	  Tfr0p=0.0;
//			  	  Te0=0.0;
//			  	  Tea=0.0;
//			  	  Teb=0.0;
//			  	  Ted=0.0;
		  
			  }	// if (Crule=1.0) or 5 or 7

			  else if (Crule==10.0)	{

				  Ter=Cstrain;
				  Tfr=Cstress;

				  Teb=Ter;
				  fb=Tfr;

				  ea1112f(Teb,espln,esplp,Teunn,Teunp);

				  if (Tea<=espln)	{

					  if (Tstrain<=Tea)	{				// Rule 12 targeting for ea on 3
						  r3f(Teunn,Tfunn,espln,Epln);	// 
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
						  fca=fc;
						  Eta=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Tea;
						  RAf(esi,fi,Ei,esf,fca,Eta);
						  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
						  Trule=12.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=espln)	{	// Rule 3
						  r3f(Teunn,Tfunn,espln,Epln);
						  Trule=3.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=Teunp)	{	// Rule 9
						  r9f(espln,Epln,Teunp,fnewp,Enewp);
						  Trule=9.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=esrep)	{	// Rule 8						  
						  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  Trule=8.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 2 and 6
						  fcEtpf(Tstrain,Te0);
					  }
				  }
				  
				  else if (Tea<=Teunp)	{			// and Tea>espln

					  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 9
						  r9f(espln,Epln,Teunp,fnewp,Enewp); 
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
						  fca=fc;
						  Eta=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Tea;
						  RAf(esi,fi,Ei,esf,fca,Eta);
						  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
						  Trule=12.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=Teunp)	{	// Rule 9
						  r9f(espln,Epln,Teunp,fnewp,Enewp);
						  Trule=9.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=esrep)	{	// Rule 8						  
						  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  Trule=8.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 2 and 6
						  fcEtpf(Tstrain,Te0);
					  }
				  }

				  else if (Tea<=esrep)	{	

					  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 8
						  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
						  fca=fc;
						  Eta=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Tea;
						  RAf(esi,fi,Ei,esf,fca,Eta);
						  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
						  Trule=12.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain<=esrep)	{	// Rule 8						  
						  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  Trule=8.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 2 and 6
						  fcEtpf(Tstrain,Te0);
					  }
				  }

				  else	{		// (Tea>esrep)

					  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 2 or 6
						  fcEtpf(Tea,Te0);
						  fca=Tstress;
						  Eta=Ttangent;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Tea;
						  RAf(esi,fi,Ei,esf,fca,Eta);
						  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
						  Trule=12.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;					  
					  }
					  else	{						// Rules 2 and 6
						  fcEtpf(Tstrain,Te0);
					  }
				  }
			  } // if (Crule=10.0)
			  
			  else if (Crule==11.0)	{

				  Ter=Cstrain;
				  Tfr=Cstress;

				  if (Teb!=Ter0p)	{

					  if (Tea<=espln)	{

						  if (Tstrain<=Tea)	{	// Rule 12 targeting for ea on 3							  
							  r3f(Teunn,Tfunn,espln,Epln);   
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
						  	  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }
						  else if (Tstrain<=espln)	{	// Rule 3
							  r3f(Teunn,Tfunn,espln,Epln);
						  	  Trule=3.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9
						  	  r9f(espln,Epln,Teunp,fnewp,Enewp);
						  	  Trule=9.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;		  
						  }
					  	  else if (Tstrain<=esrep)	{	// Rule 8						  
						 	  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  	  Trule=8.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
					  	  }		  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }

					  else if (Tea<=Teunp)	{			// and Tea>espln

					  	  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;			  
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp);
							  Trule=9.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain<=esrep)	{	// Rule 8						  
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 2 and 6			  
							  fcEtpf(Tstrain,Te0);		  
						  }
					  }

					  else if (Tea<=esrep)	{	

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 8
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else if (Tstrain<=esrep)	{	// Rule 8						  
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }					  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }				  
					  }

					  else	{		// (Tea>esrep)

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 2 or 6
							  fcEtpf(Tea,Te0);
							  fca=Tstress;
							  Eta=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 2 and 6
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }
				  }

				  else	{	// if (Teb==Ter0p)

					  if (Tea<=esrestp)	{

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 88
							  r88f(Tea,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain<esrestp)	{	// Rule 88
							  r88f(Tstrain,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
							  Trule=88.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }

					  else	{	// if (Tea>esrestp)

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 2 or 6
							  fcEtpf(Tea,Te0);
							  fca=Tstress;
							  Eta=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 2 and 6
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }
				  }
			  } // if (Crule==11.0)

			  else if (Crule==13.0 || Crule==15.0)	{

				  Ter=Cstrain;
				  Tfr=Cstress;

				  if (Crule==13.0)	{
					  Tea=Ter;
					  fa=Tfr;
					  eb1415f(Tea,fa,Esecn);
				  }
				  else if (Crule==15.0)	{
					  Tea=Cea;
					  Teb=Ceb;
				  }

				  if (Tstrain<=Teb)	{					// Rule 14
					  r14f(Ter,Tfr,Teb);
					  Trule=14.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain<Teunp)	{			// Rule 66
					  r66f(Tstrain,Te0);
					  Trule=66.0;
				  }
				  else	{								// Rule 6
					  fcEtpr6f(Tstrain,Te0);
				  }
			  } // if (Crule==13.0) or 15.0

			  else if (Crule==4.0)	{

				  Ter0p=Cstrain;
				  Tfr0p=Cstress;
				  Teb=Ter0p;

				  fnewstpf(Tfunp,delfp,Teunp,Ter0p,esplp,Te0);
				  Enewstpf(fnewstp,Tfr0p,Teunp,Ter0p);
				  esrestpf(Teunp,delep,Ter0p,esplp);
				  freErestpf(Teunp,Tfunp,Ter0p,Te0,espln);

				  if (Tstrain<esrestp)	{				// Rule 88
					  r88f(Tstrain,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
					  Trule=88.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else	{								// Rules 2 and 6
					  fcEtpf(Tstrain,Te0);
				  }
			  } // if (Crule==4.0)

			  else if (Crule==77.0)	{					// Reversal from transition 77 [Rules 12,(3,9 or 9),8,2,6 or Rules 3,9,8,2,6] 

				  if (Cstrain>=Teunn)	{			
					  Ter=Cstrain;
					  Tfr=Cstress;
					  Teb=Ter;
					  Tea=Ter0n;

					  if (Tea<=espln)	{				// Reversal from 77 by Rules [12,3,9,8,2,6]	
					  	  
						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 3
							  r3f(Teunn,Tfunn,espln,Epln);   
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
						  	  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;						  
						  }					   
						  else if (Tstrain<=espln)	{	// Rule 3
							  r3f(Teunn,Tfunn,espln,Epln);
						  	  Trule=3.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9
						  	  r9f(espln,Epln,Teunp,fnewp,Enewp);
						  	  Trule=9.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
						  }
					  	  else if (Tstrain<=esrep)	{	// Rule 8						 
						 	  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  	  Trule=8.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
					  	  }		  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }

					  else if (Tea<=Teunp)	{			// Reversal from 77 by Rules [12,9,8,2,6]

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;			  
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp);
							  Trule=9.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain<=esrep)	{	// Rule 8						  
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 2 and 6			  
							  fcEtpf(Tstrain,Te0);		  
						  }
					  }

					  else if (Tea<=esrep)	{			// Reversal from 77 by Rules [12,8,2,6]

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 8
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }
						  else if (Tstrain<=esrep)	{	// Rule 8
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }					  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }				  
					  }

					  else	{		// (Tea>esrep)		// Reversal from 88 by Rules [12,2,6]

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 2 or 6
							  fcEtpf(Tea,Te0);
							  fca=Tstress;
							  Eta=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }				
						  else	{						// Rules 2 and 6
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }				  
				  }	// if (Cstrain>=Teunn)
				  
				  else	{	// if (Cstrain<Teunn)		// Reversal from transition 77 by Rules [3,9,8,2,6] 

					  Teunn=Cstrain;				  
					  Tfunn=Cstress;
					  e0eunpfunpf(Te0,Teunp,Tfunp,Teunn,Tfunn);
					  esplnf(Teunn,Tfunn);
					  Eplnf(Teunn);			
//					  Esecpf(Te0,Teunp,Tfunp,espln);	// careful here (don't need)
					  fnewpf(Tfunp,Teunp,Te0);
					  Enewpf(Teunp,Tfunp,Te0,espln);
					  esrepf(Teunp,Te0);
					  freErepf(Teunp,Te0);

					  if (Tstrain<=espln)	{			// Rule 3	  
						  r3f(Teunn,Tfunn,espln,Epln);
						  Trule=3.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }					  
					  else if (Tstrain<=Teunp)	{		// Rule 9						  
						  r9f(espln,Epln,Teunp,fnewp,Enewp);						  	  
						  Trule=9.0;						  	  
						  RAf(esi,fi,Ei,esf,ff,Ef);					  	  
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R); 	  
						  Tstress=fc; 	  
						  Ttangent=Et;		  	  
					  }
					  else if (Tstrain<=esrep)	{		// Rule 8						  	  
						  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);  
						  Trule=8.0;	  
						  RAf(esi,fi,Ei,esf,ff,Ef);	  
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R); 	  
						  Tstress=fc;  	  
						  Ttangent=Et;						  
					  }						  
					  else	{							// Rules 2 and 6						  		  
						  fcEtpf(Tstrain,Te0);					  	  
					  } 				  
				  }	// if (Cstrain<Teunn)			  
			  }	// if (Crule==77.0)
		  } // if (Tstrain>Cstrain)

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

		  else	{	//	if (Tstrain<=Cstrain)	// Continue going to negative direction

			  if (Crule==4.0 || Crule ==10.0 || Crule==7.0 )	{	// or 10.0 or 7.0

				  if (Tstrain>=esplp)	{
					  r4f(Teunp,Tfunp,esplp,Eplp);
				  	  Trule=4.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
				  	  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain>=Teunn)	{
				  	  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
				  	  Trule=10.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain>=esren)	{
				  	  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
				      Trule=7.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else	{
				  	  fcEtnf(Tstrain);
			  	  }
			  
			  }	// if (Crule==4.0)	// or 10.0 or 7.0

			  else if (Crule==1.0 || Crule==5.0)	{	// or 5.0

				  fcEtnf(Tstrain);

			  }	// if (Crule==1.0)	// or 5.0 or 7.0

			  else if (Crule==77.0)	{			// Continue on transition 77 [Rules 77,1,5]

//				  opserr << "esrestn" << esrestn << "\n";
//				  opserr << "espln" << espln << "\n";

				  if (Tstrain>esrestn)	{		// Rule 77
					  r77f(Tstrain,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
					  Trule=77.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }

				  else	{						// Rules 1 and 5						  
					  fcEtnf(Tstrain);	
				  }
	  
			  }	// if (Crule==77.0)

			  else if (Crule==13.0)	{			// Continue on transition 13 [Rules 13,7,1,5]
				  
				  if (Tstrain>=Teunn)	{	// Rule 13
					  r13f(Ted,Teunn,fnewn,Enewn);
					  Trule=13.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain>=esren)	{	// Rule 7
					  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
					  Trule=7.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else	{						// Rules 1 and 5
					  fcEtnf(Tstrain);					
				  }

			  }	// if (Crule==13.0)

			  else if (Crule==11.0)	{			// Continue on transition 11 [[Rules 11,(4,10 or 10),7,1,5 or Rules 11,77,1,5] 

				  if (Tea!=Ter0n)	{

					  if (Teb>=esplp)	{

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 4
							  r4f(Teunp,Tfunp,esplp,Eplp);	
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esplp)	{	// Rule 4
							  r4f(Teunp,Tfunp,esplp,Eplp);
							  Trule=4.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=Teunn)	{			// and Teb<esplp

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=esren)	{	

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 7
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else	{		// (Teb<esren)

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
							  fcEtnf(Teb);
							  fcb=Tstress;
							  Etb=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Eta,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;					  
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }
				  }
				  
				  else	{	// if (Tea==Ter0n)

					  if (Teb>=esrestn)	{

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 77
							  r77f(Teb,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>esrestn)	{	// Rule 77
							  r77f(Tstrain,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
							  Trule=77.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5						  
							  fcEtnf(Tstrain);					  
						  }
					  }

					  else	{	// if (Teb<esrestn)

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
							  fcEtnf(Teb);
							  fcb=Tstress;
							  Etb=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);					  
						  }
					  }
				  }	
			  } // if (Crule==11.0)

			  else if (Crule==15.0)	{					// Continue on transition 15 [Rules 15,13,7,1,5]
				  
				  if (Tstrain>=Tea)	{					// Rule 15 targeting for ea (ed) on 13
					  r13f(Ted,Teunn,fnewn,Enewn);
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
					  fca=fc;
					  Eta=Et;
					  esi=Ter;
					  fi=Tfr;
					  Ei=Ec;
					  esf=Tea;
					  RAf(esi,fi,Ei,esf,fca,Eta);
					  r15f(Ter,Tfr,Tea,fca,Eta,A,R);
					  Trule=15.0;
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain>=Teunn)	{			// Rule 13
					  r13f(Ted,Teunn,fnewn,Enewn);
					  Trule=13.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain>=esren)	{			// Rule 7
					  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
					  Trule=7.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else	{								// Rules 1 and 5
					  fcEtnf(Tstrain);					
				  }
			  }	// if (Crule==15.0)
		  }		// if (Tstrain<=Cstrain)
	  }			// if (Cinc==-1.0)

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

	  else	{	// if (Cinc==1.0)

		  if (Tstrain<Cstrain)	{		// Start reversal from positive direction to negative direction								

			  if (Crule==2.0 || Crule==8.0)	{	// Rules [4,10,7,1,5]

				  Teunp=Cstrain;
				  fcEtpf(Teunp,Te0);	// need Tfunp only from fcEtpf
				  Tfunp=Tstress;
				  Esecpf(Te0,Teunp,Tfunp,espln);
				  esplpf(Teunp,Tfunp,Te0,espln);
				  Eplpf(Te0,Teunp);

				  if (Tstrain>=esplp)	{		// Rule 4
					  r4f(Teunp,Tfunp,esplp,Eplp);
				  	  Trule=4.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
				  	  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain>=Teunn)	{	// Rule 10
				  	  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
				  	  Trule=10.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain>=esren)	{	// Rule 7
				  	  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
				      Trule=7.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else	{						// Rules 1 and 5
				  	  fcEtnf(Tstrain);
			  	  }

			  }	// if (Crule==2.0) or 6.0 or 8.0

			  else if (Cstress==0.0)	{		// Gap Model

				  Teunp=Cstrain;
				  Tfunp=Cstress;
				  fcEtnf(Teunn);
				  Tfunn=Tstress;
				  Ter=Cstrain;
				  Tfr=Cstress;
				  Ted=Ter;
				  fnewnf(Teunn,Tfunn);
				  Enewnf(Teunn,Tfunn);

				  if (Tstrain>=Teunn)	{		// Rule 13
					  r13f(Ted,Teunn,fnewn,Enewn);
				      Trule=13.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
				  }
				  else if (Tstrain>=esren)	{	// Rule 7
					  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
				      Trule=7.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
				  }
				  else	{						// Rules 1 and 5
					  fcEtnf(Tstrain);
				  }
			  }	// if (Cstress=0.0)

			  else if (Crule==9.0)	{

				  Ter=Cstrain;
				  Tfr=Cstress;
				  Tea=Ter;
				  fa=Tfr;
				  eb1112f(Tea,espln,esplp,Teunn,Teunp);

				  if (Teb>=esplp)	{

					  if (Tstrain>=Teb)	{				// Rule 11 targeting for eb on 4
						  r4f(Teunp,Tfunp,esplp,Eplp);	// 
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
						  fcb=fc;
						  Etb=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Teb;
						  RAf(esi,fi,Ei,esf,fcb,Etb);
						  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
						  Trule=11.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=esplp)	{	// Rule 4
						  r4f(Teunp,Tfunp,esplp,Eplp);
						  Trule=4.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=Teunn)	{	// Rule 10
						  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
						  Trule=10.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=esren)	{	// Rule 7						  
						  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
						  Trule=7.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 1 and 5
						  fcEtnf(Tstrain);
					  }
				  }
				  
				  else if (Teb>=Teunn)	{			// and Teb<esplp

					  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 10
						  r10f(esplp,Eplp,Teunn,fnewn,Enewn); 
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
						  fcb=fc;
						  Etb=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Teb;
						  RAf(esi,fi,Ei,esf,fcb,Etb);
						  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
						  Trule=11.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=Teunn)	{	// Rule 10
						  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
						  Trule=10.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;;
					  }
					  else if (Tstrain>=esren)	{	// Rule 7						  
						  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
						  Trule=7.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 1 and 5
						  fcEtnf(Tstrain);
					  }
				  }

				  else if (Teb>=esren)	{	

					  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 7
						  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
						  fcb=fc;
						  Etb=Et;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Teb;
						  RAf(esi,fi,Ei,esf,fcb,Etb);
						  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
						  Trule=11.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=esren)	{	// Rule 7						  
						  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
						  Trule=7.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{						// Rules 1 and 5
						  fcEtnf(Tstrain);
					  }
				  }

				  else	{		// (Teb<esren)

					  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
						  fcEtnf(Teb);
						  fcb=Tstress;
						  Etb=Ttangent;
						  esi=Ter;
						  fi=Tfr;
						  Ei=Ec;
						  esf=Teb;
						  RAf(esi,fi,Ei,esf,fcb,Etb);
						  r11f(Ter,Tfr,Teb,fcb,Eta,A,R);
						  Trule=11.0;
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;					  
					  }
					  else	{						// Rules 1 and 5
						  fcEtnf(Tstrain);
					  }
				  }
			  } // if (Crule=9.0)
			  
			  else if (Crule==12.0)	{

				  Ter=Cstrain;
				  Tfr=Cstress;

				  if (Tea!=Ter0n)	{

					  if (Teb>=esplp)	{

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 4
							  r4f(Teunp,Tfunp,esplp,Eplp);	
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esplp)	{	// Rule 4
							  r4f(Teunp,Tfunp,esplp,Eplp);
							  Trule=4.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=Teunn)	{			// and Teb<esplp

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=esren)	{	

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 7
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else	{		// (Teb<esren)

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
							  fcEtnf(Teb);
							  fcb=Tstress;
							  Etb=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Eta,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;					  
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }
				  }
				  
				  else	{	// if (Tea==Ter0n)

					  if (Teb>=esrestn)	{

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 77
							  r77f(Teb,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>esrestn)	{	// Rule 77
							  r77f(Tstrain,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
							  Trule=77.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5						  
							  fcEtnf(Tstrain);					  
						  }
					  }

					  else	{	// if (Teb<esrestn)

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
							  fcEtnf(Teb);
							  fcb=Tstress;
							  Etb=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);					  
						  }
					  }
				  }	
			  } // if (Crule==12.0)

			  else if (Crule==14.0)	{			// Reversal from transition 14 [Rules 15,13,7,1,5]

				  Ter=Cstrain;
				  Tfr=Cstress;

				  if (Tstrain>=Tea)	{			// Rule 15 targeting for ea (ed) on 13
					  r13f(Ted,Teunn,fnewn,Enewn);
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
					  fca=fc;
					  Eta=Et;
					  esi=Ter;
					  fi=Tfr;
					  Ei=Ec;
					  esf=Tea;
					  RAf(esi,fi,Ei,esf,fca,Eta);
					  r15f(Ter,Tfr,Tea,fca,Eta,A,R);
					  Trule=15.0;
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain>=Teunn)	{	// Rule 13
					  r13f(Ted,Teunn,fnewn,Enewn);
					  Trule=13.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain>=esren)	{	// Rule 7
					  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
					  Trule=7.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else	{						// Rules 1 and 5
					  fcEtnf(Tstrain);					
				  }
			  } // if (Crule==14.0)

			  else if (Crule==3.0)	{

				  Ter0n=Cstrain;
				  Tfr0n=Cstress;
				  Tea=Ter0n;

				  fnewstnf(Tfunn,delfn,Teunn,Ter0n,espln);
				  Enewstnf(fnewstn,Tfr0n,Teunn,Ter0n);
				  esrestnf(Teunn,delen,Ter0n,espln);   
				  freErestnf(Teunn,Tfunn,Ter0n);

//				  opserr << "espln" << espln << "\n";
//				  opserr << "Tfr0n" << Tfr0n << "\n";

				  if (Tstrain>esrestn)	{		// Rule 77
					  r77f(Tstrain,Te0,Ter0n,Tfr0n,Teunn,fnewstn,Enewstn,esrestn,frestn,Erestn);
					  Trule=77.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else	{						// Rules 1 and 5						  
					  fcEtnf(Tstrain);					  
				  }
			  } // if (Crule==3.0)

			  else if (Crule==88.0)	{			// Reversal from transition 88 [Rules 11,(4,10 or 10),7,1,5 or Rules 4,10,7,1,5] 

				  if (Cstrain<=Teunp)	{			
					  Ter=Cstrain;
					  Tfr=Cstress;
					  Tea=Ter;
					  Teb=Ter0p;

					  if (Teb>=esplp)	{				// Reversal from 88 by Rules [11,4,10,7,1,5]	
					  	  
						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 4
							  r4f(Teunp,Tfunp,esplp,Eplp);	
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esplp)	{	// Rule 4
							  r4f(Teunp,Tfunp,esplp,Eplp);
							  Trule=4.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=Teunn)	{			// Reversal from 88 by Rules [11,10,7,1,5]

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Tea,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=Teunn)	{	// Rule 10
							  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
							  Trule=10.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }
					  }

					  else if (Teb>=esren)	{			// Reversal from 88 by Rules [11,7,1,5]

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 7
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Teb,esi,fi,esf,ff,Ei,Ef,A,R);
							  fcb=fc;
							  Etb=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain>=esren)	{	// Rule 7						  
							  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
							  Trule=7.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);
						  }			  
					  }

					  else	{		// (Teb<esren)		// Reversal from 88 by Rules [11,1,5]

						  if (Tstrain>=Teb)	{			// Rule 11 targeting for eb on 1 or 5
							  fcEtnf(Teb);
							  fcb=Tstress;
							  Etb=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Teb;
							  RAf(esi,fi,Ei,esf,fcb,Etb);
							  r11f(Ter,Tfr,Teb,fcb,Etb,A,R);
							  Trule=11.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 1 and 5
							  fcEtnf(Tstrain);					  
						  }
					  }				  
				  }	// if (Cstrain<=Teunp)
				  
				  else	{	// if (Cstrain>Teunp)		// Reversal from transition 88 by Rules [4,10,7,1,5] 

					  Teunp=Cstrain;				  
					  fcEtpf(Teunp,Te0);
					  Tfunp=Tstress;
					  esplpf(Teunp,Tfunp,Te0,espln);
					  Eplpf(Te0,Teunp);

					  if (Tstrain>=esplp)	{
						  r4f(Teunp,Tfunp,esplp,Eplp);
						  Trule=4.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=Teunn)	{
						  r10f(esplp,Eplp,Teunn,fnewn,Enewn);
						  Trule=10.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else if (Tstrain>=esren)	{
						  r7f(Teunn,fnewn,Enewn,esren,fren,Eren);
						  Trule=7.0;
						  RAf(esi,fi,Ei,esf,ff,Ef);
						  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  Tstress=fc;
						  Ttangent=Et;
					  }
					  else	{
						  fcEtnf(Tstrain);
					  }				  
				  }	// if (Cstrain<Teunn)			  
			  }	// if (Crule==77.0)
		  }	// if (Tstrain<Cstrain)

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

		  else	{	// if (Tstrain>=Cstrain)					// Continue going to positive direction	

			  if (Crule==3.0 || Crule==9.0 || Crule==8.0)	{	// Continue on transition 3 or 9 or 8 [Rules 3,9,8,2,6]

				  if (Tstrain<=espln)	{
					  r3f(Teunn,Tfunn,espln,Epln);
				  	  Trule=3.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
				  	  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain<=Teunp)	{
				  	  r9f(espln,Epln,Teunp,fnewp,Enewp);
				  	  Trule=9.0;
				  	  RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else if (Tstrain<=esrep)	{
				  	  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
				      Trule=8.0;
				      RAf(esi,fi,Ei,esf,ff,Ef);
				  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
				  	  Ttangent=Et;
			  	  }
			  	  else	{
				  	  fcEtpf(Tstrain,Te0);
			  	  }

			  }	// if (Crule==3)	or 9 or 8
			  
			  else if (Crule==2.0)	{		// Continue on transition 2 [Rules 2,6]
				  
				  fcEtpf(Tstrain,Te0);
			  }	// if (Crule==2)
			  
			  else if (Crule==6.0)	{		// Continue on transition 6 [Rule 6]
				  
				  fcEtpr6f(Tstrain,Te0);

			  }	// if (Crule==6)

			  else if (Crule==88.0)	{		// Continue on transition 88 [Rules 88,2,6]

				  if (Tstrain<esrestp)	{	// Rule 88							  
					  r88f(Tstrain,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
					  Trule=88.0;							  
					  RAf(esi,fi,Ei,esf,ff,Ef);							  
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);							  
					  Tstress=fc;							  
					  Ttangent=Et;						  
				  }						  
				  else	{							// Rules 2 and 6						  
					  fcEtpf(Tstrain,Te0);					  						  
				  }
			  }

			  else if (Crule==12.0)	{		// Continue on transition 12 [Rules 12,(3,9 or 9),8,2,6] or [Rules 12,88,2,6]	

				  if (Teb!=Ter0p)	{

					  if (Tea<=espln)	{

						  if (Tstrain<=Tea)	{	// Rule 12 targeting for ea on 3							  
							  r3f(Teunn,Tfunn,espln,Epln);   
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
						  	  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }
						  else if (Tstrain<=espln)	{	// Rule 3
							  r3f(Teunn,Tfunn,espln,Epln);
						  	  Trule=3.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9
						  	  r9f(espln,Epln,Teunp,fnewp,Enewp);
						  	  Trule=9.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;		  
						  }
					  	  else if (Tstrain<=esrep)	{	// Rule 8						  
						 	  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
						  	  Trule=8.0;
						  	  RAf(esi,fi,Ei,esf,ff,Ef);
						  	  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
						  	  Tstress=fc;
						  	  Ttangent=Et;
					  	  }		  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }

					  else if (Tea<=Teunp)	{			// (Tea>espln)

					  	  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp); 
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;			  
						  }
					  	  else if (Tstrain<=Teunp)	{	// Rule 9						  
							  r9f(espln,Epln,Teunp,fnewp,Enewp);
							  Trule=9.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain<=esrep)	{	// Rule 8						  
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{						// Rules 2 and 6			  
							  fcEtpf(Tstrain,Te0);		  
						  }
					  }

					  else if (Tea<=esrep)	{	

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 8
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  fca=fc;
							  Eta=Et;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else if (Tstrain<=esrep)	{	// Rule 8						  
							  r8f(Teunp,fnewp,Enewp,esrep,frep,Erep);
							  Trule=8.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }					  
						  else	{						// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }				  
					  }

					  else	{		// (Tea>esrep)

						  if (Tstrain<=Tea)	{			// Rule 12 targeting for ea on 2 or 6
							  fcEtpf(Tea,Te0);
							  fca=Tstress;
							  Eta=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{						// Rules 2 and 6
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }
				  }

				  else	{	// if (Teb==Ter0p)

					  if (Tea<=esrestp)	{

						  if (Tstrain<=Tea)	{				// Rule 12 targeting for ea on 88
							  r88f(Tea,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tea,esi,fi,esf,ff,Ei,Ef,A,R);
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else if (Tstrain<=esrestp)	{	// Rule 88
							  r88f(Tstrain,Te0,Ter0p,Tfr0p,Teunp,fnewstp,Enewstp,esrestp,frestp,Erestp);
							  Trule=88.0;
							  RAf(esi,fi,Ei,esf,ff,Ef);
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;
						  }
						  else	{							// Rules 2 and 6						  
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }

					  else	{	// if (Tea>esrestp)

						  if (Tstrain<=Tea)	{				// Rule 12 targeting for ea on 2 or 6
							  fcEtpf(Tea,Te0);
							  fca=Tstress;
							  Eta=Ttangent;
							  esi=Ter;
							  fi=Tfr;
							  Ei=Ec;
							  esf=Tea;
							  RAf(esi,fi,Ei,esf,fca,Eta);
							  r12f(Ter,Tfr,Tea,fca,Eta,A,R);
							  Trule=12.0;
							  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
							  Tstress=fc;
							  Ttangent=Et;				  
						  }					  
						  else	{							// Rules 2 and 6
							  fcEtpf(Tstrain,Te0);					  
						  }
					  }
				  }
			  } // if (Crule==12.0)

			  else if (Crule==14.0)	{						// Continue on transition 14 [Rules 14,66,6]

				  if (Tstrain<=Teb)	{						// Rule 14
					  r14f(Ter,Tfr,Teb);
					  Trule=14.0;
					  RAf(esi,fi,Ei,esf,ff,Ef);
					  fcEturf(Tstrain,esi,fi,esf,ff,Ei,Ef,A,R);
					  Tstress=fc;
					  Ttangent=Et;
				  }
				  else if (Tstrain<Teunp)	{				// Rule 66
					  r66f(Tstrain,Te0);
					  Trule=66.0;
				  }
				  else	{									// Rule 6
					  fcEtpr6f(Tstrain,Te0);
				  }
			  }

			  else if (Crule==66)	{						// Continue on transition 66 [Rules 66,6]

				  if (Tstrain<Teunp)	{
					  r66f(Tstrain,Te0);
					  Trule=66.0;
				  }
				  else	{									// Rule 6
					  fcEtpr6f(Tstrain,Te0);
				  }

			  }

		  }	// if (Tstrain>=Cstrain)

	  }	//	if (Cinc==1.0)

  }	// if monotonic or cyclic 

//Ttangent=6550.0;

return 0;	}


void Concrete05::fcEtnf(double e)
{
	xn=fabs(e/epcc);
	nn=fabs(Ec*epcc/fpcc);
	yf(xcrn,nn,rc);
	zf(xcrn,nn,rc);
	xsp=fabs(xcrn-y/(nn*z));

	if (xn<=xsp) {			//
		r1f(xn,nn,rc);
		Trule=1.0;
//		Ttangent=xsp;
	}
	else {
		r5f(xn,nn,rc);
		Trule=5.0;
//		Ttangent=xsp;
	}
}

void Concrete05::fcEtpf(double e, double e0)
{
	xp=fabs((e-e0)/et);
	np=Ec*et/ft;
	yf(xcrp,np,rt);
	zf(xcrp,np,rt);
	xcrk=fabs(xcrp-y/(np*z));

	if (xp<=xcrk) {			//
		r2f(xp,np,rt);
		Trule=2.0;
//		Ttangent=xcrk;
	}
	else {
		r6f(xp,np,rt);
		Trule=6.0;
//		Ttangent=xcrk;
	}
}

void Concrete05::fcEtpr6f(double e, double e0)
{
	xp=fabs((e-e0)/et);
	np=Ec*et/ft;
	r6f(xp,np,rt);
	Trule=6.0;
}

void Concrete05::yf(double x, double n, double r)
{
	double D;

	if (r!=1.0) {				//
		D=1.0+(n-r/(r-1.0))*x+pow(x,r)/(r-1.0);
	}
	else {
		D=1+(n-1+log10(x))*x;	//come back to log
	}
	
	y=n*x/D;
}

void Concrete05::zf(double x, double n, double r)
{
	double D;
	if (r!=1.0) {
		D=1.0+(n-r/(r-1.0))*x+pow(x,r)/(r-1.0);
	}
	else {
		D=1.0+(n-1.0+log10(x))*x;	//
	}
	
	z=(1-pow(x,r))/pow(D,2.0);
}

void Concrete05::esplnf(double eunn, double funn)
{
	Esecnf(eunn,funn);
	espln=eunn-funn/Esecn;
}
	  
void Concrete05::Eplnf(double eunn)
{
	Epln=0.1*Ec*exp(-2.0*fabs(eunn/epcc));
}	  

void Concrete05::Esecnf(double eunn, double funn)
{
	Esecn=Ec*((fabs(funn/(Ec*epcc))+0.57)/(fabs(eunn/epcc)+0.57));
}	  

void Concrete05::delenf(double eunn)
{
	delen=eunn/(1.15+2.75*fabs(eunn/epcc));
}	

void Concrete05::delfnf(double eunn, double funn)
{
	if (eunn<=epcc/10.0)	{
		delfn=0.09*funn*pow(fabs(eunn/epcc),0.5);  
//      delfn=0.0;
	}
	else	{
//      delfn=0.09*funn*pow(fabs(eunn/epcc),0.5);   
        delfn=0.0;
	}
}

void Concrete05::fnewnf(double eunn, double funn)
{
	delfnf(eunn,funn);
	fnewn=funn-delfn;
}

void Concrete05::Enewnf(double eunn, double funn)
{
	fnewnf(eunn,funn);
	esplnf(eunn,funn);
	Enewn = fnewn/(eunn-espln);

	if (Enewn < Ec)
	  Enewn=Ec;
}

void Concrete05::esrenf(double eunn)
{
	delenf(eunn);
	esren=eunn+delen;
}

void Concrete05::freErenf(double eunn)
{
	esrenf(eunn);
	
	xn=fabs(esren/epcc);
	nn=fabs(Ec*epcc/fpcc);
	yf(xcrn,nn,rc);
	zf(xcrn,nn,rc);
	xsp=fabs(xcrn-y/(nn*z));

	if (xn<=xsp) {			//
//		r1f(xn,nn,rc);
//		Trule=1.0;
		if (xn<xcrn) {
			yf(xn,nn,rc);
			zf(xn,nn,rc);
			fren=fpcc*y;
			Eren=Ec*z;
		}
		else {
			yf(xcrn,nn,rc);
			zf(xcrn,nn,rc);
			fren=fpcc*(y+nn*z*(xn-xcrn));
			Eren=Ec*z;
		}
	}
	else {
//		r5f(xn,nn,rc);
//		Trule=5.0;
		fren=0.0;
		Eren=0.0;
	}
//	fcEtnf(esren);
//	fren=Tstress;
//	Eren=Ttangent;
}

void Concrete05::fnewstnf(double funn, double delfn, double eunn, double er0n, double espln)
{
	fnewstn=funn-delfn*((eunn-er0n)/(eunn-espln));
}

void Concrete05::Enewstnf(double fnewstn, double fr0n, double eunn, double er0n)
{
	Enewstn=(fnewstn-fr0n)/(eunn-er0n);
}

void Concrete05::esrestnf(double eunn, double delen, double er0n, double espln)
{
	esrestn=eunn+delen*(eunn-er0n)/(eunn-espln);
}

void Concrete05::freErestnf(double eunn, double funn, double er0n)
{
	delenf(eunn);
	esplnf(eunn,funn);
	esrestnf(eunn,delen,er0n,espln);

	xn=fabs(esrestn/epcc);
	nn=fabs(Ec*epcc/fpcc);
	yf(xcrn,nn,rc);
	zf(xcrn,nn,rc);
	xsp=fabs(xcrn-y/(nn*z));

	if (xn<=xsp) {			//
//		r1f(xn,nn,rc);
//		Trule=1.0;
		if (xn<xcrn) {
			yf(xn,nn,rc);
			zf(xn,nn,rc);
			frestn=fpcc*y;
			Erestn=Ec*z;
		}
		else {
			yf(xcrn,nn,rc);
			zf(xcrn,nn,rc);
			frestn=fpcc*(y+nn*z*(xn-xcrn));
			Erestn=Ec*z;
		}
	}
	else {
//		r5f(xn,nn,rc);
//		Trule=5.0;
		frestn=0.0;
		Erestn=0.0;
	}

//	fcEtnf(esrestn);
//	frestn=Tstress;
//	Erestn=Ttangent;
}

void Concrete05::esplpf(double eunp, double funp, double e0, double espln)
{
	Esecpf(e0,eunp,funp,espln);
	esplp=eunp-funp/Esecp;
}
	  
void Concrete05::Eplpf(double e0, double eunp)
{
	Eplp=Ec/(pow((fabs((eunp-e0)/et)),1.1)+1.0);
//	Eplp=Eplp*0.5;
}	  

void Concrete05::Esecpf(double e0, double eunp, double funp, double espln)
{
	Esecp=Ec*((fabs(funp/(Ec*et))+0.67)/(fabs((eunp-e0)/et)+0.67));
//  Esecp=0.96*funp/(eunp-espln)+0.04*Esecp;

}	  

void Concrete05::delepf(double eunp, double e0)
{
	delep=0.22*fabs(eunp-e0);
//	delep=0.0;
}	

void Concrete05::delfpf(double funp, double eunp, double e0)
{
	if (eunp>=e0+et/2.0)	{
		delfp=0.15*funp;  
//      delfp=0.0;
	}
	else	{
//      delfp=0.15*funp;
        delfp=0.0;
	}
}

void Concrete05::fnewpf(double funp, double eunp, double e0)
{
	delfpf(funp,eunp,e0);
	fnewp=funp-delfp;
}

void Concrete05::Enewpf(double eunp, double funp, double e0, double espln)
{
	fnewpf(funp,eunp,e0);
	esplpf(eunp,funp,e0,espln);
	Enewp = fnewp/(eunp-esplp);
	if (Enewp <  Ec)
	  Enewp = Ec;	  
}

void Concrete05::esrepf(double eunp, double e0)
{
	delepf(eunp,e0);
	esrep=eunp+delep;
}

void Concrete05::freErepf(double eunp, double e0)
{
	esrepf(eunp,e0);

	xp=fabs((esrep-e0)/et);
	np=Ec*et/ft;
	yf(xcrp,np,rt);
	zf(xcrp,np,rt);
	xcrk=fabs(xcrp-y/(np*z));

	if (xp<=xcrk) {			//
//		r2f(xp,np,rt);
//		Trule=2.0;
		if (xp<xcrp) {
			yf(xp,np,rt);
			zf(xp,np,rt);
			frep=ft*y;
			Erep=Ec*z;
		}
		else {
			yf(xcrp,np,rt);
			zf(xcrp,np,rt);
			frep=ft*(y+np*z*(xp-xcrp));
//			zf(xcrp,np,rt);
			Erep=Ec*z;	
		}
	}
	else {
//		r6f(xp,np,rt);
//		Trule=6.0;
		frep=0.0;
		Erep=0.0;
	}
//	fcEtpf(esrep,e0);
//	frep=Tstress;
//	Erep=Ttangent;
}

void Concrete05::fnewstpf(double funp, double delfp, double eunp, double er0p, double esplp, double e0)
{
	fnewstp=funp-delfp*((eunp-er0p)/(eunp-esplp));
}

void Concrete05::Enewstpf(double fnewstp, double fr0p, double eunp, double er0p)
{
	Enewstp=(fnewstp-fr0p)/(eunp-er0p);
}

void Concrete05::esrestpf(double eunp, double delep, double er0p, double esplp)
{
	esrestp=eunp+delep*(eunp-er0p)/(eunp-esplp);
}

void Concrete05::freErestpf(double eunp, double funp, double er0p, double e0, double espln)
{
	delepf(eunp,e0);
	esplpf(eunp,funp,e0,espln);
	esrestpf(eunp,delep,er0p,esplp);

	xp=fabs((esrestp-e0)/et);
	np=Ec*et/ft;
	yf(xcrp,np,rt);
	zf(xcrp,np,rt);
	xcrk=fabs(xcrp-y/(np*z));
	
	if (xp<=xcrk) {			//
//		r2f(xp,np,rt);
//		Trule=2.0;
		if (xp<xcrp) {
			yf(xp,np,rt);
			zf(xp,np,rt);
			frestp=ft*y;
			Erestp=Ec*z;
		}
		else {
			yf(xcrp,np,rt);
			zf(xcrp,np,rt);
			frestp=ft*(y+np*z*(xp-xcrp));
//			zf(xcrp,np,rt);
			Erestp=Ec*z;	
		}
	}
	else {
//		r6f(xp,np,rt);
//		Trule=6.0;
		frestp=0.0;
		Erestp=0.0;
	}
	
//	fcEtpf(esrestp,e0);
//	frestp=Tstress;
//	Erestp=Ttangent;
}

void Concrete05::e0eunpfunpf(double e0,double eunp, double funp, double eunn, double funn)
{
	double xun=fabs(eunn/epcc);
	double xup=fabs((eunp-e0)/et);

	double e0ref;
	double eunpref;
	double funpref;

	if (xup<xun)	{
		xup=xun;
		e0ref=0.0;
		eunpref=xup*et;
		fcEtpf(eunpref,e0ref);
		funpref=Tstress;
	}
	else	{
		xup=xup;
		e0ref=e0;
		eunpref=eunp;
		funpref=funp;  
	}

	esplnf(eunn,funn);						//
	Eplnf(eunn);							//
	Esecpf(e0ref,eunpref,funpref,espln);	//

	double dele0=2.0*funpref/(Esecp+Epln);

	Te0=espln+dele0-xup*et;
	Teunp=xup*et+Te0;
	fcEtpf(Teunp,Te0);
	Tfunp=Tstress;

}


void Concrete05::r1f(double x, double n, double r)
{
	if (x<xcrn) {
		yf(x,n,r);
		zf(x,n,r);
		Tstress=fpcc*y;
		Ttangent=Ec*z;
	}
	else {
		yf(xcrn,n,r);
		zf(xcrn,n,r);
		Tstress=fpcc*(y+n*z*(x-xcrn));
		Ttangent=Ec*z;
	}
}


void Concrete05::r5f(double x, double n, double r)
{
	Tstress=0.0;
	Ttangent=0.0;
}


void Concrete05::r2f(double x, double n, double r)
{
	if (x<xcrp) {
		yf(x,n,r);
		zf(x,n,r);
		Tstress=ft*y;
		Ttangent=Ec*z;
	}
	else {
		yf(xcrp,n,r);
		zf(xcrp,n,r);
		Tstress=ft*(y+n*z*(x-xcrp));
//		zf(xcrp,n,r);
		Ttangent=Ec*z;
	}
}

void Concrete05::r6f(double x, double n, double r)
{
	Tstress=0.0;
	Ttangent=0.0;
}

void Concrete05::r3f(double eunn, double funn, double espln, double Epln)
{
	esi=eunn;
	fi=funn;
	Ei=Ec;
	esf=espln;
	ff=0.0;
	Ef=Epln;
}

void Concrete05::r9f(double espln, double Epln, double eunp, double fnewp, double Enewp)
{
	esi=espln;
	fi=0.0;
	Ei=Epln;
	esf=eunp;
	ff=fnewp;
	Ef=Enewp;
}

void Concrete05::r8f(double eunp, double fnewp, double Enewp, double esrep, double frep, double Erep)
{
	esi=eunp;
	fi=fnewp;
	Ei=Enewp;
	esf=esrep;
	ff=frep;
	Ef=Erep;
}

void Concrete05::r4f(double eunp, double funp, double esplp, double Eplp)
{
	esi=eunp;
	fi=funp;
	Ei=Ec;
	esf=esplp;
	ff=0.0;
	Ef=Eplp;
}

void Concrete05::r10f(double esplp, double Eplp, double eunn, double fnewn, double Enewn)
{
	esi=esplp;
	fi=0.0;
	Ei=Eplp;
	esf=eunn;
	ff=fnewn;
	Ef=Enewn;
}

void Concrete05::r7f(double eunn, double fnewn, double Enewn, double esren, double fren, double Eren)
{
	esi=eunn;
	fi=fnewn;
	Ei=Enewn;
	esf=esren;
	ff=fren;
	Ef=Eren;
}

void Concrete05::r12f(double er, double fr, double ea, double fca, double Eta, double A, double R)
{
	esi=er;
	fi=fr;
	Ei=Ec;
	esf=ea;
	ff=fca;
	Ef=Eta;
	fcEturf(ea,esi,fi,esf,ff,Ei,Ef,A,R);
	ff=fc;
	Ef=Et;
}

void Concrete05::r11f(double er, double fr, double eb, double fcb, double Etb, double A, double R)
{
	esi=er;
	fi=fr;
	Ei=Ec;
	esf=eb;
	ff=fcb;
	Ef=Etb;
	fcEturf(eb,esi,fi,esf,ff,Ei,Ef,A,R);
	ff=fc;
	Ef=Et;
}

void Concrete05::r13f(double ed, double eunn, double fnewn, double Enewn)
{
	esi=ed;
	fi=0.0;
	Ei=0.0;
	esf=eunn;
	ff=fnewn;
	Ef=Enewn;
}

void Concrete05::r14f(double er, double fr, double eb)
{
	esi=er;
	fi=fr;
	Ei=Ec;
	esf=eb;
	ff=0.0;
	Ef=0.0;
}

void Concrete05::r15f(double er, double fr, double ea, double fca, double Eta, double A, double R)
{
	esi=er;
	fi=fr;
	Ei=Ec;
	esf=ea;
	ff=fca;
	Ef=Eta;
	fcEturf(ea,esi,fi,esf,ff,Ei,Ef,A,R);
	ff=fc;
	Ef=Et;
}

void Concrete05::r66f(double e, double e0)
{
	Tstress=0.0;
	Ttangent=0.0;
}

void Concrete05::r88f(double e, double e0, double er0p, double fr0p, double eunp, double fnewstp, double Enewstp, double esrestp, double frestp, double Erestp)
{
	if ((e-e0)>=(er0p-e0) && (e-e0)<=(eunp-e0))	{
		esi=er0p;
		fi=fr0p;
		Ei=Ec;
		esf=eunp;
		ff=fnewstp;
		Ef=Enewstp;
	}
	if ((e-e0)>(eunp-e0) && (e-e0)<(esrestp-e0))	{
		esi=eunp;
		fi=fnewstp;
		Ei=Enewstp;
		esf=esrestp;
		ff=frestp;
		Ef=Erestp; 
	}
}

void Concrete05::r77f(double e, double e0, double er0n, double fr0n, double eunn, double fnewstn, double Enewstn, double esrestn, double frestn, double Erestn)
{

	if (e<=er0n && e>=eunn)	{
		esi=er0n;
		fi=fr0n;
		Ei=Ec;
		esf=eunn;
		ff=fnewstn;
		Ef=Enewstn;
	}
	if (e<eunn && e>esrestn)	{
		esi=eunn;
		fi=fnewstn;
		Ei=Enewstn;
		esf=esrestn;
		ff=frestn;
		Ef=Erestn; 
	}
}

void Concrete05::ea1112f(double eb, double espln, double esplp, double eunn, double eunp)
{
	Tea=espln+((eunn-eb)/(eunn-esplp))*(eunp-espln);
}

void Concrete05::eb1112f(double ea, double espln, double esplp, double eunn, double eunp)
{
	Teb=eunn-((ea-espln)/(eunp-espln))*(eunn-esplp);
}

void Concrete05::eb1415f(double ea, double fa, double Esecn)
{
	Teb=ea-fa/Esecn;
}

void Concrete05::RAf(double esi, double fi, double Ei, double esf, double ff, double Ef)
{
	double Esec=(ff-fi)/(esf-esi);
	R=(Ef-Esec)/(Esec-Ei);
	double check=pow(fabs(esf-esi),R);

	if (check==0.0 || check>1.797e308 || check<-1.797e308 || Esec==Ei)	{
		A=1.0e-300;
	}
	else	{
		A=(Esec-Ei)/pow(fabs(esf-esi),R);
		if (A>1.797e308 || A<-1.797e308)	{
			A=1.0e300;
		}
	}
}

void Concrete05::fcEturf(double es, double esi, double fi, double esf, double ff, double Ei, double Ef, double A, double R)
{
	double Esec=(ff-fi)/(esf-esi);
	
	if (A==1.0e300 || A==0.0)	{
		fc=fi+Esec*(es-esi);
		Et=Esec;
	}
	else if (pow(fabs(es-esi),-R)==0.0 || pow(fabs(es-esi),-R)>1.797e308 || pow(fabs(es-esi),-R)<-1.797e308)	{
		fc=fi+Esec*(es-esi);
		Et=Esec;
	}
	else if (Ei>=Esec && Ef>=Esec)	{
		fc=fi+Esec*(es-esi);
		Et=Esec;
	}
	else if (Ei<=Esec && Ef<=Esec)	{
		fc=fi+Esec*(es-esi);
		Et=Esec;
	}
	else	{
		fc=fi+(es-esi)*(Ei+A*pow(fabs(es-esi),R));
		Et=Ei+A*(R+1)*pow(fabs(es-esi),R);
		if (Et>1.797e308 || Et<-1.797e308)	{
			fc=fi+Esec*(es-esi);
			Et=Esec;
		}
	}
}




double Concrete05::getStress ()
{
   return Tstress;
}

double Concrete05::getStrain ()
{
   return Tstrain;
}

double Concrete05::getTangent ()
{
   return Ttangent;
}


int Concrete05::commitState ()
{
   // History variables
	Ceunn=Teunn;
	Cfunn=Tfunn;
	Ceunp=Teunp;
	Cfunp=Tfunp;
	Cer=Ter;
	Cfr=Tfr;
	Cer0n=Ter0n;
	Cfr0n=Tfr0n;
	Cer0p=Ter0p;
	Cfr0p=Tfr0p;
	Ce0=Te0;
	Cea=Tea;
	Ceb=Teb;
	Ced=Ted;
	Cinc=Tinc;
	Crule=Trule;
	
    // State variables
    Cstrain = Tstrain;
    Cstress = Tstress;
    Ctangent = Ttangent;	

	return 0;
}


int Concrete05::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
	Teunn=Ceunn;
	Tfunn=Cfunn;
	Teunp=Ceunp;
	Tfunp=Cfunp;
	Ter=Cer;
	Tfr=Cfr;
	Ter0n=Cer0n;
	Tfr0n=Cfr0n;
	Ter0p=Cer0p;
	Tfr0p=Cfr0p;
	Te0=Ce0;
	Tea=Cea;
	Teb=Ceb;
	Ted=Ced;
	Tinc=Cinc;
	Trule=Crule;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
//   Ttangent = 6550.0;
   
   return 0;
}


int Concrete05::revertToStart ()
{
	// Initial tangent

	double Ec0 = 6550.0; // initial stiffness compression

	// History variables

	Ceunn=0.0;
	Cfunn=0.0;
	Ceunp=0.0;
	Cfunp=0.0;
	Cer=0.0;
	Cfr=0.0;
	Cer0n=0.0;
	Cfr0n=0.0;
	Cer0p=0.0;
	Cfr0p=0.0;
	Ce0=0.0;
	Cea=0.0;
	Ceb=0.0;
	Ced=0.0;
	Cinc=0.0;
	Crule=0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}
	

UniaxialMaterial* Concrete05::getCopy ()
{
   Concrete05* theCopy = new Concrete05(this->getTag(),
									 fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp); 

   // Converged history variables
   theCopy-> Ceunn=Ceunn;
   theCopy-> Cfunn=Cfunn;
   theCopy-> Ceunp=Ceunp;
   theCopy-> Cfunp=Cfunp;
   theCopy-> Cer=Cer;
   theCopy-> Cfr=Cfr;
   theCopy-> Cer0n=Cer0n;
   theCopy-> Cfr0n=Cfr0n;
   theCopy-> Cer0p=Cer0p;
   theCopy-> Cfr0p=Cfr0p;
   theCopy-> Ce0=Ce0;
   theCopy-> Cea=Cea;
   theCopy-> Ceb=Ceb;
   theCopy-> Ced=Ced;
   theCopy-> Cinc=Cinc;
   theCopy-> Crule=Crule;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete05::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(29);	//come back to 29 (see below)
   data(0) = this->getTag();

   // Material properties
   data(1) = fpcc;
   data(2) = epcc;
   data(3) = Ec;
   data(4) = rc;
   data(5) = xcrn;
   data(6) = ft;
   data(7) = et;
   data(8) = rt;
   data(9) = xcrp;

   // History variables from last converged state

   data(10) = Ceunn;
   data(11) = Cfunn;
   data(12) = Ceunp;
   data(13) = Cfunp;
   data(14) = Cer;
   data(15) = Cfr;
   data(16) = Cer0n;
   data(17) = Cfr0n;
   data(18) = Cer0p;
   data(19) = Cfr0p;
   data(20) = Ce0;
   data(21) = Cea;
   data(22) = Ceb;
   data(23) = Ced;
   data(24) = Cinc;
   data(25) = Crule;

   // State variables from last converged state
   data(26) = Cstrain;
   data(27) = Cstress;
   data(28) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete05::sendSelf() - failed to send data\n";

   return res;
}

int Concrete05::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(29);		//come back to 29 (originally 24, data(23 max) above and below) )
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete05::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
	  fpcc = data(1);
	  epcc = data(2);
	  Ec = data(3); 
	  rc = data(4);
	  xcrn = data(5);
	  ft = data(6);
	  et = data(7);
	  rt = data(8);
	  xcrp = data(9);

      // History variables from last converged state
	  
	  Ceunn = data(10);
	  Cfunn = data(11);
	  Ceunp = data(12);
	  Cfunp = data(13);
	  Cer = data(14);
	  Cfr = data(15);
	  Cer0n = data(16);
	  Cfr0n = data(17);
	  Cer0p = data(18);
	  Cfr0p = data(19);
	  Ce0 = data(20);
	  Cea = data(21);
	  Ceb = data(22);
	  Ced = data(23);
	  Cinc = data(24);
	  Crule = data(25);

      // State variables from last converged state
      Cstrain = data(26);
      Cstress = data(27);
      Ctangent = data(28);


      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;		//come back to what this means
}

void Concrete05::Print (OPS_Stream& s, int flag)
{
   s << "Concrete05, tag: " << this->getTag() << endln;
   s << "  fpcc: " << fpcc << endln;
   s << "  epcc: " << epcc << endln;
   s << "  Ec: " << Ec << endln;
   s << "  rc: " << rc << endln;
   s << "  xcrn: " << xcrn << endln;
   s << "  ft: " << ft << endln;
   s << "  et: " << et << endln;
   s << "  rt: " << rt << endln;
   s << "  xcrp: " << xcrp << endln;

}
