// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete05.h,v $
// Created: 07/04
// Written: KO
// Description: This file contains the class implementation for 
// Concrete04. Describe..
                                                                        
#ifndef Concrete05_h
#define Concrete05_h
#include <UniaxialMaterial.h>

class Concrete05 : public UniaxialMaterial
{
 public:
//Concrete05 (int tag, double fc, double eo, double r, double k, double alphaC, double fcr, double ecr, double b, double alphaT);  // input material parameters 
  Concrete05 (int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp);  // input material parameters
  Concrete05 ();
  ~Concrete05();

      int setTrialStrain(double strain, double strainRate = 0.0); 
//      int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
      double getStrain(void);      
      double getStress(void);
      double getTangent(void);

// come back here      double getInitialTangent(void) { return fc/eo*r/(r-1.0);} //initial in compression  change

	  
	  double getInitialTangent(void) {return Ec;}


// come back here	   double getec(void) {return eo;};

      int commitState(void);
      int revertToLastCommit(void);    
      int revertToStart(void);        

      UniaxialMaterial *getCopy(void);
    
      int sendSelf(int commitTag, Channel &theChannel);  
      int recvSelf(int commitTag, Channel &theChannel, 
	  	 FEM_ObjectBroker &theBroker);    
    
      void Print(OPS_Stream &s, int flag =0);

   protected:

   private:
      /*** Material Properties ***/

//	double ecr;	// input
//	double fcr;
//	double b;
//	double fc;
//	double eo;
//	double r;
//	double k;
//	double alphaC;
//	double alphaT;

	double fpcc;
	double epcc;
	double Ec;
	double rc;
	double xcrn;
	double ft;
	double et;
	double rt;
	double xcrp;


      /*** CONVERGED History Variables ***/
//      double CminStrain;   // Smallest previous concrete strain (compression)
//      double CunloadSlope; // Unloading (reloading) slope from CminStrain
//      double CendStrain;   // Strain at the end of unloading from CminStrain

// double Cecmax;	// output converged (previous)
// double Cet;
// double CetAccum;
// double Cscmax;
// double Cet1;
// double Cet2;
// double Cstmax;
// double Cetmax;
// double CEt;
// double CEr1;
// double CEr2;

	double Ceunn;
	double Cfunn;
	double Ceunp;
	double Cfunp;
	double Cer;
	double Cfr;
	double Cer0n;
	double Cfr0n;
	double Cer0p;
	double Cfr0p;
	double Ce0;
	double Cea;
	double Ceb;
	double Ced;
	double Cinc;
	double Crule;


      /*** CONVERGED State Variables ***/
      double Cstrain;
      double Cstress;   
      double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
						// Storing it is better than recomputing it!!!

      /*** TRIAL History Variables ***/
///     double TminStrain;
///     double TunloadSlope;
///     double TendStrain;

// double ecmax;
// double et;
// double etAccum;
// double scmax;
// double et1;
// double et2;
// double stmax;
// double etmax;
// double Et;
// double Er1;
// double Er2;

	double Teunn;
	double Tfunn;
	double Teunp;
	double Tfunp;
	double Ter;
	double Tfr;
	double Ter0n;
	double Tfr0n;
	double Ter0p;
	double Tfr0p;
	double Te0;
	double Tea;
	double Teb;
	double Ted;
	double Tinc;
	double Trule;


      /*** TRIAL State Variables ***/
      //double Tstrain;
      //double Tstress;
      //double Ttangent; // Not really a state variable, but declared here
                       // for convenience

//come back here//double Eci; 
				//double Eti;
				//double ecref;
				//double scref;
				//double etref;
				//double stref;
	double xn;
	double nn;
	double xsp;
	double xp;
	double np;
	double xcrk;
	double y;
	double z;

//	double x;
//	double n;
//	double r;

	double Tstrain;
    double Tstress;   
    double Ttangent;


     // void determineTrialState (double dStrain);

//     void reload();
//     void unload();
      
	double eunn;
	double funn;
	double espln;
	double Epln;
	double Esecn;
	double delen;
	double delfn;
	double fnewn;
	double Enewn;
	double esren;
	double fren;
	double Eren;
	double fnewstn;
	double Enewstn;
	double esrestn;
	double frestn;
	double Erestn;
//	double er0n;
//	double fr0n;
	
	double eunp;
	double funp;
	double esplp;
	double Eplp;
	double Esecp;
	double delep;
	double delfp;
	double fnewp;
	double Enewp;
	double esrep;
	double frep;
	double Erep;
	double fnewstp;
	double Enewstp;
	double esrestp;
	double frestp;
	double Erestp;
//	double er0p;
//	double fr0p;

	double esi;
	double fi;
	double Ei;
	double esf;
	double ff;
	double Ef;

	double R;
	double A;

	double fc;
	double Et;

	double fca;
	double Eta;

	double fcb;
	double Etb;

	double fa;
	double fb;


//	  esplpf(Teunp,Tfunp,Te0,espln);
//	  Eplpf(Te0,Teunp);
//	  Esecpf(Te0,Teunp,Tfunp,espln);  
//	  delepf(Teunp,Te0);
//	  delfpf(Tfunp,Teunp,Te0); 
//	  fnewpf(Tfunp,Teunp,Te0);
//	  Enewpf(Teunp,Tfunp,Te0,espln);
//	  esrepf(Teunp,Te0);
//	  freErepf(Teunp,Te0);
//	  fnewstpf(Tfunp,delfp,Teunp,Ter0p,esplp,Te0);
//	  Enewstpf(fnewstp,Tfr0p,Teunp,Ter0p);
//	  esrestpf(Teunp,delep,Ter0p,esplp);
//	  freErestpf(Teunp,Tfunp,Ter0p,Te0,espln);

void fcEtnf(double e);
void fcEtpf(double e, double e0);
void fcEtpr6f(double e, double e0);
void yf(double x, double n, double r);
void zf(double x, double n, double r);

void esplnf(double eunn, double funn);
void Eplnf(double eunn);
void Esecnf(double eunn, double funn);
void delenf(double eunn);
void delfnf(double eunn, double funn);
void fnewnf(double eunn, double funn);
void Enewnf(double eunn, double funn);
void esrenf(double eunn);
void freErenf(double eunn);
void fnewstnf(double funn, double delfn, double eunn, double er0n, double espln);
void Enewstnf(double fnewstn, double fr0n, double eunn, double er0n);
void esrestnf(double eunn, double delen, double er0n, double espln);   
void freErestnf(double eunn, double funn, double er0n);






void esplpf(double eunp, double funp, double e0, double espln);
void Eplpf(double e0, double eunp);
void Esecpf(double e0, double eunp, double funp, double espln);
void delepf(double eunp, double e0);
void delfpf(double funp, double eunp, double e0);
void fnewpf(double funp, double eunp, double e0);
void Enewpf(double eunp, double funp, double e0, double espln);
void esrepf(double eunp, double e0);
void freErepf(double eunp, double e0);
void fnewstpf(double funp, double delfp, double eunp, double er0p, double esplp, double e0);
void Enewstpf(double fnewstp, double fr0p, double eunp, double er0p);
void esrestpf(double eunp, double delep, double er0p, double esplp);   
void freErestpf(double eunp, double funp, double er0p, double e0, double espln);





void e0eunpfunpf (double e0, double eunp, double funp, double eunn, double funn);



void r1f(double x, double n, double r);
void r5f(double x, double n, double r);
void r2f(double x, double n, double r);
void r6f(double x, double n, double r);
void r3f(double eunn, double funn, double espln, double Epln);
void r9f(double espln, double Epln, double eunp, double fnewp, double Enewp);
void r8f(double eunp, double fnewp, double Enewp, double esrep, double frep, double Erep);
void r4f(double eunp, double funp, double esplp, double Eplp);
void r10f(double esplp, double Eplp, double eunn, double fnewn, double Enewn);
void r7f(double eunn, double fnewn, double Enewn, double esren, double fren, double Eren);
void r12f(double er, double fr, double ea, double fca, double Eta, double A, double R);
void r11f(double er, double fr, double eb, double fcb, double Etb, double A, double R);
void r13f(double ed, double eunn, double fnewn, double Enewn);
void r14f(double er, double fr, double eb);
void r15f(double er, double fr, double ea, double fca, double Eta, double A, double R);
void r66f(double e, double e0);
void r77f(double e, double e0, double er0n, double fr0n, double eunn, double fnewstn, double Enewstn, double esrestn, double frestn, double Erestn);
void r88f(double e, double e0, double er0p, double fr0p, double eunp, double fnewstp, double Enewstp, double esrestp, double frestp, double Erestp);

void ea1112f(double eb, double espln, double esplp, double eunn, double eunp);
void eb1112f(double ea, double espln, double esplp, double eunn, double eunp);
void eb1415f(double ea, double fa, double Esecn);


void RAf(double esi, double fi, double Ei, double esf, double ff, double Ef);
void fcEturf(double es, double esi, double fi, double esf, double ff, double Ei, double Ef, double A, double R);

	
//		void envelopeC(double e);  //compressive envelope (inside cpp file)
//	  void envelopeT(double e);  //tensile envelope
//	  void DefLoop(double Er);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};


#endif



//original version



