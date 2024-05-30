// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
//
// Created: 07/2015
//
// Description: This file contains the class implementation for 
// uniaxialMaterial ConcreteCM, which is a uniaxial hysteretic 
// constitutive model for concrete developed by Chang and Mander(1994).
// This model is a refined, rule - based, generalized, and non - dimensional 
// constitutive model that allows calibration of the monotonic and hysteretic 
// material modeling parameters, and can simulate the hysteretic behavior of 
// confined and unconfined, ordinary and high - strength concrete, in both 
// cyclic compression and tension. The model addresses important behavioral 
// features, such as continuous hysteretic behavior under cyclic compression 
// and tension, progressive stiffness degradation associated with smooth 
// unloading and reloading curves at increasing strain values, and gradual 
// crack closure effects. 
//
// References:
// 1) Chang, G.A. and Mander, J.B. (1994), “Seismic Energy Based Fatigue Damage 
// Analysis of Bridge Columns: Part I – Evaluation of Seismic Capacity”, NCEER 
// Technical Report No. NCEER-94-0006, State University of New York, Buffalo.
// 2) Kutay Orakcal (2004), "Nonlinear Modeling and Analysis of Slender Reinforced 
// Concrete Walls", PhD Dissertation, Department of Civil and Environmental Engineering, 
// University of California, Los Angeles.
//
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ConcreteCM.h
//
// Rev: 1


#ifndef ConcreteCM_h
#define ConcreteCM_h
#include <UniaxialMaterial.h>

class ConcreteCM : public UniaxialMaterial
{
public:
	// Constructors
	ConcreteCM(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp);
	ConcreteCM(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int mon); 
	ConcreteCM(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int Gap, int dummy); 
	ConcreteCM(int tag);
	ConcreteCM();

        const char *getClassType(void) const {return "ConcreteCM";};

	// Destructor
	~ConcreteCM();

	int setTrialStrain(double strain, double strainRate = 0.0); 
	double getStrain(void);      
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void) {return Ec;}

	Vector getInputParameters(void);

	int commitState(void);
	int revertToLastCommit(void);    
	int revertToStart(void);        

	double getCommittedStrain(void);
	double getCommittedStress(void);
	double getCommittedCyclicCrackingStrain(void);

	UniaxialMaterial *getCopy(void);

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		FEM_ObjectBroker &theBroker);    

	Response *setResponse (const char **argv, int argc, 
		OPS_Stream &theOutputStream);
	int getResponse (int responseID, Information &matInformation);   

	void Print(OPS_Stream &s, int flag =0);

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int    setParameter             (const char **argv, int argc, Information &info);
	int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
	// AddingSensitivity:END ///////////////////////////////////////////

protected:

private:

	// Input parameters
	double fpcc;	// peak compressive stress
	double epcc;	// strain at peak compressive stress	
	double Ec;		// Young's modulus
	double rc;		// shape parameter in compression
	double xcrn;	// cracking strain in compression
	double ft;		// peak tensile stress
	double et;		// strain at peak tensile stress
	double rt;		// shape parameter in tension
	double xcrp;	// cracking strain in tension
	int mon;		// switch for monotonic vs. cyclic (internal)
	int Gap;		// switch for more versus less gradual gap closure
	int dummy;		// dummy input variable

	// CONVERGED history variables
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

	// CONVERGED State Variables
	double Cstrain;
	double Cstress;   
	double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf

	// TRIAL history variables
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

	// TRIAL state variables
	double Tstrain;
	double Tstress;   
	double Ttangent;

	// Other variables
	double xn;
	double nn;
	double xsp;
	double xp;
	double np;
	double xcrk;
	double y;
	double z;

	double eunn;
	double funn;
	double espln;
	double Epln;
	double Esecn;
	double Esectest;
	double Esectest10; //KK
	double Esectest13; //KK
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

	double esi;
	double esi10; //KK
	double esi13; //KK
	double fi; 
	double fi10; //KK
	double fi13; //KK
	double Ei;
	double Ei10; //KK
	double Ei13; //KK
	double esf;
	double esf10; //KK
	double esf13; //KK
	double ff;
	double ff10; //KK
	double ff13; //KK
	double Ef;
	double Ef10; //KK
	double Ef13; //KK

	double R;
	double R10; //KK
	double R13; //KK
	double A;
	double A10; //KK
	double A13; //KK

	double fc;
	double Et;

	double fca;
	double Eta;

	double fcb;
	double Etb;

	double fa;
	double fb;

	// Member functions
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

	// AddingSensitivity:BEGIN //////////////////////////////////////////
	int parameterID;
	Matrix *SHVs;
	// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
