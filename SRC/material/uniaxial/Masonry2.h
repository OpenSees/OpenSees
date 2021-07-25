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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-06-28 21:46:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Masonry.h,v $
                                                                        
// Written: Gonzalo Torrisi, UNCuyo
// Created: 10/2014
//
// Description: This file contains the class defination for Uniaxial material Masonry 
//			

#ifndef Masonry_h
#define Masonry_h

#include <UniaxialMaterial.h>


class Masonry : public UniaxialMaterial
{
	public:
		Masonry (int tag, double _Fm,     double _Ft,  
                     double _Um,    double _Uult,   double _Ucl,  
                   double _Emo,   double _Length, double _Area1,  
                   double _Area2, double _D1,     double _D2, 
                   double _Ach,   double _Are,    double _Ba,  
                   double _Bch,   double _Gun,    double _Gplu,   
                   double _Gplr,  double _Exp1,   double _Exp2,   
                   int _IENV);
		Masonry ();
		~Masonry();

		int setTrialStrain(double strain, double strainRate = 0.0); 
//      int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
		double getStrain(void);      
		double getStress(void);
		double getTangent(void);
		double getInitialTangent(void) ;

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
      double Fm;
	  double Ft;  
      double Um;
	  double Uult;
	  double Ucl;
	  double Emo;
	  double Length;
	  double Area1;
	  double Area2;
      double D1;
	  double D2;
	  double Ach;
	  double Are;
	  double Ba;  
      double Bch;
	  double Gun;
	  double Gplu;
	  double Gplr;
	  double Exp1;
	  double Exp2;
	  int IENV ;

      double  A1;
      double  A2; 
      double  Ko;
      double  DeltaU, Kfactor; 

	  double    K; 
double Ec;
    double Et; 
    double Uun; 
    double Sun; 
    double Eun ; 
    double Ure ; 
    double Sre ; 
    double Ere ; 
    double Uch ; 
    double Sch ; 
    double Ech ; 
    double U1 ; 
    double S1  ; 
    double  E1 ; 
    double U2  ; 
    double S2 ; 
    double E2 ; 
    double UunInt; 
    double UreInt ; 
    double Upl ; 
    double FtRed ; 
    int RuleNo ; 
    int IVIR ; 
    int InnerCycleNo ; 
    double Area ; 
     double D;
		 double U;
		 double F;
		 double S;
		 double UMAXIMA;
		 int INDIC;
		 double Uma;
		 double UMAXIMAT;
		 double UUTTT2;

   double cD  ;
   double cF  ;
   double cK  ;
   double cU  ;
   double cS   ;
   double cEt   ;
   double cUun  ;
   double cSun   ;
   double cEun ;
   double cUre  ;
   double cSre  ;
   double cEre  ;   
   double cUch   ;
   double cSch   ;
   double cEch  ;
   double cU1   ;
   double cS1   ;
   double cE1  ;
   double cU2  ;
   double cS2   ;
   double cE2    ;
   double cUunInt ;
   double cUreInt ;
   double cUpl  ;
   double cFtRed  ; 
   int cRuleNo ;
   int cIVIR   ;
   int cInnerCycleNo  ;
   double cArea   ;
   double cUMAXIMA;
   int cINDIC;



   void Comp_Envlp (double Uenv, double Um, double Fm, double Emo, double Uult, int IENV, double 
&Senv, double &Eenv ) ;

void Unload_Reload( double U, double U1, double  U2, double S1, double S2, double E1, double E2,
double &S, double &Et ) ;

void Plastic_Strain( double Uun, double Sun, double Um, double Fm, double Emo, double Ft, double 
Ba, double &Upl, double &FtRed  ) ;

void Stress_Tangent(double U, double DeltaU, double cU, double cS, double cEt,  
    double Um, double Fm, double Emo, double Ft, double Uult, double Ucl, double Ach, double Are, 
    double Ba, double Bch, double Gun, double Gplu, double Gplr, double Exp1, double Exp2,  
    double &U1, double &S1, double &E1, double &U2, double &S2, double &E2,double &S, double &Et,  
    double &FtRed, double &Upl, double &UunInt, double &UreInt, double &Uun, double &Sun, 
    double &Eun, double &Ure, double &Sre, double &Ere, double &Uch, double &Sch, double &Ech, 
    int &RuleNo, int &InnerCycleNo, int &IVIR, double &UMAXIMA, int &INDIC, double Uma);

void TRACCION(double Um, double &Upl, double Ft, double Emo, double &Et, double &S, double Uma, double U, double Ucl, double &UMAXIMA, int &INDIC);

};

#endif
