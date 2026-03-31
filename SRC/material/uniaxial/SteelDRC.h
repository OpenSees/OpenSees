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
                                                                        
// $Revision: 1.1 $
// $Date: 2008/12/09 20:00:16 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.h,v $
                                                                        
#ifndef SteelDRC_H_
#define SteelDRC_H_

// Written by: Rodrigo Carreno 
//
// Description: This file contains the class definition for the steel material originally
// developed by DOdd and Restrepo, with added corrections for efficiency, stability and 
// including post fracture behavior.
//
// What: "@(#) SteelDRC.h, revA" 
// Date: 2017/11/21

#include <UniaxialMaterial.h>

class SteelDRC : public UniaxialMaterial
{
  public:
    //SteelDRC(int tag, double Es, double fy, double eu, double fu, double esh, 
	//double Psh, double eft,double omegaFac = 1.0, int bauschType= 0, int stiffoption = 0);    
	SteelDRC(int tag, double Es, double fy, double eu, double fu, double esh,
		double Psh, double eft, double omegaFac, int bauschType, int stiffoption,
		double C_visc, double alpha, double Dfu);
	SteelDRC(int tag, double Es, double fy, double eu, double fu, double esh,
		double esh1, double fsh1, double eft, double omegaFac, int bauschType, 
		int stiffoption, double C_visc, double alpha, double Dfu);
    SteelDRC(int tag=0);
	//SteelDRC();
    ~SteelDRC();

    int setTrialStrain(double strain, double strainRate); 
    double getStrain(void);
	double getStrainRate(void);
    double getStress(void);
    double getTangent(void);

	const char *getClassType(void) const {return "SteelDRC";};
    double getInitialTangent(void);

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
	double fyEng; // Yield stress in engineering coordinates
	double eshEng; // Strain at onset of strain hardening in engineering coordinates
	double fuEng; // Ultimate strain in engineering coordinates
	// Material properties of the model (in Natural/True coordinates)
	double fyN; // Yield stress (in natural coordinates)
	double eyN; // Yield strain
	double E; // Youngs modulus (same in Natural and Engineering coordinates)
	double euN; // Ultimate strain 
	double fuN; // Ultimate stress
	double eshN; // Strain at onset of Strain hardening;
	double Psh;  // Exponent P defining the strain hardening branch of the backbone curve
	double eftN; // Post-necking strain at which stress becomes 0
	double omegaF; // Scaling factor to control area under bauschinger curve
	int bauschFlag; // Flag indicating formulation to be used to generate the bauschinger curve
					   // 0 for original iterative formulation
					   // 1 for closed-form with cubic rational Bezier curves (default)
					   // 2 for NURBS by Kim & Koutromanos (still needs fixing)
	int Etflag; // Flag to indicate stiffness to be returned by the material model
					// 0 for tangent stiffness
					// 1 for secant stiffness
					// 2 for original Young's modulus
	double C_visc; // Viscous damper constant applied to strain rate
	double alpha;	//Viscous damper exponent applied to strain rate
	double Dfu;		// Length of unloading branch in terms of fy
/* ========================================================================================================*/
	// State parameters in Committed state
	double Ceps; // Strain at committed state
	double Csig; // Stress at Ceps
	double Ctan; // Tangent at Ceps
	int Clmr; // Straining direction after last reversal from non-linear curve
	double Ce0[2]; // Strain shift of backbone curve in each direction
				   // Ce0[0] = Shift for tension curve 
				   // Ce0[1] = Shift for compression curve 
	double Ce0max; // maximum strain shift in absolute value
	double Cer;    // strain at last reversal from non-linear curve
	double Csr;    // stress at Cer;
	double CEr;    // slope at Cer;
	double Cea[2]; // Strain at the end of the unloading linear branch following 
				   // last non-linear reversals in each loading direction
				   // Cea[0] : End of unloading following reversal from compression
				   //		   to tension
				   // Cea[1]: End of unloading following reversal from tension to 
				   // 		  compression.
	double Csa[2]; // Stresses at Cea
	double Cerejoin[2]; //Strains where the stress-strain curve will rejoin shifted
						// skeleton curve in each loading direction before onset of 
						// strain hardening
						// Cerejoin[0] : rejoin backbone curve in tension
						// Cerejoin[1] : rejoin backbone curve in compression
						// NOTE: After onset of strain hardening Cerejoin[0]<0 and
						// 		 Cerejoin[1]>0
	double Csrejoin[2]; // Stress at Cerejoin
	double CErejoin[2]; // Slope at Cerejoin
	double CerejoinL[2];//Strains where the stress-strain curve will rejoin shifted skeleton
						// curve in each loading direction after the localization/uniform strain
						// point in either direction
						// CerejoinL[0] : rejoin backbone curve in tension
						// CerejoinL[1] : rejoin backbone curve in compression
						// NOTE: Before the localization strain is reached in tension or compression
						//       CerejoinL[0] = CerejoinL[1] = NaN
	double CsrejoinL[2];// Stress at CerejoinL
	double CErejoinL[2];// Slope at CerejoinL
	double Cerm[2]; 	// Strain at last major reversal and subsequent minor reversal
						// in the opposite direction.
						// Cerm[0] : Last major/minor reversal from tension to compression
						// Cerm[2] : Last major/minor reversal from compression to tension
	double Csrm[2]; 	// Stress at Cerm
	double CErm[2]; 	// Slope at Cerm
	int Cfrm[2]; 		// Flag of the reversal type in each loading direction
						// Cfrm[i] = 0 -> next reversal in direction i is not a minor reversal
						// Cfrm[i] = 1 -> next reversal in direction i is a minor reversal
						// Cfrm[i] = -1 -> next reversal in direction i is a major reversal
	double Ceam[2]; 	// Strain at the end of the unloading branches following the last major 
						// reversals in each loading direction.
						// Ceam[0] : End of unloading after major reversal from compression to
						// 			 tension
						// Ceam[1] : End of unloading after major reversal from tension to 
						//			 compression.
	double Csam[2]; 	// Stress at Ceam
	int Cfract;			// Flag indicating if material has past the ultimate strain in tension 
						// or compression
						// Cfract = 0 -> material has not yet exceeded the ultimate strain
						// Cfract = 1 -> material has started necking in tension
						// Cfract = 2 -> material has reached the ultimate strain in compression
	int CshOnset;		// Flag indicating if onset of strain hardening has occurred in the material yet 
						// CshOnset = 0 -> material has not reached onset of strain hardening yet (default)
						// CshOnset = 1 -> material has moved past the onset of strain hardening
						
/* ========================================================================================================*/
  // State parameters in Trial state
	double Teps; // Strain at trial state
	double Tsig; // Stress at Teps
	double Ttan; // Tangent at Teps
  
	int Tlmr; // Straining direction after last reversal from non-linear curve
	double Te0[2]; // Strain shift of backbone curve in each direction
				   // Te0[0] = Shift for tension curve 
				   // Te0[1] = Shift for compression curve 
	double Te0max; // maximum strain shift in absolute value
	double Ter;    // strain at last reversal from non-linear curve
	double Tsr;    // stress at Ter;
	double TEr;    // slope at Ter;
	double Tea[2]; // Strain at the end of the unloading linear branch following 
				   // last non-linear reversals in each loading direction
				   // Tea[0] : End of unloading following reversal from compression
				   //		   to tension
				   // Tea[1]: End of unloading following reversal from tension to 
				   // 		  compression.
	double Tsa[2]; // Stresses at Cea
	double Terejoin[2]; //Strains where the stress-strain curve will rejoin shifted
						// skeleton curve in each loading direction before onset of 
						// strain hardening
						// Terejoin[0] : rejoin backbone curve in tension
						// Terejoin[1] : rejoin backbone curve in compression
						// NOTE: After onset of strain hardening Terejoin[0]<0 and
						// 		 Terejoin[1]>0
	double Tsrejoin[2]; // Stress at Terejoin
	double TErejoin[2]; // Slope at Terejoin

	double TerejoinL[2];//Strains where the stress-strain curve will rejoin shifted skeleton
						// curve in each loading direction after the localization/uniform strain
						// point in either direction
						// TerejoinL[0] : rejoin backbone curve in tension
						// TerejoinL[1] : rejoin backbone curve in compression
						// NOTE: Before the localization strain is reached in tension or compression
						//       TerejoinL[0] = TerejoinL[1] = NaN

	double TsrejoinL[2];// Stress at TerejoinL
	double TErejoinL[2];// Slope at TerejoinL
	double Term[2]; 	// Strain at last major reversal and subsequent minor reversal
						// in the opposite direction.
						// Term[0] : Last major/minor reversal from tension to compression
						// Term[2] : Last major/minor reversal from compression to tension
	double Tsrm[2]; 	// Stress at Term
	double TErm[2]; 	// Slope at Term
	int Tfrm[2]; 		// Flag of the reversal type in each loading direction
						// Tfrm[i] = 0 -> next reversal in direction i is not a minor reversal
						// Tfrm[i] = 1 -> next reversal in direction i is a minor reversal
						// Tfrm[i] = -1 -> next reversal in direction i is a major reversal
	double Team[2]; 	// Strain at the end of the unloading branches following the last major 
						// reversals in each loading direction.
						// Team[0] : End of unloading after major reversal from compression to
						// 			 tension
						// Team[1] : End of unloading after major reversal from tension to 
						//			 compression.
	double Tsam[2]; 	// Stress at Team
	int Tfract;			// Flag indicating if material has past the ultimate strain in tension 
						// or compression
						// Tfract = 0 -> material has not yet exceeded the ultimate strain
						// Tfract = 1 -> material has started necking in tension
						// Tfract = 2 -> material has reached the ultimate strain in compression
	int TshOnset;		// Flag indicating if onset of strain hardening has occurred in the material yet 
						// TshOnset = 0 -> material has not reached onset of strain hardening yet (default)
						// TshOnset = 1 -> material has moved past the onset of strain hardening

/* =======================================================================================================*/
  
    double trialStrain;		 // trial strain
	double trialStrainRate;	 // trial strain rate
    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
    double commitStrain;     // last committed strain
	double commitStrainRate; // last committed strain rate
    double commitStress;     // last committed stress
    double commitTangent;    // last committed  tangent
	
	// Define internal functions used by the model
	void eng2natural(double *, int n=1);
	void natural2eng(double *, int n=1);
	void State_Determination(int S, int K, int M, int Klmr, double Eun);
	void bauschMinor(int flag, double *ptA, double *ptB, double omega,double eps_N, double &sig_N, double &tan_N);
	void bauschMajor(int flag, double *ptA, double *ptB, double eam, double sam, double e0, int S, int K, double eps_N, double &sig_N, double &tan_N);
	void bausch1 (double eps_N,double &sig_N,double &tan_N,double *pointA, double *pointB,double Pwr);
	void bauschNURBS(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB, double b);
	void bauschBezierCubic(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB, double *xi,double *wcubic);
	void nurbs(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB);
	void skeleton(double eps_N, double &sig_N, double &tan_N); // remaining parameters extracted from object
	double PowerP(double eam, double sam,double e0, int S, int K);
	double PowerP(double omega);
	//void tensionFracture(int S,int K,int M, int Klmr,double Eun);
	//void compressionFracture(int S, int K, int M, double Eun);
	void State_Reversal(int S, int K, int M, int &Klmr, double &Eun);
	//void MP_curve(double R0, double cR1, double cR2, double b, int S, int K, int M, double Eun, int flag0 = 0);
	//void MP_curve2(double R0, double cR1, double cR2, double b, double Eun, double e2, double f2,int M);
	double omegaFun(double eam, double sam, double e0, int S, int K);
	double factorb(double omega);
	double bezierWeightCubic(double omega);
};
#endif
