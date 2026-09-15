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

// $Revision: 1.5 $
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel02M.cpp,v $

// Written: SK
// Updated: 1/5/26
// Pushed to repository 14/09/2026
// Description: 
#include <math.h>

#include <stdlib.h>
#include <Steel02M.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include <ID.h>

#include <tcl.h>          
#include <OPS_Stream.h>        
#include <UniaxialMaterial.h> 


int Steel02M::globalMaxIter = 0;

void*
OPS_Steel02M()
{
    UniaxialMaterial* theMaterial = 0;

    int    iData[1];
    double dData[14];
    int numData = 1;

    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial Steel02M tag" << endln;
        return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();

    if (numData < 8 || numData>14) {
        opserr << "WARNING invalid #args, want: uniaxialMaterial Steel02M "
            << iData[0] << " E0 Fy01 Fy02 alphaP alphaN R0 cR1 cR2 <a1 a2> <fysfy01 fysfy02> <b1 b2>\n";
        return nullptr;
    }

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid double inputs for uniaxialMaterial Steel02M "
            << iData[0] << endln;
        return nullptr;
    }

    double E0 = dData[0];
    double Fy01 = dData[1];
    double Fy02 = dData[2];
    double alphaP = dData[3];
    double alphaN = dData[4];
    double R0 = dData[5];
    double cR1 = dData[6];
    double cR2 = dData[7];

    double a1 = 0.0, a2 = 0.0;
    double Fysfy01 = 1.5, Fysfy02 = 1.5;
    double b1 = 0.8, b2 = 0.8;

    int idx = 8;

    // <a1 a2>
    if (numData >= idx + 2) {
        a1 = dData[idx];
        a2 = dData[idx + 1];
    }
    idx += 2;

    // <fysfy01 fysfy02>
    if (numData >= idx + 2) {
        Fysfy01 = dData[idx];
        Fysfy02 = dData[idx + 1];
    }
    idx += 2;

    // <b1 b2>
    if (numData >= idx + 2) {
        b1 = dData[idx];
        b2 = dData[idx + 1];
    }

    theMaterial = new Steel02M(iData[0],
        dData[0],  
        dData[1],  
        dData[2],  
        dData[3],  
        dData[4], 
        dData[5],  
        dData[6],  
        dData[7],  
        a1,  
        a2, 
        Fysfy01,  
        Fysfy02,  
        b1, 
        b2
    );

    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type Steel02M Material\n";
        return 0;
    }

    return theMaterial;
}


Steel02M::Steel02M(int tag, double _E0,
    double _Fy01, double _Fy02, double _alphaP,double _alphaN,
    double _R0, double _cR1, double _cR2,
    double _a1, double _a2,
    double _fysfy01, double _fysfy02,
    double _b1, double _b2) :
    UniaxialMaterial(tag, MAT_TAG_Steel02M),
    E0(_E0), R0(_R0),
    eps(0.0), sig(0.0), Et(_E0),
    epsP(0.0), sigP(0.0), eP(_E0)
{  
    alpha[0] = _alphaP;   alpha[1] = _alphaN;
    Fy0[0]   = _Fy01;     Fy0[1]   = _Fy02;
    cR[0]    = _cR1;      cR[1]    = _cR2;
    a[0]     = _a1;       a[1]     = _a2;
    Fysfy0[0]= _fysfy01;  Fysfy0[1]= _fysfy02;
    b[0]     = _b1;       b[1]     = _b2;
	
	Et = E0;
	eP = E0;

    epsy[0] = Fy0[0] / E0;
    epsy[1] = Fy0[1] / E0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
			if (i == 0) {
                alphaL[i][j] = alpha[j]; 
                alphaLP[i][j] = alpha[j];
                epsPL[i][j] = epsy[j];
                epsPL_P[i][j] = epsy[j];
                eps0L[i][j]  = epsy[j];
                eps0LP[i][j] = epsy[j];
                RL[i][j]     = R0;
                RLP[i][j]    = R0;
                sig0L[i][j]  = Fy0[j];
                sig0LP[i][j] = Fy0[j];   
            } else {
                alphaL[i][j] = 0.0;      
                alphaLP[i][j] = 0.0;
                epsPL[i][j] = 0.0;
                epsPL_P[i][j] = 0.0;
                eps0L[i][j]  = 0.0;
                eps0LP[i][j] = 0.0;
                RL[i][j]     = 0.0;
                RLP[i][j]    = 0.0;
                sig0L[i][j]  = 0.0;
                sig0LP[i][j] = 0.0;   
            }

			epsrL[i][j]  = 0.0;
            epsrLP[i][j] = 0.0;

            sigrL[i][j]  = 0.0;
            sigrLP[i][j] = 0.0;
        }
    }
	

    for (int i = 0; i < 3; i++) {
		Flag[i] = 0;
		FlagP[i] = 0;
    }

    epsr = 0.0;
    epsr_P = 0.0;

    sigr = 0.0;
    sigr_P = 0.0;

    
    for (int i = 0; i < 2; ++i) {
		Fy[i] = Fy0[i];              
		FyP[i] = Fy0[i];             
		Fys[i] = Fysfy0[i] * Fy0[i]; 
		FysP[i] = Fysfy0[i] * Fy0[i];
    }

}


Steel02M::Steel02M()
    : UniaxialMaterial(0, MAT_TAG_Steel02M),
    E0(0.0), R0(20.0),
    eps(0.0), sig(0.0), Et(0.0),
    epsP(0.0), sigP(0.0), eP(0.0),
    epsr(0.0), sigr(0.0), sig0(0.0),
    epsr_P(0.0), sigr_P(0.0), sig0P(0.0)
{
    for (int i = 0; i < 2; i++) {
        alpha[i] = 0.0;
        Fy0[i] = 0.0;
        cR[i] = 0.0;
        a[i] = 0.0;
        b[i] = 1.0;
        Fysfy0[i] = 0.0;
		epsy[i] = 0.0;
        Fy[i] = 0.0;
        FyP[i] = 0.0;
        Fys[i] = 0.0;
        FysP[i] = 0.0;
        
    }

    for (int i = 0; i < 3; i++) {
        Flag[i] = 0;
        FlagP[i] = 0;
        for (int j = 0; j < 2; j++) {
            alphaL[i][j]  = 0.0;
            alphaLP[i][j] = 0.0;

            epsPL[i][j] = 0.0;
            epsPL_P[i][j] = 0.0;

            eps0L[i][j]   = 0.0;
            eps0LP[i][j]  = 0.0;

            epsrL[i][j]   = 0.0;
            epsrLP[i][j]  = 0.0;

            sig0L[i][j]   = 0.0;
            sig0LP[i][j]  = 0.0;

            sigrL[i][j]   = 0.0;
            sigrLP[i][j]  = 0.0;

            RL[i][j]      = R0;
            RLP[i][j]     = R0;
        }
    }
}


Steel02M::~Steel02M(void)
{
    // Does nothing
}



UniaxialMaterial*
Steel02M::getCopy(void)
{
    Steel02M* theCopy = new Steel02M(getTag(),
        E0,
        Fy0[0], Fy0[1],
        alpha[0], alpha[1], 
        R0,
        cR[0], cR[1],
        a[0], a[1],
        Fysfy0[0], Fysfy0[1],
        b[0], b[1]);         

    return theCopy;
}

double
Steel02M::getInitialTangent(void)
{
    return E0;
}


int
Steel02M::setTrialStrain(double trialStrain, double strainRate)
{
    eps = trialStrain;
    double deps = eps - epsP;
	
	double R1t, R2t, R3t;
    double eps01t, sig01t;
    double eps02t, sig02t, alpha_2t;
    double eps03t, sig03t, alpha_3t;
    double xi=0.0;
    double Fy_trial;

    Et = eP;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
			
            alphaL[i][j] = alphaLP[i][j]; 
            epsPL[i][j] = epsPL_P[i][j];
 
            eps0L[i][j]  = eps0LP[i][j];
			epsrL[i][j]  = epsrLP[i][j];

            sig0L[i][j]  = sig0LP[i][j];                    
            sigrL[i][j]  = sigrLP[i][j];

            RL[i][j]     = RLP[i][j];
        }
    }
	

    for (int i = 0; i < 3; i++) {
		Flag[i] = FlagP[i];
    }
	

    epsr = epsr_P;
	sigr = sigr_P;

    
    for (int i = 0; i < 2; ++i) {
		Fy[i] = FyP[i];              
		Fys[i] = FysP[i];
    }
   
    if (Flag[0] == 0) {
            Flag[0] = (deps > 0) ? 1 : 2;
            if (deps == 0){
                sig=0.0;
                Et=E0;
            }
    }

    if ((Flag[0] == 2 && deps > 0) || (Flag[0] == 1 && deps < 0)) {

        sigr = sigP;
        epsr = eps - deps;
        double deps_y = (Fy0[0] - Fy0[1]) / E0;

        int x, y;
        if (deps > 0) {
            Flag[0] = 1;
            x = Flag[0];
            y = 2;
        } else {
            Flag[0] = 2;
            x = Flag[0];
            y = 1;
        }
		
		int idx = x - 1; 
        int idy = y - 1; 

        R1t = RL[0][idx], R2t = RL[1][idx], R3t = RL[2][idx];
        eps01t = eps0L[0][idx], sig01t = sig0L[0][idx];
        eps02t = eps0L[1][idx], sig02t = sig0L[1][idx], alpha_2t = alphaL[1][idx];
        eps03t = eps0L[2][idx], sig03t = sig0L[2][idx], alpha_3t = alphaL[2][idx];
        Fy_trial = Fy[idx];

        double Esh = alphaL[0][idx] * E0;
	
        if (Flag[1] == 0 || (x == 1 && epsr < epsPL[1][idy]) || (x == 2 && epsr > epsPL[1][idy])) 
            {
                epsPL[1][idy] = epsr;
            }   
        epsPL[2][idy] = epsr;
		
		Flag[1] = 0; 

        if (Flag[2] == 0) { 
            if (epsr > Fy0[0]/E0 || epsr < Fy0[1]/E0) {
                Flag[2] = 1; 
            } else {
                Flag[1] = 0; 
                
                std::tie(eps01t, sig01t) = Steel02M::Steel02M_Intersection(E0, epsr, Esh, Fy_trial, sigr);

                
                xi = 2.0 * fabs(epsPL[0][idx] - eps01t) / deps_y;
                R1t = R0 * (1.0 - (cR[0] * xi) / (cR[1] + xi));
            }
		}
		
		if (Flag[2] == 1) {

            if ((x == 1 && epsr < epsPL[0][idy]) || (x == 2 && epsr > epsPL[0][idy])) {
                epsPL[0][idy] = epsr;}

			double deps_m = epsPL[0][0] - epsPL[0][1];
			double zeta = deps_m / deps_y;
        
			double shft = 1.0 + a[idx] * pow(zeta, b[idx]);
            Fy_trial = shft * Fy0[idx];

            if (fabs(Fy_trial) > fabs(Fys[idx])) {
				Fy_trial = Fys[idx];
			}
        
            double SF = Fy_trial / Fy[idx]; 

           

			
			
			std::tie(eps01t, sig01t) = Steel02M::Steel02M_Intersection(E0, epsr, Esh, Fy_trial, sigr);
            
			xi = 2.0 * fabs(epsPL[0][idx] - eps01t) / deps_y;
			R1t = R0 * (1.0 - (cR[0] * xi) / (cR[1] + xi));
			
			double epsbar_r_1 = (epsr - epsrL[0][idx]) / (eps0L[0][idx] - epsrL[0][idx]);

			if (epsbar_r_1 > 0) {
				double sigbar_r_1 = (sigr - sigrL[0][idx]) / (sig0L[0][idx] - sigrL[0][idx]);
				double sigbar_L1_1_at_epsr = alphaL[0][idx] * epsbar_r_1 + 
                (1.0 - alphaL[0][idx]) * epsbar_r_1 / pow(pow(1.0 + pow(epsbar_r_1, RL[0][idx]), (1.0 / RL[0][idx])), 1.0);

				if (sigbar_L1_1_at_epsr > sigbar_r_1) { 
					double epsbar_01t_1 = (eps01t - epsrL[0][idx]) / (eps0L[0][idx] - epsrL[0][idx]);
					if (epsbar_01t_1 < 0) {
						opserr << "Steel02M::setTrialStrain - Warning: epsbar_01t_1 < 0. Wrong calculation.\n";
					}
					double sigbar_L1_1_at_eps01t = alphaL[0][idx] * epsbar_01t_1 + 
                    (1.0 - alphaL[0][idx]) * epsbar_01t_1 / pow(1.0 + pow(epsbar_01t_1, RL[0][idx]), (1.0 / RL[0][idx]));
					
					double sig_L1_at_eps01t = sigrL[0][idx] + sigbar_L1_1_at_eps01t * (sig0L[0][idx] - sigrL[0][idx]);

					double sigbar_L1_1t_at_eps01t = (sig_L1_at_eps01t - sigr) / (sig01t - sigr);

					double sigbar_L1t_1t_at_eps01t = alphaL[0][idx] + (1.0 - alphaL[0][idx]) / pow(2.0, (1.0 / R1t));

					if (sigbar_L1t_1t_at_eps01t > sigbar_L1_1t_at_eps01t) {
						Flag[1] = 1; 
					} 
					else {
						double epsbar_1m_Lim = (epsPL[0][idx] - epsrL[0][idx]) / (eps0L[0][idx] * SF - epsrL[0][idx]);
                    
						double sigbar_L1m_1m_at_epsLim = alphaL[0][idx] * epsbar_1m_Lim + 
                        (1.0 - alphaL[0][idx]) * epsbar_1m_Lim / pow(1.0 + pow(epsbar_1m_Lim, RL[0][idx]), (1.0 / RL[0][idx]));
                    
						double sig_L1m_at_epsLim = sigrL[0][idx] + sigbar_L1m_1m_at_epsLim * (sig0L[0][idx] * SF - sigrL[0][idx]);
                    
						double sigbar_L1m_1t_at_epsLim = (sig_L1m_at_epsLim - sigr) / (sig01t - sigr);
                    
						double epsbar_Lim_1t = (epsPL[0][idx] - epsr) / (eps01t - epsr);
                    
						double sigbar_L1t_1t_at_epsLim = alphaL[0][idx] * epsbar_Lim_1t + 
                        (1.0 - alphaL[0][idx]) * epsbar_Lim_1t / pow(1.0 + pow(epsbar_Lim_1t, R1t), (1.0 / R1t));

						if (sigbar_L1t_1t_at_epsLim > sigbar_L1m_1t_at_epsLim) {
							Flag[1] = 1;
						}
					}
					

					if (Flag[1] == 1) { 
						std::tie(alpha_2t, eps02t, sig02t) = Steel02M::Steel02M_IntersectionBasic(
						alphaL[0][idx], eps0L[0][idx], epsrL[0][idx], epsr,
						RL[0][idx], sig0L[0][idx], sigrL[0][idx], sigr, eps01t, sig01t);

						double xi2 = 2.0 * E0 * (epsPL[1][idx] - eps02t) / (Fy0[idx] - Fy0[idy]);


                        double xi = (xi2 > 0.0) ? xi2 : 0.0;

						R2t = R0 * (1.0 - (cR[0] * xi) / (cR[1] + xi));
					

					    double epsbar_r_2 = (epsr - epsrL[1][idx]) / (eps0L[1][idx] - epsrL[1][idx]);
					    double epsbar_r2_1 = (epsrL[1][idx] - epsrL[0][idx]) / (eps0L[0][idx] - epsrL[0][idx]);

					    if (epsbar_r2_1 > 0 && epsbar_r_2 > 0) { 

						    double sigbar_r_2 = (sigr - sigrL[1][idx]) / (sig0L[1][idx] - sigrL[1][idx]);

						    double sig_rL2;
						    std::tie(std::ignore, sig_rL2) = Steel02M::Steel02M_MCL1_Stress(
						    alphaL, epsr, eps0L, epsrL, RL, sig0L, sigrL, idx);

            
						    double sigbar_L2_2_at_epsr = (sig_rL2 - sigrL[1][idx]) / (sig0L[1][idx] - sigrL[1][idx]);

						    if (sigbar_L2_2_at_epsr > sigbar_r_2) {

							    double Et_at_eps02t, sig_L2_at_eps02t;
							    std::tie(Et_at_eps02t, sig_L2_at_eps02t) = Steel02M::Steel02M_MCL1_Stress(alphaL, eps02t, eps0L, epsrL, RL, sig0L, sigrL, idx);
							    double sigbar_L2_2_at_eps02t = (sig_L2_at_eps02t - sigr) / (sig02t - sigr);



							    double sigbar_L2t_2t_at_eps02t = alpha_2t + (1.0 - alpha_2t) / pow(2.0, 1.0 / R2t);

                
							    double fy_diff = Fy[0] - Fy[1];
							    if (sigbar_L2t_2t_at_eps02t > sigbar_L2_2_at_eps02t || 
                                    fabs(epsrL[1][idx] - epsrL[1][idy]) * E0 / fabs(fy_diff) < 1.0) {

								    std::tie(alpha_3t, eps03t, sig03t) = Steel02M::Steel02M_IntersectionMCL1(
                                    alphaL, E0, eps0L, epsrL, epsr, RL, sig0L, sigrL, sigr, idx, eps02t, sig02t, sig_L2_at_eps02t, Et_at_eps02t);

								    Flag[1] = 2;
								    double xi3 = 2.0 * E0 * (epsPL[2][idx] - eps03t) / (Fy0[idx] - Fy0[idy]);
								    xi = (xi3 > 0.0) ? xi3 : 0.0;

								    R3t = R0 * (1.0 - (cR[0] * xi) / (cR[1] + xi));
							    }
						    }
                        }
					}
				}
			}
		}
                
    int x_idx = Flag[0] - 1;

    if (Flag[1] == 0) {
        RL[0][x_idx]    = R1t;
        eps0L[0][x_idx] = eps01t;
        sig0L[0][x_idx] = sig01t;
        epsrL[0][x_idx] = epsr;
        sigrL[0][x_idx] = sigr;

        Fy[x_idx] = Fy_trial;

        eps0L[1][x_idx] = eps01t;
        sig0L[1][x_idx] = sig01t;
        epsrL[1][x_idx] = epsr;
        sigrL[1][x_idx] = sigr;
    } 
    else if (Flag[1] == 1) {
        RL[1][x_idx]    = R2t;
        eps0L[1][x_idx] = eps02t;
        sig0L[1][x_idx] = sig02t;
        epsrL[1][x_idx] = epsr;
        sigrL[1][x_idx] = sigr;
        alphaL[1][x_idx] = alpha_2t;
    } 
    else {
        RL[2][x_idx]    = R3t;
        eps0L[2][x_idx] = eps03t;
        sig0L[2][x_idx] = sig03t;
        epsrL[2][x_idx] = epsr;
        sigrL[2][x_idx] = sigr;
        alphaL[2][x_idx] = alpha_3t;
    }
	}

	int idx = Flag[0] - 1; 
    
	if (Flag[1] == 0) {
		std::tie(Et, sig) = Steel02M::Steel02M_GMP_Stress(alphaL[0][idx], eps, eps0L[0][idx], epsr, RL[0][idx], 
        sig0L[0][idx], sigr);
		} 
	else if (Flag[1] == 1) {
		std::tie(Et, sig) = Steel02M::Steel02M_MCL1_Stress(alphaL, eps, eps0L, epsrL, RL, sig0L, sigrL, idx);
		} 
	else {
		std::tie(Et, sig) = Steel02M::Steel02M_MCL2_Stress(alphaL, eps, eps0L, epsrL, RL, sig0L, sigrL, idx);
		}

    return 0;
}


std::tuple<double, double> Steel02M::Steel02M_Intersection(double E0, double epsr, double Esh,
 double fy,double sigr)
{
    double epsy = fy / E0;
    double eps0 = (fy - Esh * epsy - sigr + E0 * epsr) / (E0 - Esh);
    double sig0 = fy + Esh * (eps0 - epsy);

    return std::make_tuple(eps0, sig0);
}


std::tuple<double, double, double> 
Steel02M::Steel02M_IntersectionBasic(double alpha, double epsb0, double epsbr, 
                                     double epsr, double R, double sigb0, 
                                     double sigbr, double sigr, double eps_m, double sig_m)
{

    double sigr_dbar = (sigr - sigbr) / (sigb0 - sigbr);
    double epsr_dbar = (epsr - epsbr) / (epsb0 - epsbr);

    double eps_m_dbar = (eps_m - epsbr) / (epsb0 - epsbr);
    double sig_m_dbar = (sig_m - sigbr) / (sigb0 - sigbr);

    double term_m = pow(1.0 + pow(fabs(eps_m_dbar), R), (1.0 / R));
    double sig_mb_dbar = alpha * eps_m_dbar + (1.0 - alpha) * eps_m_dbar / term_m;
    double alpha_mb = alpha + (1.0 - alpha) / pow(1.0 + pow(fabs(eps_m_dbar), R), (1.0 + 1.0 / R));

    double dsig_m_dbar = sig_m_dbar - sig_mb_dbar;

    int maxIter = 10;
    int iter = 0;
    while (fabs(dsig_m_dbar) > 1.0e-5 && iter <= maxIter) {

        if (iter >= 10) {
            opserr << "WARNING: Steel02M_IntersectionBasic - Slow convergence at iter: " << iter 
                   << " dsig: " << fabs(dsig_m_dbar) << " eps_m_dbar: " << eps_m_dbar << endln;
        }
        
        double denom = (1.0 - alpha_mb);
        if (fabs(denom) < 1e-12) denom = (denom >= 0) ? 1e-12 : -1e-12;

        eps_m_dbar = (sig_mb_dbar - sigr_dbar + epsr_dbar - (alpha_mb * eps_m_dbar)) / (1.0 - alpha_mb);
        sig_m_dbar = sigr_dbar + (eps_m_dbar - epsr_dbar);

        term_m = pow(1.0 + pow(fabs(eps_m_dbar), R), (1.0 / R));
        sig_mb_dbar = alpha * eps_m_dbar + (1.0 - alpha) * eps_m_dbar / term_m;
        alpha_mb = alpha + (1.0 - alpha) / pow(1.0 + pow(fabs(eps_m_dbar), R), (1.0 + 1.0 / R));

        dsig_m_dbar = sig_m_dbar - sig_mb_dbar;
        iter++;
    }

    eps_m = epsbr + eps_m_dbar * (epsb0 - epsbr);
    sig_m = sigbr + sig_m_dbar * (sigb0 - sigbr);

    return std::make_tuple(alpha_mb, eps_m, sig_m);
}


std::tuple<double, double>
Steel02M::Steel02M_MCL1_Stress(double aL[3][2], double e, double e0L[3][2], 
                               double erL[3][2], double rL[3][2], double s0L[3][2], 
                               double srL[3][2], int idx)
{
    double Et, sig;
    double alpha_t, sig_hat;

    double dEps2 = e0L[1][idx] - erL[1][idx];
    double dSig2 = s0L[1][idx] - srL[1][idx];
    double dEps1 = e0L[0][idx] - erL[0][idx];
    double dSig1 = s0L[0][idx] - srL[0][idx];

    double eps_hat = (e - erL[1][idx]) / dEps2;

    if (eps_hat <= 1.0) {
        double term = 1.0 + pow(fabs(eps_hat), rL[1][idx]);
        
        sig_hat = aL[1][idx] * eps_hat + (1.0 - aL[1][idx]) * eps_hat / pow(term, 1.0 / rL[1][idx]);

        alpha_t = aL[1][idx] + (1.0 - aL[1][idx]) / pow(term, 1.0 + 1.0 / rL[1][idx]);
    } 
    else {

        double epsb_bar = (e - erL[0][idx]) / dEps1;

        
        double term1 = 1.0 + pow(fabs(epsb_bar), rL[0][idx]);
        double sigb_bar = aL[0][idx] * epsb_bar + (1.0 - aL[0][idx]) * epsb_bar / pow(term1, 1.0 / rL[0][idx]);
        double alpha_t1bb = aL[0][idx] + (1.0 - aL[0][idx]) / pow(term1, 1.0 + 1.0 / rL[0][idx]);

        double sigb = sigb_bar * dSig1 + srL[0][idx];

        double sigb_hat = (sigb - srL[1][idx]) / dSig2;

        double term2 = 1.0 + pow(fabs(eps_hat), rL[1][idx]);
        sig_hat = sigb_hat - (1.0 - aL[1][idx]) + (1.0 - aL[1][idx]) * eps_hat / pow(term2, 1.0 / rL[1][idx]);
        double alpha_t2 = (1.0 - aL[1][idx]) / pow(term2, 1.0 + 1.0 / rL[1][idx]);

        alpha_t = alpha_t1bb + alpha_t2;
    }

    sig = sig_hat * dSig2 + srL[1][idx];
    Et = alpha_t * (dSig2 / dEps2);

    return std::make_tuple(Et, sig);
}


std::tuple<double, double, double> 
Steel02M::Steel02M_IntersectionMCL1(double alphaL[3][2], double E0, double eps0L[3][2], 
                                    double epsrL[3][2], double epsr, 
                                    double RL[3][2], double sig0L[3][2], double sigrL[3][2], 
                                    double sigr, int idx,  double eps_n,
                                    double sig_n, double sig_nb, double Et)
{

    double dEps1 = eps0L[0][idx] - epsrL[0][idx];
    double dSig1 = sig0L[0][idx] - sigrL[0][idx];

    double sigr_dbar = (sigr - sigrL[0][idx]) / dSig1;
    double epsr_dbar = (epsr - epsrL[0][idx]) / dEps1;
    
    double eps_n_bbar = (eps_n - epsrL[0][idx]) / dEps1;
    double sig_n_bbar = (sig_n - sigrL[0][idx]) / dSig1;

    double sig_nb_dbar = (sig_nb - sigrL[0][idx]) / dSig1;
    double alpha_mn = Et / E0;

    double dsig_n_bbar = sig_n_bbar - sig_nb_dbar;
    int maxIter = 10;
    int iter = 0;

    while (fabs(dsig_n_bbar) > 1.0e-5 && iter <= maxIter) {
        eps_n_bbar = (sig_nb_dbar - sigr_dbar + epsr_dbar - (alpha_mn * eps_n_bbar)) / (1.0 - alpha_mn);
        sig_n_bbar = sigr_dbar + (eps_n_bbar - epsr_dbar);

        eps_n = epsrL[0][idx] + eps_n_bbar * dEps1;
        sig_n = sigrL[0][idx] + sig_n_bbar * dSig1;

        std::tie(Et, sig_nb) = Steel02M::Steel02M_MCL1_Stress(alphaL, eps_n, eps0L, epsrL, RL, sig0L, sigrL, idx);

        sig_nb_dbar = (sig_nb - sigrL[0][idx]) / dSig1;
        alpha_mn = Et / E0;

        dsig_n_bbar = sig_n_bbar - sig_nb_dbar;
        iter++;
    }

    eps_n = epsrL[0][idx] + eps_n_bbar * dEps1;
    sig_n = sigrL[0][idx] + sig_n_bbar * dSig1;

    return std::make_tuple(alpha_mn, eps_n, sig_n);
}


std::tuple<double, double>
Steel02M::Steel02M_GMP_Stress(double alpha, double eps, double eps0, double epsr, 
                              double R, double sig0, double sigr)
{
    double Et, sig;

    double deps = eps0 - epsr;
    double dsig = sig0 - sigr;

    double eps_star = (eps - epsr) / deps;

    double term1 = 1.0 + pow(fabs(eps_star), R);

    double sig_star = alpha * eps_star + (1.0 - alpha) * eps_star / pow(term1, 1.0 / R);

    sig = sig_star * dsig + sigr;

    double alpha_t = alpha + (1.0 - alpha) / pow(term1, 1.0 + 1.0 / R);

    Et = alpha_t * (dsig / deps);

    return std::make_tuple(Et, sig);
}




std::tuple<double, double>
Steel02M::Steel02M_MCL2_Stress(double aL[3][2], double e, double e0L[3][2], 
                               double erL[3][2], double rL[3][2], double s0L[3][2], 
                               double srL[3][2], int idx)
{
    double Et, sig;
    double alpha_t, sig_hat;

    double dEps3 = e0L[2][idx] - erL[2][idx];
    double dSig3 = s0L[2][idx] - srL[2][idx];

    double eps_hat = (e - erL[2][idx]) / dEps3;

    if (eps_hat <= 1.0) {
        double term = 1.0 + pow(fabs(eps_hat), rL[2][idx]);
        
        sig_hat = aL[2][idx] * eps_hat + (1.0 - aL[2][idx]) * eps_hat / pow(term, 1.0 / rL[2][idx]);
        alpha_t = aL[2][idx] + (1.0 - aL[2][idx]) / pow(term, 1.0 + 1.0 / rL[2][idx]);
    } 
    else {
        double Et1, sigb;
        
        std::tie(Et1, sigb) = Steel02M::Steel02M_MCL1_Stress(aL, e, e0L, erL, rL, s0L, srL, idx);

        double alpha_t1 = Et1 * (dEps3 / dSig3);
        double sigb_hat = (sigb - srL[2][idx]) / dSig3;

        double term = 1.0 + pow(fabs(eps_hat), rL[2][idx]);

        sig_hat = sigb_hat - (1.0 - aL[2][idx]) + (1.0 - aL[2][idx]) * eps_hat / pow(term, 1.0 / rL[2][idx]);
        double alpha_t2 = (1.0 - aL[2][idx]) / pow(term, 1.0 + 1.0 / rL[2][idx]);

        alpha_t = alpha_t1 + alpha_t2;
    }

    sig = sig_hat * dSig3 + srL[2][idx];
    Et = alpha_t * (dSig3 / dEps3);

    return std::make_tuple(Et, sig);
}


double
Steel02M::getStrain(void)
{
    return eps;
}

double
Steel02M::getStress(void)
{
    return sig;
}

double
Steel02M::getTangent(void)
{
    return Et;
}

int Steel02M::commitState(void)
{
    epsP = eps;
    sigP = sig;
    eP = Et;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            alphaLP[i][j] = alphaL[i][j];
            epsPL_P[i][j] = epsPL[i][j]; 
            eps0LP[i][j]  = eps0L[i][j];
            epsrLP[i][j]  = epsrL[i][j];
            sig0LP[i][j]  = sig0L[i][j];                        
            sigrLP[i][j]  = sigrL[i][j];
            RLP[i][j]     = RL[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        FlagP[i] = Flag[i];
    }
    
    epsr_P = epsr;
    sigr_P = sigr;

    for (int i = 0; i < 2; i++) {
        FyP[i]  = Fy[i];   
        FysP[i] = Fys[i]; 
    }

    return 0;
}

int Steel02M::revertToLastCommit(void)
{
    eps = epsP;
    sig = sigP;
    Et = eP;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            alphaL[i][j] = alphaLP[i][j];
            epsPL[i][j] = epsPL_P[i][j];
            eps0L[i][j]  = eps0LP[i][j];
            epsrL[i][j]  = epsrLP[i][j];
            sig0L[i][j]  = sig0LP[i][j];
            sigrL[i][j]  = sigrLP[i][j];
            RL[i][j]     = RLP[i][j];
        }
    }

    for (int i = 0; i < 3; i++) {
        Flag[i] = FlagP[i];
    }

    for (int i = 0; i < 2; i++) {
        Fy[i]   = FyP[i];  
        Fys[i]  = FysP[i]; 
    }

    epsr = epsr_P;
    sigr = sigr_P;

    return 0;
}

int Steel02M::revertToStart(void)
{
    epsP = 0.0;
    sigP = 0.0;
    eP   = E0;    
    epsr_P = 0.0;
    sigr_P = 0.0;

    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            alphaLP[i][j] = (i == 0) ? alpha[j] : 0.0;
            epsPL[i][j] = (i == 0) ? epsy[j] : 0.0;
            eps0LP[i][j]  = epsy[j];
            sig0LP[i][j]  = Fy0[j];
            epsrLP[i][j]  = 0.0;
            sigrLP[i][j]  = 0.0;
            RLP[i][j]     = R0;
        }
    }

    for (int i = 0; i < 3; i++) {
        FlagP[i] = 0;
    }

    for (int i = 0; i < 2; i++) {
        FyP[i]   = Fy0[i];
        FysP[i]  = Fysfy0[i] * Fy0[i];
    }

    this->revertToLastCommit();

    return 0;
}


int
Steel02M::sendSelf(int commitTag, Channel& theChannel)
{
    static ID idData(7); 
    static Vector dData(130); 

    idData(0) = this->getTag();
    idData(1) = Flag[0];
    idData(2) = Flag[1];
    idData(3) = Flag[2];
    idData(4) = FlagP[0];
    idData(5) = FlagP[1];
    idData(6) = FlagP[2];

    if (theChannel.sendID(this->getDbTag(), commitTag, idData) < 0) {
        opserr << "Steel02M::sendSelf -- could not send ID\n";
        return -1;
    }

    int pos = 0;

    dData(pos++) = E0;
    dData(pos++) = R0;
    for (int i = 0; i < 2; i++) dData(pos++) = alpha[i];
    for (int i = 0; i < 2; i++) dData(pos++) = Fy0[i];
    for (int i = 0; i < 2; i++) dData(pos++) = cR[i];
    for (int i = 0; i < 2; i++) dData(pos++) = a[i];
    for (int i = 0; i < 2; i++) dData(pos++) = b[i];
    for (int i = 0; i < 2; i++) dData(pos++) = Fysfy0[i];

  
    dData(pos++) = eps;    dData(pos++) = epsP;
    dData(pos++) = sig;    dData(pos++) = sigP;
    dData(pos++) = Et;     dData(pos++) = eP;
    dData(pos++) = epsr;   dData(pos++) = epsr_P;
    dData(pos++) = sigr;   dData(pos++) = sigr_P;
    dData(pos++) = sig0;   dData(pos++) = sig0P;

   
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            dData(pos++) = alphaL[i][j];
            dData(pos++) = epsPL[i][j];
            dData(pos++) = eps0L[i][j];
            dData(pos++) = epsrL[i][j];
            dData(pos++) = sig0L[i][j];
            dData(pos++) = sigrL[i][j];
            dData(pos++) = RL[i][j];
        }
    }

    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            dData(pos++) = alphaLP[i][j];
            dData(pos++) = epsPL_P[i][j];
            dData(pos++) = eps0LP[i][j];
            dData(pos++) = epsrLP[i][j];
            dData(pos++) = sig0LP[i][j];
            dData(pos++) = sigrLP[i][j];
            dData(pos++) = RLP[i][j];
        }
    }

    for (int i = 0; i < 2; i++) {
        dData(pos++) = Fy[i];
        dData(pos++) = FyP[i];
        dData(pos++) = Fys[i];
        dData(pos++) = FysP[i];
    }

    if (theChannel.sendVector(this->getDbTag(), commitTag, dData) < 0) {
        opserr << "Steel02M::sendSelf -- could not send Vector\n";
        return -2;
    }

    return 0;
}


int
Steel02M::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    static ID idData(7); 
    static Vector dData(130);

    if (theChannel.recvID(this->getDbTag(), commitTag, idData) < 0) {
        opserr << "Steel02M::recvSelf -- could not receive ID\n";
        return -1;
    }

    this->setTag(idData(0));
    Flag[0]  = idData(1);
    Flag[1]  = idData(2);
    Flag[2]  = idData(3);
    FlagP[0] = idData(4);
    FlagP[1] = idData(5);
    FlagP[2] = idData(6);

    if (theChannel.recvVector(this->getDbTag(), commitTag, dData) < 0) {
        opserr << "Steel02M::recvSelf -- could not receive Vector\n";
        return -2;
    }

    int pos = 0;

    E0 = dData(pos++);
    R0 = dData(pos++);
    for (int i = 0; i < 2; i++) alpha[i]    = dData(pos++);
    for (int i = 0; i < 2; i++) Fy0[i]      = dData(pos++);
    for (int i = 0; i < 2; i++) cR[i]       = dData(pos++);
    for (int i = 0; i < 2; i++) a[i]        = dData(pos++);
    for (int i = 0; i < 2; i++) b[i]        = dData(pos++);
    for (int i = 0; i < 2; i++) Fysfy0[i]   = dData(pos++);

    eps     = dData(pos++);    epsP    = dData(pos++);
    sig     = dData(pos++);    sigP    = dData(pos++);
    Et      = dData(pos++);    eP      = dData(pos++);
    epsr    = dData(pos++);    epsr_P  = dData(pos++);
    sigr    = dData(pos++);   sigr_P  = dData(pos++);
    sig0    = dData(pos++);    sig0P   = dData(pos++);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            alphaL[i][j] = dData(pos++);
            epsPL[i][j]=dData(pos++);
            eps0L[i][j]  = dData(pos++);
            epsrL[i][j]  = dData(pos++);
            sig0L[i][j]  = dData(pos++);
            sigrL[i][j]  = dData(pos++);
            RL[i][j]     = dData(pos++);
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            alphaLP[i][j] = dData(pos++);
            epsPL_P[i][j]=dData(pos++);
            eps0LP[i][j]  = dData(pos++);
            epsrLP[i][j]  = dData(pos++);
            sig0LP[i][j]  = dData(pos++);
            sigrLP[i][j]  = dData(pos++);
            RLP[i][j]     = dData(pos++);
        }
    }

    for (int i = 0; i < 2; i++) {
        Fy[i]    = dData(pos++);
        FyP[i]   = dData(pos++);
        Fys[i]   = dData(pos++);
        FysP[i]  = dData(pos++);
    }

    return 0;
}

void
Steel02M::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Steel02M tag: " << this->getTag() << endln;
        s << "  E0: " << E0 << ", R0: " << R0 << endln;
        
        s << "  alpha: [" << alpha[0] << ", " << alpha[1] << "]" << endln;
        
        s << "  Fy0: [" << Fy0[0] << ", " << Fy0[1] << "], ";
        s << "  cR: [" << cR[0] << ", " << cR[1] << "], ";
        
        s << "  fysfy0: [" << Fysfy0[0] << ", " << Fysfy0[1] << "], ";
        
        s << "  a: [" << a[0] << ", " << a[1] << "], ";
        s << "  b: [" << b[0] << ", " << b[1] << "]" << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"tag\": " << this->getTag() << ", ";
        s << "\"type\": \"Steel02M\", ";
        s << "\"E0\": " << E0 << ", ";
        
        s << "\"alpha\": [" << alpha[0] << ", " << alpha[1] << "], ";
        
        s << "\"R0\": " << R0 << ", ";
        s << "\"Fy0\": [" << Fy0[0] << ", " << Fy0[1] << "], ";
        s << "\"cR\": [" << cR[0] << ", " << cR[1] << "], ";
        
        s << "\"fysfy0\": [" << Fysfy0[0] << ", " << Fysfy0[1] << "], ";
        
        s << "\"a\": [" << a[0] << ", " << a[1] << "], ";
        s << "\"b\": [" << b[0] << ", " << b[1] << "]";
        s << "}";
    }
}

int
Steel02M::setParameter(const char** argv, int argc, Parameter& param)
{
    if (strcmp(argv[0], "Fy0") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(Fy0[idx]);
            return param.addObject(100 + idx, this);
        }
    }

    if (strcmp(argv[0], "alpha") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(alpha[idx]);
            return param.addObject(102 + idx, this);
        }
    }

    if (strcmp(argv[0], "b") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(b[idx]);
            return param.addObject(200 + idx, this);
        }
    }

    if (strcmp(argv[0], "a") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(a[idx]);
            return param.addObject(300 + idx, this);
        }
    }

    if (strcmp(argv[0], "R0") == 0) {
        param.setValue(R0);
        return param.addObject(400, this);
    }

    if (strcmp(argv[0], "cR") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(cR[idx]);
            return param.addObject(500 + idx, this);
        }
    }

    if (strcmp(argv[0], "fysfy0") == 0 && argc == 2) {
        int idx = atoi(argv[1]);
        if (idx >= 0 && idx < 2) {
            param.setValue(Fysfy0[idx]);
            return param.addObject(600 + idx, this);
        }
    }

    if (strcmp(argv[0], "E") == 0) {
        param.setValue(E0);
        return param.addObject(700, this);
    }

    return -1;
}

int
Steel02M::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case -1:
        return -1;

    case 100:
        Fy0[0] = info.theDouble;  
        epsy[0] = Fy0[0] / E0;
        break;
    case 101:
        Fy0[1] = info.theDouble;  
        epsy[1] = Fy0[1] / E0;
        break;

    case 102:
        alpha[0] = info.theDouble; 
        break;
    case 103:
        alpha[1] = info.theDouble; 
        
        break;

    case 200:
        b[0] = info.theDouble;
        break;
    case 201:
        b[1] = info.theDouble;
        break;
    case 300:
        a[0] = info.theDouble;
        break;
    case 301:
        a[1] = info.theDouble;
        break;
    case 400:
        R0 = info.theDouble;
        break;
    case 500:
        cR[0] = info.theDouble;
        break;
    case 501:
        cR[1] = info.theDouble;
        break;
    case 600:
        Fysfy0[0] = info.theDouble;
        break;
    case 601:
        Fysfy0[1] = info.theDouble;
        break;
    case 700:
        E0 = info.theDouble;
        epsy[0] = Fy0[0] / E0;
        epsy[1] = Fy0[1] / E0;
        break;

    default:
        return -1;
    }

    return 0;
}
