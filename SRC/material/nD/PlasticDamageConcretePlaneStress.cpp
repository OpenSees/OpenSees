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

 // Written in Matlab: Thanh Do
 // Created: 07/16

 #include <PlasticDamageConcretePlaneStress.h>           
 #include <Channel.h>
 #include <cmath>
 #include <elementAPI.h>

//#define _DEBUG_PDC_PlaneStress 1

 static Vector Iv6(6); 
 static Matrix Ivp(6,6); 
 static Matrix Idp(6,6); 
 static Matrix I(6,6);
 static Matrix Id(6,6); 

static const signed char b_A[3] = { -1, 1, 0 };
static const signed char c_a[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
static const signed char iv1[3] = { 0, 0, 1 };

void *
OPS_NewPlasticDamageConcretePlaneStress(void)
{
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5 || numArgs > 9) {
    opserr << "Want: nDMaterial PlasticDamageConcretePlaneStress $tag $E $nu $ft $fc <$beta $Ap $An $Bn>\n";
    return 0;
  }

  int iData[1];
  double dData[8];
  dData[4] = 0.6;
  dData[5] = 0.5;
  dData[6] = 2.0;
  dData[7] = 0.75;

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
    return 0;
  }

  numData = numArgs - 1;;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] <<"\n";
    return 0;
  }

  theMaterial = new PlasticDamageConcretePlaneStress(iData[0],
					    dData[0], dData[1], dData[2],dData[3],
					    dData[4], dData[5], dData[6],dData[7]);

  return theMaterial;
}

 PlasticDamageConcretePlaneStress::PlasticDamageConcretePlaneStress(int tag,
						  double _e, 
						  double _nu, 
						  double _ft,
						  double _fc, 
						  double _beta, 
						  double _Ap, 
						  double _An, 
						  double _Bn)
   :NDMaterial(tag,100),
    E(_e), nu(_nu), ft(_ft), fc(_fc), beta(_beta), Ap(_Ap), An(_An), Bn(_Bn),
    Ce(3,3), Ce0(3,3),CeCommitted(3,3),
    stress(3),strain(3), Cstress(3), Cstrain(3)
 {
   stress.Zero();
   strain.Zero();
   Cstress.Zero();
   Cstrain.Zero();
   
   double f2c = 1.16 * fc;
   double k = 1.4142135623730951 * (f2c - fc) / (2.0 * f2c - fc);

   //  initial damage threshold
   double rn0 = (-k + 1.4142135623730951) * fc / 1.7320508075688772;

   //  some useful constants
   //  shear modulus
   double G = E / 2.0 / (1.0 + nu);

   // initial tangent
   Ce0(0,0) = E / (1.0 - nu * nu);
   Ce0(0,1) = nu * E / (1.0 - nu * nu);
   Ce0(0,2) = 0.0;
   Ce0(1,0) = nu * E / (1.0 - nu * nu);
   Ce0(1,1) = E / (1.0 - nu * nu);
   Ce0(1,2) = 0.0;
   Ce0(2,0) = 0.0;
   Ce0(2,1) = 0.0;
   Ce0(2,2) = G;  

   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
       Ce(i,j) = Ce0(i,j);
   
   for (int i=0; i<3; i++) {
     sig[i]=0.;
     eps[i]=0.;
     Deps[i] = 0;
     eps_p[i] = 0;
   }
   eps_p[3] = 0;
   rn = rn0;
   rp = ft;
   dp = 0.;
   dn = 0.;

     this->commitState();
 }

 PlasticDamageConcretePlaneStress::PlasticDamageConcretePlaneStress()
   :NDMaterial (0, 0),
    stress(3),strain(3), Cstress(3), Cstrain(3)
 {

 }

 PlasticDamageConcretePlaneStress::~PlasticDamageConcretePlaneStress ()
 {

 }

 int
 PlasticDamageConcretePlaneStress::setTrialStrain(Vector const&v1, Vector const&v2){
   return this->setTrialStrain(v1);
 }

 int
 PlasticDamageConcretePlaneStress::setTrialStrain (const Vector &strain)
 {
   eps[0] = strain(0);
   eps[1] = strain(1);
   eps[2] = strain(2);

   double f2c = 1.16 * fc;
   double k = 1.4142135623730951 * (f2c - fc) / (2.0 * f2c - fc);

   //  initial damage threshold
   double rn0 = (-k + 1.4142135623730951) * fc / 1.7320508075688772;

   //  some useful constants
   //  shear modulus
   double G = E / 2.0 / (1.0 + nu);

   Ce = Ce0;

   // initial tangent
   double invCe[9];

   invCe[0] = 1.0 / E;
   invCe[3] = -nu / E;
   invCe[6] = 0.0;
   invCe[1] = -nu / E;
   invCe[4] = 1.0 / E;
   invCe[7] = 0.0;
   invCe[2] = 0.0;
   invCe[5] = 0.0;
   invCe[8] = 1.0 / G;

   //  tolerance for function sign evaluation
   //  tolerance for rotation angle
   //   double tol = 1.0e-14;
   //   double tolFi = 1.0e-14;
   
   //       %% state determination
   //  retrieve history variables from last convergence
   rn = Committed_rn;
   rp = Committed_rp;
   dp = Committed_dp;
   dn = Committed_dn;

   //  incremental strain
   //  trial effective stress
   //   double eps_p[3];
   eps_p[0] = Committed_eps_p[0];
   eps_p[1] = Committed_eps_p[1];
   eps_p[2] = Committed_eps_p[2];
   
   //   double b_MatState[3];
   
   for (int i0 = 0; i0 < 3; i0++) {
     Deps[i0] = eps[i0] - Committed_eps[i0];
   }

   double sige_tr[3];
   for (int i0 = 0; i0 < 3; i0++) {
     sige_tr[i0] = 0.0;
     for (int kstr = 0; kstr < 3; kstr++) {
       sige_tr[i0] += Ce(i0,kstr) * (eps[kstr] - eps_p[kstr]);
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector EPS(eps,3 );
   Vector EPS_P(eps_p,3 );
   Vector sige_TRIAL(sige_tr,3 );
   opserr << "rn: " << rn << "\nrp: " << rp << endln;
   opserr << "dp: " << dp << "\ndn: " << dn << endln;
   opserr << "EPS: " << EPS;
   opserr << "EPS_P: " << EPS_P;
   opserr << "sige_TRIAL: " << sige_TRIAL;
#endif

   double fi;
   //  rotation angle: principal effective stress directions
   if (fabs(sige_tr[0] - sige_tr[1]) > 1.0E-14) {
     fi = 0.5 * atan(2.0 * sige_tr[2] / (sige_tr[0] - sige_tr[1]));
   } else if (fabs(sige_tr[2]) > 1.0E-14) {
     fi = 0.78539816339744828;  // fi = pi/4.0
   } else {
     fi = 0.0;
   }
   
   if (fabs(fi) < 1.0E-14) {
     fi = 0.0;
   }
   
   //  principal effective stresses, sigp_effective
   double c = cos(fi);
   double s = sin(fi);
   double sigpe[2];
   sigpe[0] = (sige_tr[0] * (c * c) + sige_tr[1] * (s * s)) + 2.0 *
     sige_tr[2] * c * s;
   sigpe[1] = (sige_tr[0] * (s * s) + sige_tr[1] * (c * c)) - 2.0 *
     sige_tr[2] * c * s;
   
   //  negative and positive part of principal stresses
   double sigpeneg[2], sigpepos[2];
   for (int i0=0; i0<2; i0++) {
     sigpeneg[i0]=(sigpe[i0]-fabs(sigpe[i0]))/2.;
     sigpepos[i0]=(sigpe[i0]+fabs(sigpe[i0]))/2.;
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector sigPE(sigpe,2);
   opserr << "sigpe: " << sigPE;   
   Vector sigpeNEG(sigpeneg, 2);
   opserr << "sigpeneg: " << sigpeNEG;
   Vector sigpePOS(sigpepos, 2);
   opserr << "sigpepos: " << sigpePOS;
#endif

   //  octahedral stresses for negative part of effective stresses
   double sigoct = (sigpeneg[0] + sigpeneg[1])/3.0;
   double a = sigpeneg[0] - sigpeneg[1];
   double tauoct = sqrt(a*a + sigpeneg[1] * sigpeneg[1] + 
			sigpeneg[0] * sigpeneg[0]) / 3.0;


   //  negative equivalent stress
   double taunt = 1.7320508075688772 * (k * sigoct + tauoct);

   double taun = taunt;

   //  positive equivalent stress
   double taupt = sqrt((sigpepos[0] * sigpepos[0] + 
			sigpepos[1] * sigpepos[1]) - 2.0 *
		       nu * sigpepos[0] * sigpepos[1]);
   double taup = taupt;

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "sigoct: " << sigoct << " tauoct: " << tauoct << endln;
   opserr << "taunt: " << taunt << " taup: " << taup << endln;
   opserr << "rn: " << rn << endln;
#endif

   //  trial damage surface
   double gnegt = taunt - rn;

   //       %% Plasticity
   //  check if plastic strain coefficient is positive
   // THANH
   double nrm = 0.0;
   //   int L_tr_size_idx_0 = 1;
   double L_tr[3];
   double LDeps = 0.0;
   double gp = 0.0;
   double gn = 0.0;
   // THANH

   double lam = 0;
   double drp_dtaup = 0;
   double ddp_drp = 0;
   double drn_dtaun = 0;
   double ddn_drn = 0;
   //   double d_sige_tr[9];
   double dlam_dLDeps;
   double dlam_dnrm;
   //   double b_dnrm_ds;
   //   double dnrm_ds[3];

   double Deps_p[3];
   for (int i=0; i<3; i++) 
     Deps_p[i] = 0.0;


#ifdef _DEBUG_PDC_PlaneStress
     opserr << "gnegt: " << gnegt << endln;
#endif

   if (gnegt > 1.0E-14) {
     //  norm of the trial stress
     nrm = sqrt((sige_tr[0] * sige_tr[0] + sige_tr[1] * sige_tr[1]) + 2.0 *
		(sige_tr[2] * sige_tr[2]));

#ifdef _DEBUG_PDC_PlaneStress
     opserr << "NRM: " << nrm << endln;
#endif

     //  compute product of stress-strain, LDeps
     if (fabs(nrm) < 1.0E-14) {
       LDeps = 0;
       for (int i0 = 0; i0 < 3; i0++) {
	 L_tr[i0] = 0;
       }
     } else {
       //double L_tr_size_idx_0 = 3;
       for (int i0 = 0; i0 < 3; i0++) {
	 L_tr[i0] = sige_tr[i0] / nrm;
       }
       
       LDeps = 0;
       for (int i0 = 0; i0 < 3; i0++) {
	 LDeps += L_tr[i0] * Deps[i0];
       }
     }
#ifdef _DEBUG_PDC_PlaneStress
     Vector L_TR(L_tr,3);
     opserr << "L_TR: " << L_TR;
     Vector DEPS(Deps,3);
     opserr << "DEPS: " << DEPS;
     opserr << "LDEPS: " << LDeps << endln;
#endif

     //  check if the internal product is positive
     if (LDeps > 0.0) {
       //  compute ratio to update trial stress
       lam = 1.0 - (beta * E * LDeps) / nrm;

       //  update tension and compression equivalent stresses
       taun = taunt * lam;
       taup = taupt * lam;

       //  update damage function
       gp = taup - rp;
       gn = taun - rn;
       
       //  check if damage function is positive
       if ((gp > 1.0E-14) || (gn > 1.0E-14)) {
	 //  update plastic strain
	 for (int i0=0; i0< 3; i0++) {
	   // Deps_p   = beta*E*LDeps*invCe*L_tr;
	   Deps_p[i0] = beta*E*LDeps*(invCe[i0*3]*L_tr[0]+ invCe[i0*3+1]*L_tr[1]+ invCe[i0*3+2]*L_tr[2]);
	   // 
	   eps_p[i0] = eps_p[i0] + Deps_p[i0];
	 }

	 //  compute incremental plastic strain
       } else {
	 //  update variable without plastic strain
	 for (int i0=0; i0< 3; i0++) {
	   Deps_p[i0] = 0;
	   // eps_p = eps_p + Deps_p; does nothing if Deps_p == 0
	 }
	 lam = 0.0;
       }
            
     } else {
       //  update variable without plastic strain
       for (int i0=0; i0< 3; i0++) {
	 Deps_p[i0] = 0;
       }
       lam = 0.0;
     }
   } else {
     //  update variable without plastic strain
     for (int i0=0; i0< 3; i0++) {
       Deps_p[i0] = 0;
     }
     lam = 0.0;
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DEPS_P(Deps_p,3);
   opserr << "Deps_p: " << DEPS_P;
#endif

   //       %% Damage
   //  check if damage function is positive
   double gpos = taup - rp;
   double gneg = taun - rn;
   
   //  positive damage
   if (gpos > 1.0E-14) {
     //  update threshold
     //  update damage variable
     rp = taup; // rp0 = ft
     dp = 1.0 - (ft /rp) * exp(Ap * (1.0 - rp/ft));

     //  derivative terms
     drp_dtaup = 1;
     ddp_drp = (ft/(rp*rp) + Ap/rp) * exp(Ap * (1.0 - rp/ft));
#ifdef _DEBUG_PDC_PlaneStress
     opserr << "ft: " << ft << " Ap: " << Ap << endln;
     opserr << "rp: " << rp << " dp: " <<  dp << " ddp_drp: " << ddp_drp << endln;
#endif     
     //  limits for damage variable
     if (dp > 0.99999999999999) {
       dp = 0.99999999999999;
       ddp_drp = 0.0;
     }
   } else {
     //  no damage
     drp_dtaup = 0;
     ddp_drp = 0.0;
   }

   //  negative damage
   if (gneg > 1.0E-14) {
     //  update threshold
     //  update damage variable
     rn = taun;
     dn = (1.0 - sqrt(rn0/rn) * (1.0 - An)) - An * exp
       (Bn * (1.0 - sqrt(rn/rn0)));

     //  derivative terms
     drn_dtaun = 1;
     ddn_drn = (1.0 - An) * sqrt(rn0) / (2.0 * pow(rn,1.5))
       + An * Bn * exp(Bn * (1.0 - sqrt(rn/rn0)))
       / (2.0 * sqrt(rn*rn0));
     
     //  limits for damage variable
     if (dn > 0.99999999999999) {
       dn = 0.99999999999999;
       ddn_drn = 0.0;
     }
   } else {
     //  no damage
     drn_dtaun = 0;
     ddn_drn = 0.0;
   }

   //  Update stresses
   //  positive and negative effective principal stresses
   double Dpp[2];
   double Dnp[2];
   if (lam > 0.0) {
     for (int i0 = 0; i0 < 2; i0++) {
       Dpp[i0] = sigpepos[i0] * lam;
       Dnp[i0] = sigpeneg[i0] * lam;
     }
   } else {
     for (int i0 = 0; i0 < 2; i0++) {
       Dpp[i0] = sigpepos[i0];
       Dnp[i0] = sigpeneg[i0];
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DPP(Dpp,2);
   Vector DNP(Dnp,2);
   opserr << "dnp: " << DNP << " dpp: " << DPP;
   opserr << "dn: " << dn << " dp: " << dp << endln;
#endif


   //  principal stresses with damage
   double sigp[2];
   for (int i0 = 0; i0 < 2; i0++) {
     sigp[i0] = (1.0 - dn) * Dnp[i0] + (1.0 - dp) * Dpp[i0];
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector SIGP(sigp,2);
   opserr << "sigp: " << SIGP;   
#endif

   //  Cauchy stresses
   c = cos(fi);
   s = sin(fi);
   sig[0] = sigp[0]*(c*c) + sigp[1]*s*s;
   sig[1] = sigp[0]*(s*s) + sigp[1]*c*c;
   sig[2] = (sigp[0]-sigp[1])*(c*s);

#ifdef _DEBUG_PDC_PlaneStress
   Vector SIG(sig,3);
   opserr << "sig: " << SIG;
#endif

   //  Tangent
   //  dfi_ds := derivative of rotation angle fi w.r.t. the trial effective stresses 
   double dfi_ds[3];
   if (fabs(sige_tr[0] - sige_tr[1]) > 1.0E-14) {
     a = sige_tr[0] - sige_tr[1];
     double y = a * a + 4.0 * (sige_tr[2] * sige_tr[2]);
     dfi_ds[0] = -sige_tr[2]/y;
     dfi_ds[1] = sige_tr[2]/y;
     dfi_ds[2] = (sige_tr[0] - sige_tr[1])/y;
   } else if (fabs(sige_tr[2]) > 1.0E-14) {
     double y = 4.0 * sige_tr[2];
     dfi_ds[0]=-1.0/y;
     dfi_ds[1]=1.0/y;
     dfi_ds[2]=0.0;
   } else {
     for (int i0 = 0; i0 < 3; i0++) {
       dfi_ds[i0] = 0.0;
     }
   }


   //  no rotation
   if (fabs(fi) < 1.0E-14) {
     for (int i0 = 0; i0 < 3; i0++) {
       dfi_ds[i0] = 0.0;
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DFI_DS(dfi_ds,3);
   opserr << "dfi_ds: " << DFI_DS;
#endif
   
   //  dsigp_ds := derivative of principal effective stresses w.r.t. effective stresses 
   //      dsigp_ds = [(cos(fi))^2 (sin(fi))^2  2*cos(fi)*sin(fi)
   //                  (sin(fi))^2 (cos(fi))^2 -2*cos(fi)*sin(fi)];   

   //FMK
   static Matrix dsigp_ds(2,3); // 2X3
   c = cos(fi);
   s = sin(fi);
#ifdef _DEBUG_PDC_PlaneStress
   opserr << "FI: " << fi << endln;
#endif

   dsigp_ds(0,0) = c*c; dsigp_ds(0,1) = s*s; dsigp_ds(0,2) = 2*c*s;
   dsigp_ds(1,0) = s*s; dsigp_ds(1,1) = c*c; dsigp_ds(1,2) = -2*c*s;

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dsigp_ds: " << dsigp_ds;
#endif

   // dsigp_dfi = [2*(sige_tr(2)-sige_tr(1))*cos(fi)*sin(fi)+2*sige_tr(3)*((cos(fi))^2-(sin(fi))^2)
   //              2*(sige_tr(1)-sige_tr(2))*cos(fi)*sin(fi)+2*sige_tr(3)*((sin(fi))^2-(cos(fi))^2)];

#ifdef _DEBUG_PDC_PlaneStress
   Vector SIGE_TR(sige_tr,3);
   opserr << "sige_tr: " << SIGE_TR;
#endif

   static Vector dsigp_dfi(2); // 2X1
   dsigp_dfi(0) = 2*(sige_tr[1]-sige_tr[0])*c*s+2*sige_tr[2]*(c*c-s*s);
   dsigp_dfi(1) = 2*(sige_tr[0]-sige_tr[1])*c*s+2*sige_tr[2]*(s*s-c*c);

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dsigp_dfi: " << dsigp_dfi;
#endif

   // dsigp_ds = dsigp_ds + dsigp_dfi*dfi_ds;  2X3
   dsigp_ds(0,0) += dsigp_dfi(0)*dfi_ds[0];
   dsigp_ds(0,1) += dsigp_dfi(0)*dfi_ds[1];
   dsigp_ds(0,2) += dsigp_dfi(0)*dfi_ds[2];

   //   dsigp_ds = dsigp_ds + dsigp_dfi*dfi_ds;
   dsigp_ds(1,0) += dsigp_dfi(1)*dfi_ds[0];
   dsigp_ds(1,1) += dsigp_dfi(1)*dfi_ds[1];
   dsigp_ds(1,2) += dsigp_dfi(1)*dfi_ds[2];

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dsigp_ds: " << dsigp_ds;
#endif
   //  dsigpn_dD, dsigpp_dD := der. of positive and negative principal stresses w.r.t. effective stresses 
   static Matrix dsigpn_dsigp(2,2); dsigpn_dsigp.Zero();
   static Matrix dsigpp_dsigp(2,2); dsigpp_dsigp.Zero();

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "sigpe: " << sigPE;   
#endif

   if (sigpe[0] >= 0.0) {
     dsigpp_dsigp(0,0) = 1;
   } else {
     dsigpn_dsigp(0,0) = 1;
   }
   
   if (sigpe[1] >= 0.0) {
     dsigpp_dsigp(1,1) = 1;
   } else {
     dsigpn_dsigp(1,1) = 1;
   }

   //  dsigpn_ds, dsigpp_ds := der. of positive and negative principal stresses w.r.t. effective stresses     
   static Matrix dsigpn_ds(2,3);
   static Matrix dsigpp_ds(2,3);

   dsigpn_ds = dsigpn_dsigp*dsigp_ds;
   dsigpp_ds = dsigpp_dsigp*dsigp_ds;

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dsigpn_ds: " << dsigpn_ds;
   opserr << "dsigpp_ds: " << dsigpp_ds;
#endif

   // dsigoct_dsigpn = [1/3 1/3];
   // dsigoct_ds = dsigoct_dsigpn*dsigpn_ds;
   double dsigoct_ds[3];
   for (int i=0; i<3; i++) {
     dsigoct_ds[i] = .333333333*dsigpn_ds(0,i) + .333333333*dsigpn_ds(1,i);
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DSIGOCT_DS(dsigoct_ds,3);
   opserr << "DSIGOCT_DS" << DSIGOCT_DS;
#endif

   //  dsigoct_ds := der. of normal octahedral stresses w.r.t. effective stresses    
   double dtauoct_dsigpn[2];
   if (tauoct > 0.0) {
     //  dtauoct_dsigpn = [2*sigpeneg(1)-sigpeneg(2) 2*sigpeneg(2)-sigpeneg(1)]/(9*tauoct);      
     double y = 9.0 * tauoct;
     dtauoct_dsigpn[0] = (2.0 * sigpeneg[0] - sigpeneg[1])/y;
     dtauoct_dsigpn[1] = (2.0 * sigpeneg[1] - sigpeneg[0])/y;
   } else {
     for (int i0 = 0; i0 < 2; i0++) {
       dtauoct_dsigpn[i0] = 0.0;
     }
   }
   // dtauoct_ds = dtauoct_dsigpn*dsigpn_ds;
   double dtauoct_ds[3];
   for (int i=0; i<3; i++) {
     dtauoct_ds[i] = dtauoct_dsigpn[0]*dsigpn_ds(0,i) + 
       dtauoct_dsigpn[1]*dsigpn_ds(1,i);
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DTAUOCT_DS(dtauoct_ds,3);
   opserr << "DTAUOCT_DS" << DTAUOCT_DS;
#endif

   //  dtaupt_ds, dtaunt_ds := der. of pos. and neg. equivalent stresses w.r.t. effective stress 
   double dtaunt_dsigoct = 0.0;
   double dtaunt_dtauoct = 0.0;
   if (taunt > 0.0) {
     dtaunt_dsigoct = 1.7320508075688772 * k;
     dtaunt_dtauoct = 1.7320508075688772;
   } else {
     dtaunt_dsigoct = 0.0;
     dtaunt_dtauoct = 0.0;
   }

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "DTAUNT_DSIGOCT" << dtaunt_dsigoct << endln;
   opserr << "DTAUNT_DTAUOCT" << dtaunt_dtauoct << endln;
#endif

   // dtaunt_ds = dtaunt_dsigoct*dsigoct_ds + dtaunt_dtauoct*dtauoct_ds;
   double dtaunt_ds[3];
   double dtaun_ds[3];
   for (int i=0; i<3; i++) {
     dtaun_ds[i] = 0;
     dtaunt_ds[i] = dtaunt_dsigoct*dsigoct_ds[i] + dtaunt_dtauoct*dtauoct_ds[i];
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DTAUNT_DS(dtaunt_ds,3);
   opserr << "DTAUNT_DS" << DTAUNT_DS;
#endif


   double dtaupt_dsigpp[2];
   if (taupt > 0) {
     // dtaupt_dsigpp = [sigpepos(1)-nu*sigpepos(2) sigpepos(2)-nu*sigpepos(1)]/taupt;    
     dtaupt_dsigpp[0] = (sigpepos[0]-nu*sigpepos[1])/taupt;
     dtaupt_dsigpp[1] = (sigpepos[1]-nu*sigpepos[0])/taupt;
   } else {
     dtaupt_dsigpp[0] = 0;
     dtaupt_dsigpp[1] = 0;
   }

   Vector DTAUPT_DSIGPP(dtaupt_dsigpp,2);
   
   // dtaupt_ds = dtaupt_dsigpp*dsigpp_ds; (1x2) X (2x3)
   double dtaupt_ds[3];
   double dtaup_ds[3];
   for (int i=0; i<3; i++) {
     dtaup_ds[i] = 0;
     // dtaupt_ds = dtaupt_dsigpp*dsigpp_ds;
     dtaupt_ds[i] = dtaupt_dsigpp[0]*dsigpp_ds(0,i) + 
       dtaupt_dsigpp[1]*dsigpp_ds(1,i);
   }
#ifdef _DEBUG_PDC_PlaneStress
   Vector DTAUPT_DS(dtaupt_ds,3);
   opserr << "dtaupt_ds: " << DTAUPT_DS;
#endif
   //  PLASTICITY
   double dLDeps_ds[3];
   double dlam_ds[3];   
   for (int i=0; i<3; i++) {
     dLDeps_ds[i] = 0;
     dlam_ds[i] = 0;
   }

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "nrm: " << nrm << endln;
   opserr << "gnegt: " << gnegt << endln;
#endif

   if (gnegt > 1.0E-14) {
     
     //  THANH
     double y = 1.0 / nrm;

     // dnrm_ds = [sige_tr(1) sige_tr(2) 2*sige_tr(3)]*(1/nrm);
     double dnrm_ds[3];
     dnrm_ds[0] = sige_tr[0]*y;
     dnrm_ds[1] = sige_tr[1]*y;
     dnrm_ds[2] = 2*sige_tr[2]*y;

#ifdef _DEBUG_PDC_PlaneStress
     Vector SIGE_TR1(sige_tr,3);
     opserr << "SIGE_TR1: " << SIGE_TR1;
     Vector DNRM_DS(dnrm_ds,3);
     opserr << "DNRM_DS: " << DNRM_DS;
#endif

     if (fabs(nrm) >= 1.0E-14) {
       static Matrix dL_ds(3,3); 
       dL_ds.Zero(); 
       double *dDeps_ds = invCe; 

       //  dL_ds = I*(1/nrm) - sige_tr*dnrm_ds*(1/nrm^2);
       for (int i=0; i<3; i++) {
	   dL_ds(i,i) = 1/nrm;
	   for (int j=0; j<3; j++) {
	       dL_ds(i,j) = dL_ds(i,j)-(sige_tr[i]*dnrm_ds[j])/(nrm*nrm);
	     }
	 }

#ifdef _DEBUG_PDC_PlaneStress
       opserr << "DL_DS: " << dL_ds;
#endif

       double *dLDeps_dL = Deps;
       double *dLDeps_dDeps = L_tr;

       // dLDeps_ds = dLDeps_dL*dL_ds + dLDeps_dDeps*dDeps_ds;
       Matrix DNRM_DS(dnrm_ds,3,1);

       Matrix DLDEPS_dL(dLDeps_dL,1,3);
       Matrix DLDEPS_dDeps(dLDeps_dDeps,1,3);
       Matrix DDEPS_ds(dDeps_ds,3,3);
       Matrix DLDEPS_ds(dLDeps_ds,1,3);


       DLDEPS_ds = DLDEPS_dL*dL_ds + DLDEPS_dDeps*DDEPS_ds;
       /*
	 for (int i=0; i<3; i++) {
	 dLDeps_ds[i] = 0;
	 for (int j=0; j<3; j++) {
	   dLDeps_ds[i] += dLDeps_dL[i]*dL_ds(j,i) + dLDeps_dL[i]*dDeps_ds[i*3+j];
	 }
       }
       */
     }


     if (LDeps > 0.0) {
       // dlam_dLDeps = -beta*E/nrm;
       // dlam_dnrm = beta*E*LDeps/(nrm*nrm);
       // dlam_ds = dlam_dLDeps*dLDeps_ds + dlam_dnrm*dnrm_ds;
       // dtaun_ds = lam*dtaunt_ds + taun*dlam_ds/lam;
       // dtaup_ds = lam*dtaupt_ds + taup*dlam_ds/lam;

       dlam_dLDeps = -beta * E / nrm;
       dlam_dnrm = beta * E * LDeps / (nrm * nrm);
#ifdef _DEBUG_PDC_PlaneStress
       opserr << "dlam_dLDeps: " << dlam_dLDeps << endln;
       opserr << "dlam_dnrm: " << dlam_dnrm << endln;
       opserr << "DNRM_DS: " << DNRM_DS;
       Matrix DLDEPS_ds(dLDeps_ds,1,3);
       opserr << "DLDEPS_DS: " << DLDEPS_ds;
#endif
       for (int i0 = 0; i0 < 3; i0++) {
	 double b_dlam_ds = dlam_dLDeps * dLDeps_ds[i0] + dlam_dnrm * dnrm_ds[i0];
	 dlam_ds[i0] = b_dlam_ds;
	 dtaun_ds[i0] = lam * dtaunt_ds[i0] + taun * b_dlam_ds / lam;
	 dtaup_ds[i0] = lam * dtaupt_ds[i0] + taup * b_dlam_ds / lam;
       }
#ifdef _DEBUG_PDC_PlaneStress
       Vector DLAM_DS(dlam_ds,3);
       opserr << "DLAM_DS: " << DLAM_DS;
       Vector DTAUN_DS(dtaun_ds,3);
       opserr << "DTAUN_DS" << DTAUN_DS;
       Vector DTAUP_DS(dtaup_ds,3);
       opserr << "DTAUP_DS" << DTAUP_DS;
#endif
       if ((gp <= 1.0E-14) && (gn <= 1.0E-14)) {
	 for (int i0 = 0; i0 < 3; i0++) {
	   dlam_ds[i0] = 0.0;
	 }
       }
     } else {
       //  no plasticity, i.e. no adjustment with lambda
       for (int i0 = 0; i0 < 3; i0++) {
	 dlam_ds[i0] = 0.0;
	 dtaup_ds[i0] = dtaupt_ds[i0];
	 dtaun_ds[i0] = dtaunt_ds[i0];	 
       }
     }

   } else {
     //  no plasticity, i.e. no adjustment with lambda
     for (int i0 = 0; i0 < 3; i0++) {
       dlam_ds[i0] = 0.0;
       dtaup_ds[i0] = dtaupt_ds[i0];
       dtaun_ds[i0] = dtaunt_ds[i0];	 
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DTAUPT1_DS(dtaupt_ds,3);
   opserr << "DTAUPT_DS" << DTAUPT1_DS;
   opserr << "DTAUNT_DS" << DTAUNT_DS;

   Vector DTAUP_DS(dtaup_ds,3);
   Vector DTAUN_DS(dtaun_ds,3);
   opserr << "DTAUP_DS" << DTAUP_DS;
   opserr << "DTAUN_DS" << DTAUN_DS;

   opserr <<  "PlasticDamageConcretePlaneStress::sTS -12\n";     
   opserr << "gpos: " << gpos << endln;
#endif

   //  DAMAGE
   //  positive damage terms
   double ddp_ds[3];
   if (gpos > 1.0E-14) {
     for (int i0 = 0; i0 < 3; i0++) {
       ddp_ds[i0] = ddp_drp * ((double)drp_dtaup * dtaup_ds[i0]);
     }
   } else {
     for (int i0 = 0; i0 < 3; i0++) {
       ddp_ds[i0] = 0.0;
     }
   }
#ifdef _DEBUG_PDC_PlaneStress
   Vector DDP_DS(ddp_ds,3);
   opserr << "DDP_DS: " << DDP_DS;
#endif
   //  negative damage terms
   double ddn_ds[3];
   if (gneg > 1.0E-14) {
     for (int i0 = 0; i0 < 3; i0++) {
       ddn_ds[i0] = ddn_drn * ((double)drn_dtaun * dtaun_ds[i0]);
     }
   } else {
     for (int i0 = 0; i0 < 3; i0++) {
       ddn_ds[i0] = 0.0;
     }
   }
#ifdef _DEBUG_PDC_PlaneStress
   Vector DDN_DS(ddn_ds,3);
   opserr << "DDN_DS: " << DDN_DS;
#endif
   static Matrix dDpp_ds(2,3);
   static Matrix dDnp_ds(2,3);

   dDpp_ds.Zero();
   dDnp_ds.Zero();

   //FMK
#ifdef _DEBUG_PDC_PlaneStress
   opserr << "lam: " << lam << endln;
   Vector DLAM_DS(dlam_ds,3);
   opserr << "DLAM_DS: " << DLAM_DS;
#endif

   if (lam > 0.0) {
     // dDpp_ds = lam*dsigpp_ds + sigpepos'*dlam_ds;
     // dDnp_ds = lam*dsigpn_ds + sigpeneg'*dlam_ds;          
     for (int i0 = 0; i0 < 2; i0++) {
       for (int j0 = 0; j0 < 3; j0++) {
	 dDpp_ds(i0,j0)=lam*dsigpp_ds(i0,j0) + sigpepos[i0]*dlam_ds[j0];
	 dDnp_ds(i0,j0)=lam*dsigpn_ds(i0,j0) + sigpeneg[i0]*dlam_ds[j0];
       }
     }
   } else {
     // dDpp_ds = dsigpp_ds;      
     // dDnp_ds = dsigpn_ds;      
     for (int i0 = 0; i0 < 2; i0++) {
       for (int j0 = 0; j0 < 3; j0++) {
	 dDpp_ds(i0,j0) = dsigpp_ds(i0,j0);
	 dDnp_ds(i0,j0) = dsigpn_ds(i0,j0);
       }
     }
   }
#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dnp: " << DNP << " dpp: " << DPP;
   opserr << "dn: " << dn << " dp: " << dp << endln;
   opserr << "dDpp_ds: " << dDpp_ds;
   opserr << "dDnp_ds: " << dDnp_ds;
#endif
   // dsn_ds = -Dnp'*ddn_ds -Dpp'*ddp_ds + (1-dn)*dDnp_ds + (1-dp)*dDpp_ds;        
   Matrix dsn_ds(2,3);

   double dsig1_ds[3];
   double dsig2_ds[3];
   double dsig3_ds[3];

   for (int i = 0; i<2; i++) {
     for (int j=0; j<3; j++) {
       dsn_ds(i,j) = -Dnp[i]*ddn_ds[j] - Dpp[i] * ddp_ds[j] + (1-dn)*dDnp_ds(i,j) + (1-dp)*dDpp_ds(i,j);
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   opserr << "DDN_DS: " << DDN_DS;
   opserr << "DSN_DS: " << dsn_ds;
#endif
   

   c = cos(fi);
   s = sin(fi);
   for (int i=0; i<3; i++) {
     dsig1_ds[i] = c*c*dsn_ds(0,i) + s*s*dsn_ds(1,i) + 2*(sigp[1]-sigp[0])*s*c*dfi_ds[i];    
     dsig2_ds[i] = s*s*dsn_ds(0,i) + c*c*dsn_ds(1,i) + 2*(sigp[0]-sigp[1])*s*c*dfi_ds[i];    
   }
 

   if (fabs(sige_tr[2]) < 10e-14) {
     dsig3_ds[0] = 0;
     dsig3_ds[1] = 0;
     dsig3_ds[2] = 1;
   } else {
     // dsig3_ds=(sin(fi)*cos(fi))*[1 -1]*dsn_ds+ ...
     //             (sigp(1)-sigp(2))*((cos(fi))^2-(sin(fi))^2)*dfi_ds;    
     for (int i=0; i<3; i++) {

       dsig3_ds[i] = s*c*(dsn_ds(0,i)-dsn_ds(1,i)) + 
	 (sigp[0]-sigp[1])*(c*c-s*s)*dfi_ds[i];
     }
   }

#ifdef _DEBUG_PDC_PlaneStress
   Vector DSIG1_DS(dsig1_ds,3);
   Vector DSIG2_DS(dsig2_ds,3);
   Vector DSIG3_DS(dsig3_ds,3);

   opserr << "DSIG1_DS: " << DSIG1_DS;
   opserr << "DSIG2_DS: " << DSIG2_DS;
   opserr << "DSIG3_DS: " << DSIG3_DS;
#endif
    
   // compute tangent
   static Matrix dsig_ds(3,3);
   for (int i=0; i<3; i++) {
     dsig_ds(0,i) = dsig1_ds[i];
     dsig_ds(1,i) = dsig2_ds[i];
     dsig_ds(2,i) = dsig3_ds[i];
   }
#ifdef _DEBUG_PDC_PlaneStress
   opserr << "dsig_ds: " << dsig_ds;
#endif
   Ce = dsig_ds*Ce + 1.0e-12*Ce;
   for (int i = 0; i < 3; i++) {
     stress(i) = sig[i];
   }

   //   opserr <<  "PlasticDamageConcretePlaneStress::sTS -18\n";         
   return 0;
 }


int
PlasticDamageConcretePlaneStress::setTrialStrainIncr (const Vector &dStrain)
{
  strain += dStrain;
  this->setTrialStrain(strain);
  return 0;
}

int
PlasticDamageConcretePlaneStress::setTrialStrainIncr (const Vector &dStrain, const Vector &rate)
{
  strain += dStrain;
  this->setTrialStrain(strain);
  return 0;
}

const Matrix&
PlasticDamageConcretePlaneStress::getTangent (void)
{
  return Ce;
}

const Matrix&
PlasticDamageConcretePlaneStress::getInitialTangent (void)
{
  return Ce0;
}

const Vector&
PlasticDamageConcretePlaneStress::getStress (void)
{
  return stress;
}

const Vector&
PlasticDamageConcretePlaneStress::getStrain (void)
{
  return strain;
}

int
PlasticDamageConcretePlaneStress::commitState (void)
{
  CeCommitted = Ce;

  for (int i=0; i<4; i++) {
    Committed_eps_p[i] = eps_p[i];
  }
  for (int i=0; i<3; i++) {
    Committed_sig[i] = sig[i];
    Committed_eps[i] = eps[i];
  }

  Committed_rn = rn;
  Committed_rp = rp;
  Committed_dp = dp;
  Committed_dn = dn;

  Cstress = stress;
  Cstrain = strain;
  
  return 0;
}

int
PlasticDamageConcretePlaneStress::revertToLastCommit (void)
{
  Ce = CeCommitted; 

  for (int i=0; i<4; i++) {
    eps_p[i] = Committed_eps_p[i];
  }
  for (int i=0; i<3; i++) {
    sig[i] = Committed_sig[i];
    eps[i] = Committed_eps[i];
  }
  rn = Committed_rn;
  rp = Committed_rp;
  dp = Committed_dp;
  dn = Committed_dn;
  
  stress = Cstress;
  strain = Cstrain;

  return 0;
}

int
PlasticDamageConcretePlaneStress::revertToStart (void)
{
  return 0;
}

NDMaterial*
PlasticDamageConcretePlaneStress::getCopy(const char *type)
{
  if (strcmp(type,"PlaneStress") == 0 || strcmp(type,"PlaneStress2D") == 0) {
    PlasticDamageConcretePlaneStress *theCopy =
      new PlasticDamageConcretePlaneStress (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
    return theCopy;  
  } else {
    return 0;
  }
}

NDMaterial*
PlasticDamageConcretePlaneStress::getCopy (void)
{
  PlasticDamageConcretePlaneStress *theCopy =
    new PlasticDamageConcretePlaneStress (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
  return theCopy;
}

const char*
PlasticDamageConcretePlaneStress::getType (void) const
{
  return "PlaneStress2d";
}

int
PlasticDamageConcretePlaneStress::getOrder (void) const
{
  return 6;
}


int 
PlasticDamageConcretePlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcretePlaneStress::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PlasticDamageConcretePlaneStress::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcretePlaneStress::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
PlasticDamageConcretePlaneStress::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_CURRENTSTATE) {
    opserr << "PlasticDamageConcretePlaneStress: " << this->getTag();
    opserr << "strain: " << eps;
    opserr << "strain: " << sig;
    opserr << "tangent: " << this->getTangent();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"PlasticDamageConcretePlaneStress\", "; 
    s << " \"E\":" << E << ","; 
    s << " \"nu\":" << nu << ","; 
    s << " \"ft\":" << ft << ","; 
    s << " \"fc\":" << fc << ","; 
    s << " \"Ap\":" << Ap << ","; 
    s << " \"An\":"<< An << ","; 
    s << " \"Bn\":" << Bn << ","; 
    s << " \"beta\":" << beta << "}"; 
  }

}       
