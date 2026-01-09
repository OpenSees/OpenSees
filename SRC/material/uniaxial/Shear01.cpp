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
//
// $Date: 2014/04/13 12:26
// 
// Written: Erkan Bicici 
//
// Description: The file contains the application of shear model 
//              developed by Sezen (2002). The material is defined 
//              with sezen point.  
//
//		
//


#include <elementAPI.h>
#include "Shear01.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>



static int numShear01 = 0;

void *
OPS_Shear01()
{
  // print out some KUDO's
  if (numShear01 == 0) {
    opserr << "Shear01 unaxial material - Written by Erkan 2014\n";
    numShear01 =1;
  }
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[12];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Shear01 tag" << endln;
    return 0;
  }
  numData = 12;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Parameters\n";
    return 0;	
  }
  // 
  // create a new material
  //
  theMaterial = new Shear01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);       
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of Shear01\n";
    return 0;
  }
  // return the material
  return theMaterial;
}

Shear01::Shear01(int tag, double ep, double crkp, double e2p, double eypp, double eyep, double eafp,  double en, double crkn, double e2n, double eypn, double eyen, double eafn)
:UniaxialMaterial(tag, 0),
 ezero(0.0), Ep(ep), ecrp(crkp), E2p(e2p), epp(eypp), eyp(eyep), eup(eafp),En(en), ecrn(crkn), E2n(e2n), epn(eypn), eyn(eyen), eun(eafn),
 commitStrain(0.0), commitStress(0.0), commitTangent(Ep), Old_Dif_Strain(0.0)
{
  	// Set envelope slopes
	this->setParameters();
	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

Shear01::~Shear01()
{
  // does nothing
}

int 
Shear01::setTrialStrain(double strain, double strainRate)
{

 trialStrain = strain;//recording trialstrain
 Dif_Strain = trialStrain - commitStrain; // difference from the previous one
 Cur_Dif_Strain = Dif_Strain; // Current difference. The program stores the previous and current differences
 NumberCycle();// The number of cyclec calculation

// Store Maximum and minimum
 if (commitStress > CmaxStress)		{TmaxStress = commitStress;}
 if (commitStress < CminStress)		{TminStress = commitStress;} 

 // Calculation of stress
 
	if (trialStrain == 0.0){ trialTangent = 0.0; trialStress  = 0.0; }
		
	if ((trialStrain >= CmaxStrain)&&(TEnvIndicator == 0)) {
	TmaxStrain = trialStrain;
	TenvmaxStress = PosEnvStress(TmaxStrain);
				//if (TEnvIndicator == 0) {
					trialTangent = PosEnvTangent(trialStrain);
					trialStress  = PosEnvStress(trialStrain);
				/*}
				else {
					reloading(Dif_Strain);
				}*/

		}
	else if ((trialStrain <= CminStrain) && (TEnvIndicator == 0)) {
	TminStrain = trialStrain;
	TenvminStress = NegEnvStress(TminStrain);

				//if (TEnvIndicator == 0) {
					trialTangent = NegEnvTangent(trialStrain);
					trialStress  = NegEnvStress(trialStrain);
					/*					}
				else {
					unloading(Dif_Strain);
					}*/


		}
	else{ 	
				if (Dif_Strain < 0.0){
					unloading(Dif_Strain);
					Envelope_Check();
										}
				else if (Dif_Strain > 0.0){
					reloading(Dif_Strain);
					Envelope_Check();
										}
				else if (Dif_Strain == 0.0) {
					trialTangent = 0.0;
					trialStress = commitStress;
											}
		}	

return 0;
}

double
Shear01::PosEnvTangent(double trialStrain)
	{
		if (trialStrain <= ecrp)
			return Ep;
		else if (trialStrain <= epp) 
			return E2p;
		else if (trialStrain <= eyp)
			return 0.0;
		else if (trialStrain <= eup)
			return -Edegp;
		else 
			return 0.0;}

double
Shear01::PosEnvStress(double trialStrain)
	{
		if (trialStrain <= ecrp) 
			return Ep*trialStrain;
		else if (trialStrain <= epp)
			return (Vcrp + (E2p*(trialStrain-ecrp)));
		else if (trialStrain <= eyp)
			return Vmaxp;
		else if (trialStrain <= eup)
			return (Vmaxp+(-Edegp*(trialStrain-eyp)));
		else 
			return 0.0;}

double
Shear01::NegEnvTangent(double trialStrain)
	{
		if (trialStrain >= ecrn)
			return En;
		else if (trialStrain >= epn)
			return E2n;
		else if (trialStrain >= eyn)
			return 0.0;
		else if (trialStrain >= eun)
			return -Edegn;
		else 
			return 0.0;}

double
Shear01::NegEnvStress(double trialStrain)
	{
		if (trialStrain >= ecrn)
			return En*trialStrain;
		else if (trialStrain >= epn)
			return (Vcrn + (E2n*(trialStrain-ecrn)));
		else if (trialStrain >= eyn)
			return Vmaxn;
		else if (trialStrain >= eun)
			return (Vmaxn+(-Edegn*(trialStrain-eyn)));
		else 
			return 0.0;}

void
Shear01::reloading(double Dif_Strain)
{
	if (Old_Dif_Strain < 0.0 &&  commitStress > 0.0) {
		TIndicatorNcycle = 1;
		PositiveIncycle(Dif_Strain);
	} 
	else if (TIndicatorNcycle == 1.0 && Tncycle == Cncycle) {
		//TIndicatorNcycle = 1.0;
		PositiveIncycle(Dif_Strain);
	}
	else { //This part is the cyclic rules!
		TIndicatorNcycle = 0;

		//Stage I // Cracking is not happended!!
		if (TminStrain > ecrn) {
			trialTangent = Ep;
			trialStress = Ep*trialStrain;
		}

		//StageII // Cracking happended but peak is not reached!!
		else if (TminStrain > epn) {
			// Other side is in Stage I and II
			if (TmaxStrain < epp) {
				trialTangent = (Vcrp - commitStress) / (ecrp - commitStrain);
				trialStress = commitStress + (trialTangent * Dif_Strain);
			}
			// Otherside is on Stage III
			else if (TmaxStrain < eyp) {
				if (commitStress <= Vcrn) {
					trialTangent = kneg;
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
				/*else if (commitStress <= Vcrp) {
					trialTangent = (Vcrp - commitStress) / (ecrp - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}*/
				else {
					trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
			}
			// Otherside is in Stage VI
			else {
				if (commitStress < 0.0) {
					trialTangent = kneg;
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
				else {
					trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
			}
		}
		// Stage III // Peak is reached but degradation is not onset.
		else if (TminStrain > eyn) {
			// Other side is in Stage I and II
			if (TmaxStrain < epp) {
				if (commitStress <= Vcrn) {
					trialTangent = kneg;
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
				/*else if (commitStress < Vcrp) {
					trialTangent = (Vcrp - commitStress) / (ecrp - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}*/
				else {
					trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
			}
			// Otherside is on Stage III
			else if (TmaxStrain < eyp) {

				if (commitStress <= Vcrn) {
					trialTangent = kneg;
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
				else if (commitStress <= Vcrp) {
					// If the peak is not reach in positive direction.
					/*if (TmaxStrain < epp) {
						trialTangent = (Vcrp - commitStress) / (ecrp - commitStrain);
						trialStress = commitStress + (trialTangent * Dif_Strain);
					}
					else {*/
					// peak is reached in both direction
					double PinchXP = TmaxStrain - (Vmaxp - Vcrp) / (k2neg);
					trialTangent = (Vcrp - commitStress) / (PinchXP - commitStrain);
					trialStress = commitStress + (trialTangent * Dif_Strain);
				}
			//}
			else {
				trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);
				trialStress = commitStress + (trialTangent * Dif_Strain);;
			}
		}
		// Otherside is in Stage VI
		else {
			if (commitStress < 0.0) {
				trialTangent = k2neg;
				trialStress = commitStress + (trialTangent * Dif_Strain);
			}
			else {
				if (fabs(TminStrain) > fabs(TmaxStrain)) {
					trialTangent = (-TenvminStress - commitStress) / (-TminStrain - commitStrain);//
					trialStress = (commitStress + (trialTangent * Dif_Strain));
					TEnvIndicator = 1;
					Envelope_Check();
				}
				else {
					trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);//
					trialStress = (commitStress + (trialTangent * Dif_Strain));
				}
			}
		}
		}
		// Stage VI. // Degradation is stared!!
		else {
			if (commitStress < 0.0) {
				trialTangent = k2neg;
				trialStress = commitStress + (trialTangent * Dif_Strain);
			}
			else {
				if (TmaxStrain < eyp) {
					trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);//
					trialStress = (commitStress + (trialTangent * Dif_Strain));
				}
				else {
					if (fabs(TminStrain) > fabs(TmaxStrain)) {
						trialTangent = (-TenvminStress - commitStress) / (-TminStrain - commitStrain);//
						trialStress = (commitStress + (trialTangent * Dif_Strain));
						TEnvIndicator = 1;
						Envelope_Check();
					}
					else {
						trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);//
						trialStress = (commitStress + (trialTangent * Dif_Strain));
					}
				}
			}
		}
	}
}


void
Shear01::unloading(double Dif_Strain)
{
	if (Old_Dif_Strain > 0.0 && commitStress < 0.0) {
		TIndicatorNcycle = 1;
		NegativeIncycle(Dif_Strain);
	}
	else if (TIndicatorNcycle == 1.0 && Tncycle == Cncycle) {
		//TIndicatorNcycle = 1.0;
		NegativeIncycle(Dif_Strain);
	}
	else { //This part is the cyclic rules!
		TIndicatorNcycle = 0;
		
		// Stage I // Cracking is not happended!!
		if (TmaxStrain < ecrp) { 
			trialTangent = En;
			trialStress = En*trialStrain;
		}
		//StageII // Cracking happended but peak is not reached!!
		else if (TmaxStrain < epp) {
							//other side is in Stage I and II
							if (TminStrain > epn) {
								trialTangent = (Vcrn - commitStress) / (ecrn - commitStrain);
								trialStress = commitStress + (trialTangent * Dif_Strain);
							}
							// Otherside is on Stage III
							else if (TminStrain > eyn) { 
								if (commitStress >= Vcrp) {
									trialTangent = kpos;
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
								/*else if (commitStress >= Vcrn) {
									trialTangent = (Vcrn - commitStress) / (ecrn - commitStrain);
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}*/
								else {
									trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
							}
							// Otherside is in Stage VI
							else {
								if (commitStress > 0.0) {
									trialTangent = kpos;
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
								else {
									trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);//
									trialStress = (commitStress + (trialTangent * Dif_Strain));
								}
							}
						}

		// Stage III
		else if (TmaxStrain < eyp) { 
							// Otherside is in Stage I and II
							if (TminStrain > epn) { 
								if (commitStress >= Vcrp) {
									trialTangent = kpos;
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
								/*else if (commitStress > Vcrn) {
									trialTangent = (Vcrn - commitStress) / (ecrn - commitStrain);
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}*/
								else {
									trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
							}
							// Otherside is in Stage III
							else if (TminStrain > eyn) {
								if (commitStress >= Vcrp) {
									trialTangent = kpos;
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
								else if (commitStress >= Vcrn) {

									/*if (TminStrain > epn) {
										trialTangent = (Vcrn - commitStress) / (ecrn - commitStrain);
										trialStress = commitStress + (trialTangent * Dif_Strain);

									}
									else {*/
										double PinchXN = TminStrain - (Vmaxn - Vcrn) / (k2neg);
										trialTangent = (Vcrn - commitStress) / (PinchXN - commitStrain);
										trialStress = commitStress + (trialTangent * Dif_Strain);
									//}
								}
								else {
									trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
							}
							// otherside is in stage VI
							else {
								if (commitStress > 0.0) {
									trialTangent = k2pos;
									trialStress = commitStress + (trialTangent * Dif_Strain);
								}
								else {
									if (fabs(TmaxStrain) > fabs(TminStrain)) {
										trialTangent = (-TenvmaxStress - commitStress) / (-TmaxStrain - commitStrain);//
										trialStress = (commitStress + (trialTangent * Dif_Strain));
										TEnvIndicator = 1;
										Envelope_Check();
									}
									else {
										trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);//
										trialStress = (commitStress + (trialTangent * Dif_Strain));
								}
							}
						}
					}
		// Stage VI. // Degradation is stared!!
		else { 

			if (commitStress > 0.0) {
				trialTangent = k2pos;
				trialStress = commitStress + (trialTangent * Dif_Strain);
			}
			else {
				if (TminStrain > eyn) {
					trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);//
					trialStress = (commitStress + (trialTangent * Dif_Strain));
				}else{
					if (fabs(TmaxStrain) > fabs(TminStrain)) {
						trialTangent = (-TenvmaxStress - commitStress) / (-TmaxStrain - commitStrain);//
						trialStress = (commitStress + (trialTangent * Dif_Strain));
						TEnvIndicator = 1;
						Envelope_Check();
					}
					else {
						trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);//
						trialStress = (commitStress + (trialTangent * Dif_Strain));
					}
				}
			}
		}
	}
}


double 
Shear01::getStrain(void)
{
  return trialStrain;
}

double 
Shear01::getStress(void)
{
  return trialStress;
}

double 
Shear01::getTangent(void)
{
  return trialTangent;
}

int 
Shear01::commitState(void)
{
  commitStrain = trialStrain;
  commitTangent=trialTangent;
  commitStress = trialStress;
  Old_Dif_Strain = Cur_Dif_Strain;
  CmaxStrain = TmaxStrain;
  CminStrain = TminStrain;
  CmaxStress = TmaxStress;
  CminStress = TminStress;
  CenvmaxStress = TenvmaxStress;
  CenvminStress = TenvminStress;
  Cncycle = Tncycle;
  CEnvIndicator = TEnvIndicator;
  CIndicatorNcycle = TIndicatorNcycle;
    return 0;
}	

int 
Shear01::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;
  Cur_Dif_Strain = Old_Dif_Strain;
  TmaxStrain = CmaxStrain;
  TminStrain = CminStrain;
  TmaxStress = CmaxStress;
  TminStress = CminStress;
  TenvmaxStress = CenvmaxStress;
  TenvminStress = CenvminStress;
  Tncycle = Cncycle;
  TEnvIndicator = CEnvIndicator;
  TIndicatorNcycle = CIndicatorNcycle;
  return 0;
}

int 
Shear01::revertToStart(void)
{
commitStrain = 0.0;
commitTangent = Ep;
commitStress = 0.0;
Old_Dif_Strain = 0.0;
CmaxStrain=0.0;
CminStrain=0.0;
CmaxStress=0.0;
CminStress=0.0;
CenvmaxStress = 0.0;
CenvminStress = 0.0;
Cncycle  = 0.0;
CEnvIndicator = 0;
CIndicatorNcycle = 0;

// this->revertToLastCommit();
    return 0;
}

UniaxialMaterial *
Shear01::getCopy(void)
{
  Shear01 *theCopy =
  new Shear01(this->getTag(),Ep,ecrp,E2p,epp,eyp,eup,En,ecrn,E2n,epn,eyn,eun);
  theCopy->Ep = Ep;
  theCopy->ecrp = ecrp;
  theCopy->E2p = E2p;
  theCopy->epp = epp;
  theCopy->eyp = eyp;
  theCopy->eup = eup;

  theCopy->En = En;
  theCopy->ecrn = ecrn;
  theCopy->E2n = E2n;
  theCopy->epn = epn;
  theCopy->eyn = eyn;
  theCopy->eun = eun;

  theCopy->commitStrain = commitStrain;
  theCopy->commitStress = commitStress;
  theCopy->commitTangent = commitTangent;

  theCopy->Cncycle = Cncycle ;
  theCopy->CEnvIndicator = CEnvIndicator;

  return theCopy;
}

int 
Shear01::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(15);
  data(0) = this->getTag();
  data(1) = Ep;
  data(2) = ecrp;
  data(3) = E2p;
  data(4) = epp;
  data(5) = eyp;
  data(6) = eup;
  data(7) = En;
  data(8) = ecrn;
  data(9) = E2n;
  data(10) = epn;
  data(11) = eyn;
  data(12) = eun;
  data(13) = commitStrain;
  data(14) = commitStress;
  data(15) = commitTangent;
  

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Shear01::sendSelf() - failed to send data\n";

  return res;
}

int 
Shear01::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(15);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Shear01::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    Ep     = data(1);
    ecrp    = data(2);
    E2p    = data(3);
    epp    = data(4);
    eyp    = data(5);  
    eup    = data(6);
    En     = data(7);
    ecrn    = data(8);
    E2n    = data(9);
    epn    = data(10);
    eyn    = data(11);  
    eun    = data(12);
    commitStrain=data(13);
    commitStress=data(14);
    commitTangent=data(15);
    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
	TEnvIndicator = CEnvIndicator;
  }
  return res;
}

void 
Shear01::Print(OPS_Stream &s, int flag)
{
 
}

void
Shear01::setParameters(void)
{
//For positive values
  Vcrp  = ecrp * Ep;    //the cracking force
  Vmaxp = Vcrp + E2p * (epp-ecrp); //the maximum force
  Edegp = Vmaxp/(eup-eyp); //the tangent after onset of shear degredation
  Edegp_2=Vmaxp/(eup-epp); // degredation of type 2
//For negative values
  Vcrn  = ecrn * En;    //the cracking force
  Vmaxn = Vcrn + E2n * (epn-ecrn); //the maximum force
  Edegn = Vmaxn/(eun-eyn); //the tangent after onset of shear degredation
  Edegn_2=Vmaxn/(eun-epn); // degredation of type 2
  
  kpos = (Vmaxp - Vcrn)/(epp - ecrn);
  k2pos = (Vmaxp - Vcrn) / (eyp - ecrn);
  if (k2pos < E2n) { k2pos = E2n; }
  //bu kisim onemli dikkat et !!

  kneg = (Vmaxn - Vcrp)/(epn - ecrp);
  k2neg = (Vmaxn - Vcrp) / (eyn - ecrp);
  if (k2neg < E2p) { k2neg = E2p; }
}

void
Shear01::NumberCycle(void){
if ( Old_Dif_Strain > 0.0 && Cur_Dif_Strain < 0.0 ) {  Tncycle = Cncycle+0.5;}
if ( Old_Dif_Strain < 0.0 && Cur_Dif_Strain > 0.0 ) {  Tncycle = Cncycle+0.5;}
beta = 1 - Tncycle * 0.0100;
}

void
Shear01::PositiveIncycle(double Dif_Strain) {
	trialTangent = (TenvmaxStress - commitStress) / (TmaxStrain - commitStrain);
	trialStress = commitStress + (trialTangent * Dif_Strain);
}

void
Shear01::NegativeIncycle(double Dif_Strain) {
	trialTangent = (TenvminStress - commitStress) / (TminStrain - commitStrain);
	trialStress = commitStress + (trialTangent * Dif_Strain);
}

void
Shear01::Envelope_Check(void) {

	if (trialStrain >= ecrp) {
		if (trialStress > PosEnvStress(trialStrain)) { 
			trialStress = PosEnvStress(trialStrain); 
			trialTangent = (trialStress - commitStress) / (trialStrain - commitStrain);
			TEnvIndicator = 0;
		}
	}
	if (trialStrain <= ecrn) {
		if (trialStress < NegEnvStress(trialStrain)) { 
			trialStress = NegEnvStress(trialStrain); 
			trialTangent = (trialStress - commitStress) / (trialStrain - commitStrain);
			TEnvIndicator = 0;
		}
	}


}
