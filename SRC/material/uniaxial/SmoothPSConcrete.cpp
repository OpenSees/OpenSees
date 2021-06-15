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
                                                                        
// Function contributed by Quan Gu & Michele Barbato

// $Revision: 1.3 $
// $Date: 2009-07-23 23:43:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SmoothPSConcrete.cpp,v $
                                                                        
// Description: This file contains the class definition for 
// SmoothPSConcrete.h 
//   - Popovics-Saenz envelope
//   - No tension
//   - Smoooth transition to  a linear unloading/reloading
//
// What: "@(#) SmoothPSConcrete.C"


#include <SmoothPSConcrete.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <string.h>
#include <elementAPI.h>
#include <float.h>
#include <MaterialResponse.h>
# define MAT_TAG_SmoothPSConcrete 35457

void* OPS_SmoothPSConcrete()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 6 || argc > 9) {
	opserr << "WARNING invalid number of arguments\n";
	opserr << "Want: uniaxialMaterial SmoothPSConcrete tag? fc? fu? Ec? <eps0?> <epsu?> <eta?>\n";
	return 0;
    }    
      
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete tag\n";
	return 0;
    }
    
    // double fu, Ec, fc;
    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete double inputs\n";
	return 0;
    }
    
    // double eps0=0.002;
    // double epsu=0.005;
    // double eta=0.2;
    double opt[3] = {0.002, 0.005, 0.2};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 3) numdata = 3;
    if (OPS_GetDoubleInput(&numdata, opt) < 0) {
	opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete double inputs\n";
	return 0;
    }

    // Parsing was successful, allocate the material
    return new SmoothPSConcrete( tag, data[0], data[1], data[2],
				 opt[0], opt[1], opt[2]);       
}

SmoothPSConcrete::SmoothPSConcrete
(int tag, double FPC,double FPCU, double EC, double EPSC0,  double EPSCU,  double ETA)
  :UniaxialMaterial(tag, MAT_TAG_SmoothPSConcrete),
   fc(FPC), eps0(EPSC0), fcu(FPCU), epsu(EPSCU),Ec(EC),eta(ETA)  
{
     
///// refer to the same function revertToStart();
	this->revertToStart ();
	this->compute_epsmax(&epsmax,&sigmax);
	
	SHVs =0;
	parameterID =0;
}


SmoothPSConcrete::~SmoothPSConcrete ()
{
   if (SHVs != 0)  delete SHVs;
}

int SmoothPSConcrete::setTrialStrain (double strain, double strainRate)
{


  // Retrieve history variables from Past
            Tsig  = Csig;
            TEt    = CEt;
            Tepsr  = Cepsr;
            Tsigr  = Csigr;
            Tepsp  = Cepsp;
            TEur   = CEur;
            Tflag  = Cflag;
            Tepsr1 = Cepsr1;
            Tepsr2 = Cepsr2;
            Tsigr1 = Csigr1;
            Tsigr2 = Csigr2;
            TEt2   = CEt2;
            Tepsc  = strain;
            Tdepsc = Tepsc-Cepsc ; // total strain increment

   		if (fabs(Tdepsc) < DBL_EPSILON) {return 0; };
           
                
            
                
                if (Tflag == -1) {// virgin material
                    if (Tepsc >= 0.0001) {
                        Tflag = -2;
                        Tepsr1 = 0.0001;
                        Tsigr1 = 0;
                        Tepsr2 = -0.0001;
                        Monotonic_Envelope (Tepsr2, &Tsigr2,&TEt2);
                        Tsig = 0;
                        TEt  = 0;
					}
                    else if ((Tepsc >= 0) && (Tepsc < 0.0001)) {
                        Tflag = -3; 
                        Tsig = 0;
                        TEt  = 0;
					}
                    else if ((Tepsc < 0) && (Tepsc > -0.0001)) {
                        Tflag = -3; 
                        Monotonic_Envelope (Tepsc, &Tsig,&TEt);
					}
					
                    else if (Tepsc <= -0.0001){
                        Tflag = 0;
                        Monotonic_Envelope (Tepsc, &Tsig,&TEt);
					}
				} //if (Tflag == -1)
                    
                else if (Tflag == -2){ //material never been in compression
                    if (Tepsc >= 0.0001){
                        //flag = -2;
                        Tsig = 0;
                        TEt  = 0;
					}
                    else if ((Tepsc > -0.0001) && (Tepsc < 0.0001))  //flag = -2; //must define smooth transition curve
                          Transition_r(Tepsc, Tepsr1,Tepsr2,Tsigr1,Tsigr2,0.001,-0.001, 0, TEt2, &Tsig,&TEt);

                    else if (Tepsc <= -0.0001) {
                        Tflag = 0;
                        Monotonic_Envelope (Tepsc, &Tsig,&TEt);
					}
                    
				}  // else if (Tflag==-2 )
				
                else if (Tflag == -3) {//temporary curve at beginning

                    if ((Tepsc >= 0) && (Tepsc < 0.0001)) { // flag = -3;
                       
                        if (Tdepsc > 0){ //still on temporary curve
                            Tsig = 0;
                            TEt  = 0; 
						}
                        else {       // depsc < 0 //smooth transition curve
                            Tflag = -2;
                            Tepsr1 = Tepsc-Tdepsc;
                            Tsigr1 = 0;
                            Tepsr2 = -Tepsr1;
                            Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
                            Transition_r (Tepsc, Tepsr1,Tepsr2,Tsigr1,Tsigr2,0.001,-0.001,0,TEt2,&Tsig,&TEt);
                        }
					} //  if ((Tepsc >= 0) && (Tepsc < 0.0001)) 

				    else if ((Tepsc < 0) && (Tepsc > -0.0001)) {     // flag = -3;
                   
						if (Tdepsc < 0)  //still on temporary curve
							Monotonic_Envelope (Tepsc ,&Tsig,&TEt);
						else {// depsc > 0 //smooth transition curve
							Tflag = -2;
							Tepsr2 = Tepsc-Tdepsc;
							Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
							Tepsr1 = -Tepsr2;
							Tsigr1 = 0;
							Transition_r (Tepsc,Tepsr1,Tepsr2,Tsigr1,Tsigr2,0.001,-0.001, 0, TEt2,&Tsig,&TEt);
						};
					} //else if ((Tepsc < 0) && (Tepsc > -0.0001))

                    else if (Tepsc >=0.0001) {
                        Tflag = -2;
                        Tepsr1 = 0.0001;
                        Tsigr1 = 0;
                        Tepsr2 = -0.0001;
                        Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
                        Tsig = 0;
                        TEt  = 0;
					}  //else if (Tepsc >=0.0001)
                    else if (Tepsc <= -0.0001){
                        Tflag = 0;
                        Monotonic_Envelope (Tepsc,&Tsig,&TEt);
                    } //else if (epsc <= -0.0001)

				}  //   else if (flag == -3)
				
                else if (Tflag == 0) {// monotonic envelope
                    if (Tdepsc < 0) // keep loading on monotonic envelope
                        //flag = 0;
                        Monotonic_Envelope (Tepsc,&Tsig,&TEt);
                    else {// depsc > 0 // unloading
                        Tepsr = (Tepsc-Tdepsc);
                        Tsigr = Csig;
                        Compute_epsp ();
                        TEur = fabs(Tsigr/(Tepsr-Tepsp));

                        if (Tepsc < Tepsr + eta*(Tepsp-Tepsr)){ //temporary curve
                            Tflag = 3;
                            Tsig = TEur*(Tepsc-Tepsp);
                            TEt  = TEur;
						}
                        else {
                            Tepsr1 = Tepsr +eta*(Tepsp-Tepsr);
                            Tepsr2 = Tepsr -eta*(Tepsp-Tepsr);
                            Tsigr1 = TEur*(Tepsr1-Tepsp);
                            Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
                            if ((Tepsc >= Tepsr + eta*(Tepsp-Tepsr)) && (Tepsc <= Tepsp - eta*(Tepsp-Tepsr))) {
                                Tflag = 1;
                                Tsig = TEur*(Tepsc-Tepsp);
                                TEt  = TEur;
							}
                            else if ((Tepsc > Tepsp - eta*(Tepsp-Tepsr)) && (Tepsc < Tepsp + eta*(Tepsp-Tepsr))){
                                Tflag = 1;  
                                Transition_p (eta*(Tepsp-Tepsr));
							}
                              
                            else if (Tepsc >= Tepsp + eta*(Tepsp-Tepsr)){
                                Tflag = 2;
                                Tsig = 0;
                                TEt  = 0;
                            }
                        } //else

                    } //else
				} //else if (flag == 0)
				

                else if (Tflag == 1) {// unloading curve with fixed end-points
                    if (Tepsc >= Tepsp + eta*(Tepsp-Tepsr)){
                        Tflag = 2;
                        Tsig = 0;
                        TEt  = 0;
					}
                    else if ((Tepsc < Tepsp + eta*(Tepsp-Tepsr)) && (Tepsc > Tepsp - eta*(Tepsp-Tepsr)))//flag = 1;
                        Transition_p (eta*(Tepsp-Tepsr));

                    else if ((Tepsc <= Tepsp - eta*(Tepsp-Tepsr)) && (Tepsc >= Tepsr + eta*(Tepsp-Tepsr))){     //flag = 1;
                        Tsig = TEur*(Tepsc-Tepsp);
                        TEt  = TEur;
					}
                    else if ((Tepsc < Tepsr + eta*(Tepsp-Tepsr)) && (Tepsc > Tepsr - eta*(Tepsp-Tepsr)))    //flag = 1; 
                        Transition_r(Tepsc,Tepsr1,Tepsr2,Tsigr1,Tsigr2,Tepsr+eta*(Tepsp-Tepsr),Tepsr-eta*(Tepsp-Tepsr),TEur,TEt2,&Tsig,&TEt);
                    else if (Tepsc <= Tepsr - eta*(Tepsp-Tepsr)){
                        Tflag = 0;
                        Monotonic_Envelope (Tepsc,&Tsig,&TEt);
					}
                    
				
				}
                //else if (flag == 1) 
				

                else if (Tflag == 2) {// friction state
                    if (Tepsc >= Tepsp + eta*(Tepsp-Tepsr)){ //flag = 2;
                       
                        Tsig = 0;
                        TEt  = 0;
					}

                    else if ((Tepsc < Tepsp + eta*(Tepsp-Tepsr)) && (Tepsc > Tepsp - eta*(Tepsp-Tepsr)))  {
                        Tflag = 1;  
                        Transition_p (eta*(Tepsp-Tepsr));
					}

                    else if ((Tepsc <= Tepsp - eta*(Tepsp-Tepsr)) && (Tepsc >= Tepsr + eta*(Tepsp-Tepsr))){
                        Tflag = 1;
                        Tsig = TEur*(Tepsc-Tepsp);
                        TEt  = TEur;
					}

                    else if ((Tepsc < Tepsr + eta*(Tepsp-Tepsr)) && (Tepsc > Tepsr - eta*(Tepsp-Tepsr))){
                        Tflag = 1;  
                        Transition_r(Tepsc,Tepsr1,Tepsr2,Tsigr1,Tsigr2,Tepsr+eta*(Tepsp-Tepsr),Tepsr-eta*(Tepsp-Tepsr),TEur,TEt2,&Tsig,&TEt);
					}
                    else if (Tepsc <= Tepsr - eta*(Tepsp-Tepsr)){
                        Tflag = 0;
                        Monotonic_Envelope (Tepsc,&Tsig,&TEt);
					
					
					}
                     
				
				}
                //    elseif (flag == 2) // friction state
			 


                else if (Tflag == 3) {// unloading curve with variable end-points
                    if ((Tepsc > Tepsr) && (Tepsc < Tepsr + eta*(Tepsp-Tepsr))){ // flag = -3;
                       
                        if (Tdepsc > 0){ //still on temporary curve
                            Tsig = TEur*(Tepsc-Tepsp);
                            TEt  = TEur; 
						}
                        else { // depsc < 0 //smooth transition curve
                            Tflag = 1;
                            Tepsr1 = Tepsc-Tdepsc;
                            Tsigr1 = TEur*(Tepsr1-Tepsp);
                            Tepsr2 = 2*Tepsr-Tepsr1;
                            Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
                            Transition_r (Tepsc,Tepsr1,Tepsr2,Tsigr1,Tsigr2,Tepsr+eta*(Tepsp-Tepsr),Tepsr-eta*(Tepsp-Tepsr),TEur,TEt2,&Tsig,&TEt);
                        }
					}
					else if (Tepsc <= Tepsr - eta*(Tepsp-Tepsr)){ 
						Tflag = 0;
                        Monotonic_Envelope (Tepsc,&Tsig,&TEt); 
					}
                    else if ((Tepsc <= Tepsr) && (Tepsc > Tepsr - eta*(Tepsp-Tepsr))){
                        Tflag = 1;
                        Tepsr1 = Tepsc-Tdepsc;
                        Tsigr1 = TEur*(Tepsr1-Tepsp);
                        Tepsr2 = 2*Tepsr-Tepsr1;
                        Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);
                        Transition_r (Tepsc,Tepsr1,Tepsr2,Tsigr1,Tsigr2,Tepsr+eta*(Tepsp-Tepsr),Tepsr-eta*(Tepsp-Tepsr),TEur,TEt2,&Tsig,&TEt);
                    }
					else {// epsc >= epsr + eta*(epsp-epsr))
                        Tepsr1 = Tepsr + eta*(Tepsp-Tepsr);
                        Tsigr1 = TEur*(Tepsr1-Tepsp);
                        Tepsr2 = Tepsr - eta*(Tepsp-Tepsr);
                        Monotonic_Envelope (Tepsr2,&Tsigr2,&TEt2);

                        if (((Tepsc >= Tepsr + eta*(Tepsp-Tepsr))) && (Tepsc <= Tepsp - eta*(Tepsp-Tepsr))){
                            Tflag = 1;
                            Tsig = TEur*(Tepsc-Tepsp);
                            TEt  = TEur;
						}
                        else if( (Tepsc > Tepsp - eta*(Tepsp-Tepsr)) && (Tepsc < Tepsp + eta*(Tepsp-Tepsr))){
                            Tflag = 1;
                            Transition_p (eta*(Tepsp-Tepsr));
						}
                        else if (Tepsc >= Tepsp + eta*(Tepsp-Tepsr)){
                            Tflag = 2;
                            Tsig = 0;
                            TEt = 0;
						}
					
					}//else
				
				
				} // else if (Tflag == 3)
            

  
  return 0;
}



int 
SmoothPSConcrete::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
  // Set trial strain
  setTrialStrain (strain, strainRate);
  stress =Tsig;
  tangent=TEt;
  
  return 0;
}



double SmoothPSConcrete::getStress ()
{
   return Tsig;
}

double SmoothPSConcrete::getStrain ()
{
   return Tepsc;
}

double SmoothPSConcrete::getTangent ()
{
   return TEt;
}

int SmoothPSConcrete::commitState ()
{
   // History variables
   // save history variables   
   Csig   = Tsig;
   CEt    = TEt;
   Cepsr  = Tepsr;
   Csigr  = Tsigr;
   Cepsp  = Tepsp;
   CEur   = TEur;
   Cflag  = Tflag;
   Cepsr1 = Tepsr1;
   Cepsr2 = Tepsr2;
   Csigr1 = Tsigr1;
   Csigr2 = Tsigr2;
   CEt2   = TEt2;
   Cepsc  = Tepsc;  
   


   return 0;
}

int SmoothPSConcrete::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   Tepsc  = Cepsc; 
   Tsig   = Csig;
   TEt    = CEt;
   Tepsr  = Cepsr;
   Tsigr  = Csigr;
   Tepsp  = Cepsp;
   TEur   = CEur;
   Tflag  = Cflag;
   Tepsr1 = Cepsr1;
   Tepsr2 = Cepsr2;
   Tsigr1 = Csigr1;
   Tsigr2 = Csigr2;
   TEt2   = CEt2;


   return 0;
}

int SmoothPSConcrete::revertToStart ()
{

	    	Cepsc = 0.0;
	        Csig  = 0;
            CEt   = Ec;      
			Cepsr = 0;
            Csigr = 0;
            Cepsp = 0;
            CEur  = CEt;
            Cflag  = -1;
            Cepsr1 = 0;
            Cepsr2 = 0;
            Csigr1 = 0;
            Csigr2 = 0;
            CEt2   = CEt;

			Tepsc  = Cepsc; 
			Tsig   = Csig;
			TEt    = CEt;
            Tepsr = 0;
            Tsigr = 0;
            Tepsp = 0;
            TEur  = CEt;
            Tflag  = -1;
            Tepsr1 = 0;
            Tepsr2 = 0;
            Tsigr1 = 0;
            Tsigr2 = 0;
            TEt2   = CEt;



// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

   return 0;
}

UniaxialMaterial* SmoothPSConcrete::getCopy (void)
{
// Note: it will not work if called after running analysis since the stress state is not copied herein.

	SmoothPSConcrete* theCopy = new SmoothPSConcrete(this->getTag(),
                                    fc,fcu,  Ec, eps0, epsu,eta);

   theCopy->revertToStart ();

   return theCopy;
}





int SmoothPSConcrete::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(20);
   data(0) = this->getTag();

   // Material properties
   data(1) = fc;
   data(2) = eps0;
   data(3) = fcu;
   data(4) = epsu;
   data(5) = eta;

   // History variables from last converged state
   data(6) = Cepsc; 
   data(7) = Csig;
   data(8) = CEt;
   data(9) = Cepsr;
    
   data(10) = Csigr;
   data(11) = Cepsp;
   data(12) = CEur;
   data(13) = Cflag;
   data(14) = Cepsr1;
   data(15) = Cepsr2;
   data(16) = Csigr1;
   data(17) = Csigr2;
   data(18) = CEt2;

   data(19) = Ec;
  
   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "SmoothPSConcrete::sendSelf() - failed to send data\n";

   return res;
}

int SmoothPSConcrete::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(20);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "SmoothPSConcrete::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {

   this->setTag(int(data(0)));
   fc=data(1);
   eps0=data(2) ;
   fcu= data(3);
   epsu = data(4);
   eta =  data(5);

   // History variables from last converged state
   Cepsc = data(6)  ; 
   Csig = data(7);
   CEt = data(8);
   Cepsr = data(9) ;
    
   Csigr = data(10);
   Cepsp = data(11);
   CEur  = data(12);
   Cflag = data(13);
   Cepsr1 = data(14);
   Cepsr2 = data(15) ;
   Csigr1 = data(16);
   Csigr2 = data(17);
   CEt2 = data(18);

   Ec =  data(19);
      // Set trial state variables
   revertToLastCommit();
   }

   return res;
}

void SmoothPSConcrete::Print (OPS_Stream& s, int flag)
{
   s << "SmoothPSConcrete, tag: " << this->getTag() << endln;
   s << "  fc: " << fc << endln;
   s << "  eps0: " << eps0 << endln;
   s << "  fcu: " << fcu << endln;
   s << "  epsu: " << epsu << endln;
   s << "  eta: " << eta << endln;
   s << "  Ec: " << Ec << endln;
}



/*Response* 
SmoothPSConcrete::setResponse(const char **argv, int argc, Information &matInfo)
{


  if (strcmp(argv[0],"stressSensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+100, this->getStress());  // use only size of matrix

  }
  else if (strcmp(argv[0],"strainSensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+500, this->getStrain());
  }

  
 
  //by default, See if the response is one of the defaults
  Response *res = UniaxialMaterial::setResponse(argv, argc, matInfo);
  if (res != 0)      return res;
  else	return 0;
}

int 
SmoothPSConcrete::getResponse(int responseID, Information &matInfo)
{
  
	if (responseID>100 && responseID<500) {
		return matInfo.setDouble(this->getStressSensitivity(responseID-100, false));
	} // if

	else if (responseID>500) {
		return matInfo.setDouble(this->getStrainSensitivity(responseID-500));

	} // if

	else

	  // Just call the base class method ... don't need to define
	  // this function, but keeping it here just for clarity
	  return UniaxialMaterial::getResponse(responseID, matInfo);	
	
}
*/

// AddingSensitivity:BEGIN ///////////////////////////////////
int
SmoothPSConcrete::setParameter(const char **argv, int argc, Parameter &param)
{
	                
                    
	if (strcmp(argv[0],"fc") == 0) {// Compressive strength
		return param.addObject(1, this);
	}
	if ((strcmp(argv[0],"epsco") == 0)||(strcmp(argv[0],"epso") == 0)) {// Strain at compressive strength
		return param.addObject(2, this);
	}
	if ((strcmp(argv[0],"epsu") == 0)||(strcmp(argv[0],"epscu") == 0)) {// Strain at crushing strength
		return param.addObject(3, this);
	}
	if (strcmp(argv[0],"fcu") == 0) {// Crushing strength
		return param.addObject(4, this);
	}
	if (strcmp(argv[0],"Ec") == 0) {// initial stiffness
		return param.addObject(5, this);
	}
	if (strcmp(argv[0],"eta") == 0) {// smoothing parameter
		return param.addObject(6, this);
	}

	else {
		opserr << "WARNING: Could not set parameter in SmoothPSConcrete! " << endln;
		return -1;
	}
}
    
                            

int
SmoothPSConcrete::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->fc = info.theDouble;
		break;
	case 2:
		this->eps0 = info.theDouble;
		break;
	case 3:
		this->epsu = info.theDouble;
		break;
	case 4:
		this->fcu = info.theDouble;
		break;
	case 5:
		this->Ec = info.theDouble;
		break;
	case 6:
		this->eta = info.theDouble;
		break;
	default:
		break;
	}
    this->revertToStart();    
	this->compute_epsmax(&epsmax,&sigmax);
	return 0;
}




int
SmoothPSConcrete::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double
SmoothPSConcrete::getStressSensitivity(int gradNumber, bool conditional){


	//	double TstressSensitivity =0.0;
	
		if (conditional==false){
			if (SHVs==0) {
				opserr<<"warning: recordering SmoothPSConcrete::getStressSensitivity! SHVs=0";
				return 0.0;
			}	
		
		return (*SHVs)(1,gradNumber);
		}


     
            
      //  % extract history variables' values at previous step
            double sigp  = Csig;
            double epsr  = Cepsr;
            double sigr  = Csigr;
            double epsp  = Cepsp;
            double Eur   = CEur;
            double flag  = Cflag;
            double epsr1 = Cepsr1;
            double epsr2 = Cepsr2;
            double sigr1 = Csigr1;
            double sigr2 = Csigr2;
            double Et2   = CEt2;
           
            double epsc  = Tepsc; // % total strain
            double depsc = Tdepsc;   //% total strain increment
            
       
            double depsdh  = 0; // for unconditional 

            double depsdhp = 0;
			double dsigdhp  = 0;
            double depsrdh  = 0;
            double dsigrdh = 0;
            double depspdh  = 0;
            double dEurdh   = 0;
            double depsr1dh = 0;
            double depsr2dh = 0;
            double dsigr1dh = 0;
            double dsigr2dh = 0;
            double dEt2dh   = 0;
			
			if (SHVs !=0) {  // Retrieve history variables
				depsdhp   = (*SHVs)(0,(gradNumber));
				dsigdhp   = (*SHVs)(1,(gradNumber));
				depsrdh   = (*SHVs)(2,(gradNumber));
				dsigrdh   = (*SHVs)(3,(gradNumber));
				depspdh   = (*SHVs)(4,(gradNumber));
				dEurdh    = (*SHVs)(5,(gradNumber));
				depsr1dh  = (*SHVs)(6,(gradNumber));
				depsr2dh  = (*SHVs)(7,(gradNumber));
				dsigr1dh  = (*SHVs)(8,(gradNumber));
				dsigr2dh  = (*SHVs)(9,(gradNumber));
				dEt2dh    = (*SHVs)(10,(gradNumber));		
            }


           // % setup the parameter derivatives
            double dfcdh   = 0;
            double deps0dh = 0;
            double depsudh = 0;
            double dfudh   = 0;
            double dEcdh   = 0;
            double detadh  = 0;

			if (parameterID == 1) {
				dfcdh = 1.0;
			}
			else if (parameterID == 2) {
				deps0dh = 1.0;
			}
			else if (parameterID == 3) {
				depsudh = 1.0;
			}
			else if (parameterID == 4) {
				dfudh = 1.0;
			}
			else if (parameterID == 5) {
				dEcdh = 1.0;
			}
			else if (parameterID == 6) {
				detadh = 1.0;
			}

           double dsigdh = 0.0; 
			
           if (fabs(Tdepsc) < DBL_EPSILON) {//  % ------------------------------ total strain is not changing
                dsigdh = dsigdhp - CEt*depsdhp;
			}
                
            else { //%(depsc ~= 0)
                
                if (flag == -1) { // % virgin material
                    if (epsc >= 0.0001) {
                        depsr1dh = 0.0;
                        dsigr1dh = 0.0;
                        depsr2dh = 0.0;
                        dsigr2dh = Monotonic_Envelope_sens (-0.0001,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
                        dEt2dh = Monotonic_Envelope_Et_sens (-0.0001,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh); 
                        dsigdh = 0;
					}
                    else if ((epsc >= 0) & (epsc < 0.0001)) {
                        dsigdh = 0;
					}
                    else if ((epsc < 0) & (epsc > -0.0001)) {
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else if (epsc <= -0.0001) {
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				} //  if (flag == -1)
                else if (flag == -2) { //%material never been in compression
                    if (epsc >= 0.0001){
                        dsigdh = 0;
					}
                    else if ((epsc > -0.0001) & (epsc < 0.0001)){ //%must define smooth transition curve
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= -0.0001){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
                    }
				} // elseif (flag == -2)
                else if (flag == -3) { //%temporary curve at beginning
                    if ((epsc >= 0) & (epsc < 0.0001)) {
                        if (depsc > 0 ) { //%still on temporary curve
                            dsigdh = 0; 
						}
                        else { //% depsc < 0 %smooth transition curve
                            epsr1 = epsc-depsc;
                            sigr1 = 0;
                            epsr2 = -epsr1;
                            depsr1dh = depsdhp;
                            dsigr1dh = 0;
                            depsr2dh = -depsdhp;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }


					} //if ((epsc >= 0) & (epsc < 0.0001))
                    else if ((epsc < 0) & (epsc > -0.0001)) {
                        //% flag = -3;
                        if (depsc < 0) { // %still on temporary curve
                            dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
						}
                        else{ // % depsc > 0 %smooth transition curve
                            epsr2 = epsc-depsc;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            epsr1 = -epsr2;
                            sigr1 = 0;
                            depsr2dh = depsdhp;
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh   = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            depsr1dh = -depsdhp;
                            dsigr1dh = 0;
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }
					}//else if ((epsc < 0) & (epsc > -0.0001))
                    else if (epsc >=0.0001) {
                        epsr1 = 0.0001;
                        sigr1 = 0;
                        epsr2 = -0.0001;
                        depsr1dh = 0;
                        depsr2dh = 0;
                        dsigr1dh = 0;
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        dsigdh = 0;
					}//else if (epsc >=0.0001)
                    else if (epsc <= -0.0001){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
    			
				}//else if (flag == -3) 
                else if (flag == 0) { // % monotonic envelope
                    if (depsc < 0) { // % keep loading on monotonic envelope
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else { //% depsc > 0 % unloading
                        epsr = (epsc-depsc);
                        sigr = sigp;

                        //epsp = Compute_epsp (epsr,sigr,eps0,fc,Ec);
						    // Modified Popovics-Saenz concrete stress-strain relationship: compute epsp

						if (fabs(epsr)< eps0)	epsp = -(fabs(epsr) - fabs(sigr)/Ec);
						else  epsp = -(eps0 - fc/Ec);
                      
						Eur = fabs(sigr/(epsr-epsp));
                        depsrdh = depsdhp;
                        dsigrdh = dsigdhp;
                        depspdh = Compute_depspdh(epsr,sigr,depsrdh,dsigrdh,deps0dh,dfcdh,dEcdh);
						double sign = 1.0;
						if ((sigr/(epsr-epsp))<0) sign = -1.0;
						else if ((sigr/(epsr-epsp))==0) sign = 0.0;

                        dEurdh = sign*(dsigrdh*(epsr-epsp)-sigr*(depsrdh-depspdh))/pow((epsr-epsp),2); 
                        
						if (epsc < epsr + eta*(epsp-epsr)) { // %temporary curve
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
						}
                        else {
                            epsr1 = epsr + eta*(epsp-epsr);
                            epsr2 = epsr - eta*(epsp-epsr);
                            sigr1 = Eur*(epsr1-epsp);
                            Monotonic_Envelope(epsr2,&sigr2,&Et2);
                            depsr1dh = depsrdh + eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                            depsr2dh = depsrdh - eta*(depspdh-depsrdh) - detadh*(epsp-epsr);
                            dsigr1dh = Eur*(depsr1dh-depspdh)+dEurdh*(epsr1-epsp);
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);    
                            if ((epsc >= epsr + eta*(epsp-epsr)) & (epsc <= epsp - eta*(epsp-epsr))){
                                dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
							}
                            else if ((epsc > epsp - eta*(epsp-epsr)) & (epsc < epsp + eta*(epsp-epsr))){
                                double dDeltadh = eta*(depspdh-depsrdh)+ detadh*(epsp-epsr);
                                dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
							}
                            else if (epsc >= epsp + eta*(epsp-epsr)){
                                dsigdh = 0;
							}
                            
                        }
                    } //else //% depsc > 0 % unloading
				} //else if (flag == 0) {   
                else if (flag == 1) {//% unloading curve with fixed end-points
                    if (epsc >= epsp + eta*(epsp-epsr)){
                        dsigdh = 0;
					}
                    else if ((epsc < epsp + eta*(epsp-epsr)) & (epsc > epsp - eta*(epsp-epsr))){
                        double dDeltadh = eta*(depspdh-depsrdh) + detadh*(epsp-epsr);                        
                        dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
					}
                    else if ((epsc <= epsp - eta*(epsp-epsr)) & (epsc >= epsr + eta*(epsp-epsr))){
                        dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
					}
                    else if ((epsc < epsr + eta*(epsp-epsr)) & (epsc > epsr - eta*(epsp-epsr))) {
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= epsr - eta*(epsp-epsr)){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				}//else if (flag == 1) 
                else if (flag == 2){ // % friction state
                    if (epsc >= epsp + eta*(epsp-epsr)){
                        dsigdh = 0;
					}
                    else if ((epsc < epsp + eta*(epsp-epsr)) & (epsc > epsp - eta*(epsp-epsr))){
                        double dDeltadh=eta*(depspdh-depsrdh) + detadh*(epsp-epsr); 
                        dsigdh = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
					}
                    else if ((epsc <= epsp - eta*(epsp-epsr)) & (epsc >= epsr + eta*(epsp-epsr))){
                        dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
					}
                    else if ((epsc < epsr + eta*(epsp-epsr)) & (epsc > epsr - eta*(epsp-epsr))){
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= epsr - eta*(epsp-epsr)){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				} //else if (flag == 2){    
                else if (flag == 3) { //% unloading curve with variable end-points
                    if ((epsc > epsr) & (epsc < epsr + eta*(epsp-epsr))){
                        if (depsc > 0){// %still on temporary curve
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp); 
						}
                        else {// % depsc < 0 %smooth transition curve
                            epsr1 = epsc-depsc;
                            sigr1 = Eur*(epsr1-epsp);
                            epsr2 = 2*epsr-epsr1;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            depsr1dh = depsdhp;
                            dsigr1dh = dsigdhp;
                            depsr2dh = 2*depsrdh-depsr1dh;
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }
					} //if ((epsc > epsr) & (epsc < epsr + eta*(epsp-epsr)))
                    else if (epsc <= epsr - eta*(epsp-epsr)){ 
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else if ((epsc <= epsr) & (epsc > epsr - eta*(epsp-epsr))){
                        epsr1 = epsc-depsc;
                        sigr1 = Eur*(epsr1-epsp);
                        epsr2 = 2*epsr-epsr1;
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        depsr1dh = depsdhp;
                        dsigr1dh = dsigdhp;
                        depsr2dh = 2*depsrdh-depsr1dh;
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else { // % epsc >= epsr + eta*(epsp-epsr))
                        epsr1 = epsr + eta*(epsp-epsr);
                        sigr1 = Eur*(epsr1-epsp);
                        epsr2 = epsr - eta*(epsp-epsr);
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        depsr1dh = depsrdh + eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                        dsigr1dh = dEurdh*(epsr1-epsp)+Eur*(depsr1dh-depspdh);
                        depsr2dh = depsrdh - eta*(depspdh-depsrdh) - detadh*(epsp-epsr);
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        if (((epsc >= epsr + eta*(epsp-epsr))) & (epsc <= epsp - eta*(epsp-epsr))){
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
						}
                        else if ((epsc > epsp - eta*(epsp-epsr)) & (epsc < epsp + eta*(epsp-epsr))){
                            double dDeltadh = eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                            dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
						}
                        else if ((epsc >= epsp + eta*(epsp-epsr))){
                            dsigdh = 0;
                        }
                    }
				} //else if (flag == 3)  
        } //else { //%(depsc ~= 0)
		
			return dsigdh;
			
	
};
		


int
SmoothPSConcrete::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{

	
	
   
            
      //  % extract history variables' values at previous step
            double sigp  = Csig;
            double epsr  = Cepsr;
            double sigr  = Csigr;
            double epsp  = Cepsp;
            double Eur   = CEur;
            double flag  = Cflag;
            double epsr1 = Cepsr1;
            double epsr2 = Cepsr2;
            double sigr1 = Csigr1;
            double sigr2 = Csigr2;
            double Et2   = CEt2;
           
            double epsc  = Tepsc; // % total strain
            double depsc = Tdepsc;   //% total strain increment
            
       
            double depsdh  = TstrainSensitivity; // for unconditional 

            double depsdhp = 0;
			double dsigdhp  = 0;
            double depsrdh  = 0;
            double dsigrdh = 0;
            double depspdh  = 0;
            double dEurdh   = 0;
            double depsr1dh = 0;
            double depsr2dh = 0;
            double dsigr1dh = 0;
            double dsigr2dh = 0;
            double dEt2dh   = 0;
			
			if (SHVs ==0) {  // Retrieve history variables
				SHVs= new Matrix(11,numGrads);
            }

			else if (SHVs !=0) {  // Retrieve history variables
				depsdhp   = (*SHVs)(0,(gradNumber));
				dsigdhp   = (*SHVs)(1,(gradNumber));
				depsrdh   = (*SHVs)(2,(gradNumber));
				dsigrdh   = (*SHVs)(3,(gradNumber));
				depspdh   = (*SHVs)(4,(gradNumber));
				dEurdh    = (*SHVs)(5,(gradNumber));
				depsr1dh  = (*SHVs)(6,(gradNumber));
				depsr2dh  = (*SHVs)(7,(gradNumber));
				dsigr1dh  = (*SHVs)(8,(gradNumber));
				dsigr2dh  = (*SHVs)(9,(gradNumber));
				dEt2dh    = (*SHVs)(10,(gradNumber));		
            }


           // % setup the parameter derivatives
            double dfcdh   = 0;
            double deps0dh = 0;
            double depsudh = 0;
            double dfudh   = 0;
            double dEcdh   = 0;
            double detadh  = 0;

			if (parameterID == 1) {
				dfcdh = 1.0;
			}
			else if (parameterID == 2) {
				deps0dh = 1.0;
			}
			else if (parameterID == 3) {
				depsudh = 1.0;
			}
			else if (parameterID == 4) {
				dfudh = 1.0;
			}
			else if (parameterID == 5) {
				dEcdh = 1.0;
			}
			else if (parameterID == 6) {
				detadh = 1.0;
			}

           double dsigdh = 0.0; 
			
           if (fabs(Tdepsc) < DBL_EPSILON) {//  % ------------------------------ total strain is not changing
                dsigdh = dsigdhp - CEt*depsdhp;
			}
                
            else { //%(depsc ~= 0)
                
                if (flag == -1) { // % virgin material
                    if (epsc >= 0.0001) {
                        depsr1dh = 0.0;
                        dsigr1dh = 0.0;
                        depsr2dh = 0.0;
                        dsigr2dh = Monotonic_Envelope_sens (-0.0001,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
                        dEt2dh = Monotonic_Envelope_Et_sens (-0.0001,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh); 
                        dsigdh = 0;
					}
                    else if ((epsc >= 0) & (epsc < 0.0001)) {
                        dsigdh = 0;
					}
                    else if ((epsc < 0) & (epsc > -0.0001)) {
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else if (epsc <= -0.0001) {
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				} //  if (flag == -1)
                else if (flag == -2) { //%material never been in compression
                    if (epsc >= 0.0001){
                        dsigdh = 0;
					}
                    else if ((epsc > -0.0001) & (epsc < 0.0001)){ //%must define smooth transition curve
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= -0.0001){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
                    }
				} // elseif (flag == -2)
                else if (flag == -3) { //%temporary curve at beginning
                    if ((epsc >= 0) & (epsc < 0.0001)) {
                        if (depsc > 0 ) { //%still on temporary curve
                            dsigdh = 0; 
						}
                        else { //% depsc < 0 %smooth transition curve
                            epsr1 = epsc-depsc;
                            sigr1 = 0;
                            epsr2 = -epsr1;
                            depsr1dh = depsdhp;
                            dsigr1dh = 0;
                            depsr2dh = -depsdhp;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }


					} //if ((epsc >= 0) & (epsc < 0.0001))
                    else if ((epsc < 0) & (epsc > -0.0001)) {
                        //% flag = -3;
                        if (depsc < 0) { // %still on temporary curve
                            dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
						}
                        else{ // % depsc > 0 %smooth transition curve
                            epsr2 = epsc-depsc;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            epsr1 = -epsr2;
                            sigr1 = 0;
                            depsr2dh = depsdhp;
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh   = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            depsr1dh = -depsdhp;
                            dsigr1dh = 0;
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,0.001,-0.001, 0,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,0,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }
					}//else if ((epsc < 0) & (epsc > -0.0001))
                    else if (epsc >=0.0001) {
                        epsr1 = 0.0001;
                        sigr1 = 0;
                        epsr2 = -0.0001;
                        depsr1dh = 0;
                        depsr2dh = 0;
                        dsigr1dh = 0;
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        dsigdh = 0;
					}//else if (epsc >=0.0001)
                    else if (epsc <= -0.0001){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
    			
				}//else if (flag == -3) 
                else if (flag == 0) { // % monotonic envelope
                    if (depsc < 0) { // % keep loading on monotonic envelope
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else { //% depsc > 0 % unloading
                        epsr = (epsc-depsc);
                        sigr = sigp;

                        //epsp = Compute_epsp (epsr,sigr,eps0,fc,Ec);
						    // Modified Popovics-Saenz concrete stress-strain relationship: compute epsp

						if (fabs(epsr)< eps0)	epsp = -(fabs(epsr) - fabs(sigr)/Ec);
						else  epsp = -(eps0 - fc/Ec);
                      
						Eur = fabs(sigr/(epsr-epsp));
                        depsrdh = depsdhp;
                        dsigrdh = dsigdhp;
                        depspdh = Compute_depspdh(epsr,sigr,depsrdh,dsigrdh,deps0dh,dfcdh,dEcdh);
						double sign = 1.0;
						if ((sigr/(epsr-epsp))<0) sign = -1.0;
						else if ((sigr/(epsr-epsp))==0) sign = 0.0;

                        dEurdh = sign*(dsigrdh*(epsr-epsp)-sigr*(depsrdh-depspdh))/pow((epsr-epsp),2); 
                        
						if (epsc < epsr + eta*(epsp-epsr)) { // %temporary curve
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
						}
                        else {
                            epsr1 = epsr + eta*(epsp-epsr);
                            epsr2 = epsr - eta*(epsp-epsr);
                            sigr1 = Eur*(epsr1-epsp);
                            Monotonic_Envelope(epsr2,&sigr2,&Et2);
                            depsr1dh = depsrdh + eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                            depsr2dh = depsrdh - eta*(depspdh-depsrdh) - detadh*(epsp-epsr);
                            dsigr1dh = Eur*(depsr1dh-depspdh)+dEurdh*(epsr1-epsp);
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);    
                            if ((epsc >= epsr + eta*(epsp-epsr)) & (epsc <= epsp - eta*(epsp-epsr))){
                                dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
							}
                            else if ((epsc > epsp - eta*(epsp-epsr)) & (epsc < epsp + eta*(epsp-epsr))){
                                double dDeltadh = eta*(depspdh-depsrdh)+ detadh*(epsp-epsr);
                                dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
							}
                            else if (epsc >= epsp + eta*(epsp-epsr)){
                                dsigdh = 0;
							}
                            
                        }
                    } //else //% depsc > 0 % unloading
				} //else if (flag == 0) {   
                else if (flag == 1) {//% unloading curve with fixed end-points
                    if (epsc >= epsp + eta*(epsp-epsr)){
                        dsigdh = 0;
					}
                    else if ((epsc < epsp + eta*(epsp-epsr)) & (epsc > epsp - eta*(epsp-epsr))){
                        double dDeltadh = eta*(depspdh-depsrdh) + detadh*(epsp-epsr);                        
                        dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
					}
                    else if ((epsc <= epsp - eta*(epsp-epsr)) & (epsc >= epsr + eta*(epsp-epsr))){
                        dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
					}
                    else if ((epsc < epsr + eta*(epsp-epsr)) & (epsc > epsr - eta*(epsp-epsr))) {
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= epsr - eta*(epsp-epsr)){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				}//else if (flag == 1) 
                else if (flag == 2){ // % friction state
                    if (epsc >= epsp + eta*(epsp-epsr)){
                        dsigdh = 0;
					}
                    else if ((epsc < epsp + eta*(epsp-epsr)) & (epsc > epsp - eta*(epsp-epsr))){
                        double dDeltadh=eta*(depspdh-depsrdh) + detadh*(epsp-epsr); 
                        dsigdh = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
					}
                    else if ((epsc <= epsp - eta*(epsp-epsr)) & (epsc >= epsr + eta*(epsp-epsr))){
                        dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
					}
                    else if ((epsc < epsr + eta*(epsp-epsr)) & (epsc > epsr - eta*(epsp-epsr))){
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else if (epsc <= epsr - eta*(epsp-epsr)){
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    
				} //else if (flag == 2){    
                else if (flag == 3) { //% unloading curve with variable end-points
                    if ((epsc > epsr) & (epsc < epsr + eta*(epsp-epsr))){
                        if (depsc > 0){// %still on temporary curve
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp); 
						}
                        else {// % depsc < 0 %smooth transition curve
                            epsr1 = epsc-depsc;
                            sigr1 = Eur*(epsr1-epsp);
                            epsr2 = 2*epsr-epsr1;
                            Monotonic_Envelope (epsr2,&sigr2,&Et2);
                            depsr1dh = depsdhp;
                            dsigr1dh = dsigdhp;
                            depsr2dh = 2*depsrdh-depsr1dh;
                            dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                            dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                            dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        }
					} //if ((epsc > epsr) & (epsc < epsr + eta*(epsp-epsr)))
                    else if (epsc <= epsr - eta*(epsp-epsr)){ 
                        dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);
					}
                    else if ((epsc <= epsr) & (epsc > epsr - eta*(epsp-epsr))){
                        epsr1 = epsc-depsc;
                        sigr1 = Eur*(epsr1-epsp);
                        epsr2 = 2*epsr-epsr1;
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        depsr1dh = depsdhp;
                        dsigr1dh = dsigdhp;
                        depsr2dh = 2*depsrdh-depsr1dh;
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        dsigdh = Transition_r_sens(epsc,epsr1,epsr2,sigr1,sigr2,epsr+eta*(epsp-epsr),epsr-eta*(epsp-epsr),Eur,Et2,depsdh,depsr1dh,depsr2dh,dsigr1dh,dsigr2dh,dEurdh,dEt2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
					}
                    else { // % epsc >= epsr + eta*(epsp-epsr))
                        epsr1 = epsr + eta*(epsp-epsr);
                        sigr1 = Eur*(epsr1-epsp);
                        epsr2 = epsr - eta*(epsp-epsr);
                        Monotonic_Envelope (epsr2,&sigr2,&Et2);
                        depsr1dh = depsrdh + eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                        dsigr1dh = dEurdh*(epsr1-epsp)+Eur*(depsr1dh-depspdh);
                        depsr2dh = depsrdh - eta*(depspdh-depsrdh) - detadh*(epsp-epsr);
                        dsigr2dh = Monotonic_Envelope_sens (epsr2,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsr2dh);
                        dEt2dh = Monotonic_Envelope_Et_sens (epsr2,depsr2dh,dfcdh,deps0dh,depsudh,dfudh,dEcdh);
                        if (((epsc >= epsr + eta*(epsp-epsr))) & (epsc <= epsp - eta*(epsp-epsr))){
                            dsigdh = Eur*(depsdh-depspdh)+dEurdh*(epsc-epsp);
						}
                        else if ((epsc > epsp - eta*(epsp-epsr)) & (epsc < epsp + eta*(epsp-epsr))){
                            double dDeltadh = eta*(depspdh-depsrdh) + detadh*(epsp-epsr);
                            dsigdh   = Transition_p_sens (epsc,epsp,eta*(epsp-epsr),Eur,depsdh,depspdh,dDeltadh,dEurdh); 
						}
                        else if ((epsc >= epsp + eta*(epsp-epsr))){
                            dsigdh = 0;
                        }
                    }
				} //else if (flag == 3)  
        } //else { //%(depsc ~= 0)
		

				(*SHVs)(0,(gradNumber)) = depsdh;
				(*SHVs)(1,(gradNumber)) = dsigdh ;
				(*SHVs)(2,(gradNumber)) = depsrdh;
				(*SHVs)(3,(gradNumber)) = dsigrdh;
				(*SHVs)(4,(gradNumber)) = depspdh;
				(*SHVs)(5,(gradNumber)) = dEurdh;
				(*SHVs)(6,(gradNumber)) = depsr1dh;
				(*SHVs)(7,(gradNumber)) = depsr2dh;
				(*SHVs)(8,(gradNumber)) = dsigr1dh;
				(*SHVs)(9,(gradNumber)) = dsigr2dh;
				(*SHVs)(10,(gradNumber)) = dEt2dh;

	return 0;
}

double SmoothPSConcrete::getStrainSensitivity (int gradNumber) {
	if (SHVs==0) {
		opserr<<"warning recordering SmoothPSConcrete::getStrainSensitivity! SHVs=0";
		return 0.0;
	}  
	return (*SHVs)(0,(gradNumber));
}



int SmoothPSConcrete::Monotonic_Envelope (double epsc, double * sig, double * Et){
    // Modified Popovics-Saenz concrete stress-strain relationship: monotonic envelope
    
 
  
    double K = Ec*eps0/fc;
    
    if (epsc > -eps0)  {       // ascending branch in compression (Popovics)
        double D = K-1; 
        double r = K/D;
        double Eta = epsc/(-eps0);
        (*sig) = -fc*K*Eta/(1+D*(pow(Eta,r)));
        (*Et)  = -(1./eps0)*fc*K*(-1.-D*pow(Eta,r) + D*pow(Eta,r)*r)/pow((1+D*pow(Eta,r)),2);
	}
    else if ((epsc <= -eps0) && (epsc>-epsmax)){        // descending branch in compression (Saenz)
        double Eta = epsc/(-eps0);
        double Ksig = fc/fcu;
        double Keps = epsu/eps0;
        double C = K * (Ksig -1)/pow((Keps - 1),2) - 1./Keps;
        double A = C+K-2;
        double B = 1-2*C;
        (*sig) = -fc*K*Eta/(1+A*Eta+B*Eta*Eta+C*pow(Eta,3));
        (*Et)  = -(1/eps0)*fc*K*(-1+B*Eta*Eta+2*C*pow(Eta,3))/pow((1+A*Eta+B*Eta*Eta+C*Eta*Eta*Eta),2);
	}
	else {
		(*sig) = sigmax; 
		(*Et) = 0;
	}

	return 0;
}
    
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int SmoothPSConcrete::Compute_epsp (void){
    // Modified Popovics-Saenz concrete stress-strain relationship: compute epsp
    double epsr = fabs(Tepsr);
    if (epsr < eps0)
        this->Tepsp = -(epsr - fabs(Tsigr)/Ec);
    else
        this->Tepsp = -(eps0 - fc/Ec);
	return 0;
}
    
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int SmoothPSConcrete::Transition_r(double epsc, double e1,double e2,double s1,double s2,double e1th,double e2th,double Et1,double Et2, double * sig, double * Et){

    // e2th < epsc < e1th

	int result=0;

    double b = (-(Et2+2*Et1)*(e2-e1)+3*(s2-s1))/((e2-e1)*(e2-e1));
    double a = (Et2-Et1-2*b*(e2-e1))/3./((e2-e1)*(e2-e1));
    //c = Et1;
    if (fabs(e1-e1th)<1e-16){ //e1th=e1 & e2th=e2
        this->Tsig = a*pow((epsc-e1),3) + b*(epsc-e1)*(epsc-e1)+ Et1*(epsc-e1) + s1;
        this->TEt  = 3*a*(epsc-e1)*(epsc-e1) + 2*b*(epsc-e1) + Et1;
	}
    else {
        if (epsc >= e1) { // & epsc < e1th
            (*sig) = s1 + Et1*(epsc-e1);
            (*Et) = Et1;
		}
        else if ((epsc > e2) && (epsc < e1)) {
            (*sig) = a*pow((epsc-e1),3) + b*(epsc-e1)*(epsc-e1) + Et1*(epsc-e1) + s1;
            (*Et)  = 3*a*(epsc-e1)*(epsc-e1)+ 2*b*(epsc-e1) + Et1;
		}
        else // epsc<=e2 & epsc > e2th
            result=Monotonic_Envelope (epsc, sig, Et);
	}
	return result;
}
            
        
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int SmoothPSConcrete::Transition_p (double Delta){
    //Et1 = 0; Et2 = Eur; s1 = 0; c = Et1;
        
     double e1 = Tepsp + Delta;
     double e2 = Tepsp - Delta;
     double s2 = - Delta*TEur;
     double a  = (TEur*(e2-e1)-2*s2)/pow((e2-e1),3);
     double b  = (-TEur*(e2-e1)+3*s2)/pow((e2-e1),2); 

     Tsig = a*pow((Tepsc-e1),3) + b*(Tepsc-e1)*(Tepsc-e1);
     TEt  = 3*a*(Tepsc-e1)*(Tepsc-e1) + 2*b*(Tepsc-e1);
	 return 0;
}
 
double  SmoothPSConcrete::Monotonic_Envelope_sens(double epsilon,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh,double depsdh){
    // Modified Popovics-Saenz concrete stress-strain relationship: monotonic envelope 
    
  
    double K = Ec*eps0/fc;
    double dsigdh=0;
    if (epsilon >- eps0) { //         % ascending branch in compression (Popovics)
        double D = K-1; 
        double r = K/D;
        double Eta = epsilon/(-eps0);
        double Et  = -(1/eps0)*fc*K*(-1-D*pow(Eta,r) + D*pow(Eta,r)*r)/(pow((1+D*pow(Eta,r)),2));
        
        dsigdh = ( Ec*Ec*epsilon*pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*
            eps0*(Ec*eps0-fc-log(-epsilon/eps0)*fc)/(pow((fc+pow(-epsilon/eps0,
            Ec*eps0/(Ec*eps0-fc))*Ec*eps0-pow(-epsilon/eps0,Ec*eps0/
            (Ec*eps0-fc))*fc),2)*(Ec*eps0-fc)) )*dfcdh+
            (log(-epsilon/eps0)*pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*
            fc*fc*epsilon*Ec*Ec/(pow((fc+pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*
            Ec*eps0-pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*fc),2)*(Ec*eps0-fc)) )*deps0dh+
            (epsilon*fc*fc*(Ec*eps0-fc-pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*
            Ec*eps0+pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*fc+pow(-epsilon/eps0,
            Ec*eps0/(Ec*eps0-fc))*Ec*eps0*log(-epsilon/eps0))/(pow((fc+
            pow(-epsilon/eps0,Ec*eps0/(Ec*eps0-fc))*Ec*eps0-pow(-epsilon/eps0
            ,Ec*eps0/(Ec*eps0-fc))*fc),2)*(Ec*eps0-fc)) )*dEcdh+Et*depsdh;
	}
        
    else if ((epsilon <= -eps0) && (epsilon>-epsmax)){        // descending branch in compression (Saenz)
        double Eta = epsilon/(-eps0);
        double Ksig = fc/fcu;
        double Keps = epsu/eps0;
        double C = K * (Ksig -1)/pow((Keps - 1),2) - 1./Keps;
        double A = C+K-2;
        double B = 1-2*C;
        double Et  = -(1./eps0)*fc*K*(-1+B*pow(Eta,2)+2*C*pow(Eta,3))/pow((1+A*Eta+B*Eta*Eta+C*pow(Eta,3)),2);
         
		dsigdh = (-Ec*epsilon*(Ec*epsilon/(fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1))+
            2*Ec*epsilon*epsilon/(eps0*fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1))+
            Ec*epsilon*epsilon*epsilon/(eps0*eps0*fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1)))/
            pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*
            epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/(eps0*eps0)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*epsilon*epsilon*epsilon/pow(eps0,3)),2) )*dfudh+
            (-Ec*epsilon*(-(-Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))+Ec*eps0/(fc*fcu*(epsu/eps0-1)*(epsu/eps0-1))
            -Ec*eps0/(fc*fc))*epsilon/eps0+(2*Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))-2*Ec*eps0/
            (fc*fcu*(epsu/eps0-1)*(epsu/eps0-1)))*epsilon*epsilon/(eps0*eps0)-(-Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))+
            Ec*eps0/(fc*fcu*(epsu/eps0-1)*(epsu/eps0-1)))*pow(epsilon,3)/pow(eps0,3))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))
            -eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/(eps0*eps0)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*epsilon*epsilon*epsilon/pow(eps0,3)),2))*dfcdh+
            (-Ec*epsilon*(-(-2*Ec*(fc/fcu-1)/(fc*pow((epsu/eps0-1),3))+eps0/(epsu*epsu))*epsilon/eps0+(4*Ec*(fc/fcu-1)/
            (fc*pow((epsu/eps0-1),3))-2*eps0/(epsu*epsu))*epsilon*epsilon/(eps0*eps0)-(-2*Ec*(fc/fcu-1)/(fc*pow((epsu/eps0-1),3))
            +eps0/(epsu*epsu))*pow(epsilon,3)/pow(eps0,3))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)
            *epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/pow(eps0,2)-
            (Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*depsudh+
            (-Ec*epsilon*(-(Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))
            -1/epsu+Ec/fc)*epsilon/eps0+(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*
            epsilon/pow(eps0,2)+(-2*Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-4*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))+
            2*1/epsu)*epsilon*epsilon/pow(eps0,2)-2*(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/pow(eps0,3)
            -(Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))-1/epsu)*pow(epsilon,3)/pow(eps0,3)+
            3*(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,4))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))
            -eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/(eps0*eps0)
            -(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*deps0dh+
            (epsilon/(1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+
            (1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/(eps0*eps0)-
            (Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3))-Ec*epsilon*(-(eps0*(fc/fcu-1)/
            (fc*(epsu/eps0-1)*(epsu/eps0-1))+eps0/fc)*epsilon/eps0-2*(fc/fcu-1)*epsilon*epsilon/(eps0*fc*(epsu/eps0-1)*(epsu/eps0-1))-
            (fc/fcu-1)*pow(epsilon,3)/(eps0*eps0*fc*(epsu/eps0-1)*(epsu/eps0-1)))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-
            eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/pow(eps0,2)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*dEcdh+
            Et*depsdh;

	}
	else {
		epsilon = -epsmax;
	    double Eta = epsilon/(-eps0);
        double Ksig = fc/fcu;
        double Keps = epsu/eps0;
        double C = K * (Ksig -1)/pow((Keps - 1),2) - 1./Keps;
        double A = C+K-2;
        double B = 1-2*C;
        //double Et  = -(1./eps0)*fc*K*(-1+B*pow(Eta,2)+2*C*pow(Eta,3))/pow((1+A*Eta+B*Eta*Eta+C*pow(Eta,3)),2);
         
		dsigdh = (-Ec*epsilon*(Ec*epsilon/(fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1))+
            2*Ec*epsilon*epsilon/(eps0*fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1))+
            Ec*epsilon*epsilon*epsilon/(eps0*eps0*fcu*fcu*(epsu/eps0-1)*(epsu/eps0-1)))/
            pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*
            epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/(eps0*eps0)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*epsilon*epsilon*epsilon/pow(eps0,3)),2) )*dfudh+
            (-Ec*epsilon*(-(-Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))+Ec*eps0/(fc*fcu*(epsu/eps0-1)*(epsu/eps0-1))
            -Ec*eps0/(fc*fc))*epsilon/eps0+(2*Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))-2*Ec*eps0/
            (fc*fcu*(epsu/eps0-1)*(epsu/eps0-1)))*epsilon*epsilon/(eps0*eps0)-(-Ec*eps0*(fc/fcu-1)/(fc*fc*(epsu/eps0-1)*(epsu/eps0-1))+
            Ec*eps0/(fc*fcu*(epsu/eps0-1)*(epsu/eps0-1)))*pow(epsilon,3)/pow(eps0,3))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))
            -eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/(eps0*eps0)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*epsilon*epsilon*epsilon/pow(eps0,3)),2))*dfcdh+
            (-Ec*epsilon*(-(-2*Ec*(fc/fcu-1)/(fc*pow((epsu/eps0-1),3))+eps0/(epsu*epsu))*epsilon/eps0+(4*Ec*(fc/fcu-1)/
            (fc*pow((epsu/eps0-1),3))-2*eps0/(epsu*epsu))*epsilon*epsilon/(eps0*eps0)-(-2*Ec*(fc/fcu-1)/(fc*pow((epsu/eps0-1),3))
            +eps0/(epsu*epsu))*pow(epsilon,3)/pow(eps0,3))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)
            *epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/pow(eps0,2)-
            (Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*depsudh+
            (-Ec*epsilon*(-(Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))
            -1/epsu+Ec/fc)*epsilon/eps0+(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*
            epsilon/pow(eps0,2)+(-2*Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-4*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))+
            2*1/epsu)*epsilon*epsilon/pow(eps0,2)-2*(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/pow(eps0,3)
            -(Ec*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*Ec*(fc/fcu-1)*epsu/(eps0*fc*pow((epsu/eps0-1),3))-1/epsu)*pow(epsilon,3)/pow(eps0,3)+
            3*(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,4))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))
            -eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/(eps0*eps0)
            -(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*deps0dh+
            (epsilon/(1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+
            (1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*epsilon*epsilon/(eps0*eps0)-
            (Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3))-Ec*epsilon*(-(eps0*(fc/fcu-1)/
            (fc*(epsu/eps0-1)*(epsu/eps0-1))+eps0/fc)*epsilon/eps0-2*(fc/fcu-1)*epsilon*epsilon/(eps0*fc*(epsu/eps0-1)*(epsu/eps0-1))-
            (fc/fcu-1)*pow(epsilon,3)/(eps0*eps0*fc*(epsu/eps0-1)*(epsu/eps0-1)))/pow((1-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-
            eps0/epsu+Ec*eps0/fc-2)*epsilon/eps0+(1-2*Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))+2*eps0/epsu)*
            epsilon*epsilon/pow(eps0,2)-(Ec*eps0*(fc/fcu-1)/(fc*(epsu/eps0-1)*(epsu/eps0-1))-eps0/epsu)*pow(epsilon,3)/pow(eps0,3)),2))*dEcdh;
			

	
	}

	return dsigdh;
};


double SmoothPSConcrete::Compute_depspdh (double epsr,double sigr,double depsrdh,double dsigrdh,double deps0dh,double dfcdh,double dEcdh){
    // % Modified Popovics-Saenz concrete stress-strain law: compute depspdh
    double az=1.0;
	if (epsr<0) az=-1.0;
	else if (epsr==0) az=0.0;
	
	double depspdh =0.0;

	epsr = fabs(epsr);
	if (epsr < eps0)  {
		double sign = 1.0;
		if 	(sigr<0)   sign = -1.0;
		else if (sigr==0)   sign = 0.0;
			
		depspdh = -(az*depsrdh - sign*dsigrdh/Ec + fabs(sigr)*dEcdh/(Ec*Ec));
	}
	else {
		depspdh = -(deps0dh - dfcdh/Ec + fc*dEcdh/(Ec*Ec));
	}
	return depspdh;
}

double SmoothPSConcrete::Transition_r_sens(double epsc,double e1,double e2,double s1,double s2,double e1th,double e2th,double Et1,double Et2,
										   double depsdh,double de1dh,double de2dh,double ds1dh,double ds2dh,double dEt1dh,double dEt2dh,
										   double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh){
        
        double b = (-(Et2+2*Et1)*(e2-e1)+3*(s2-s1))/pow((e2-e1),2);
        double a = 1./3*(Et2-Et1-2*b*(e2-e1))/pow((e2-e1),2);

        double dbdh = (3*(ds2dh-ds1dh)-(dEt2dh+2*dEt1dh)*(e2-e1)-(Et2+2*Et1)*(de2dh-de1dh))/pow((e2-e1),2) 
			-2*((-(Et2+2*Et1)*(e2-e1)+3*(s2-s1))/pow((e2-e1),3))*(de2dh-de1dh);
        double dadh = 1./3*((dEt2dh-dEt1dh-2*dbdh*(e2-e1)-2*b*(de2dh-de1dh))/pow((e2-e1),2) -2*(Et2-Et1-2*b*(e2-e1))/pow((e2-e1),3)*(de2dh-de1dh));
        // %c = Et1;
		double dsigdh = 0.0;

        if (fabs(e1-e1th)<1e-16) { //%e1th=e1 & e2th=e2
            double Et     = 3*a*(epsc-e1)*(epsc-e1) + 2*b*(epsc-e1) + Et1;
            dsigdh = dadh*pow((epsc-e1),3) + dbdh*(epsc-e1)*(epsc-e1) + dEt1dh*(epsc-e1) + ds1dh + Et*(depsdh-de1dh); 
		}
        else {
            if (epsc >= e1) {  // % & epsc < e1th
                dsigdh = ds1dh + dEt1dh*(epsc-e1) + Et1*(depsdh-de1dh);
			}
            else if ((epsc > e2) & (epsc < e1)){
                double Et     = 3*a*(epsc-e1)*(epsc-e1) + 2*b*(epsc-e1) + Et1;
                dsigdh = dadh*pow((epsc-e1),3) + dbdh*(epsc-e1)*(epsc-e1) + dEt1dh*(epsc-e1) + ds1dh + Et*(depsdh-de1dh);
			}
            else { //% epsc<=e2 & epsc > e2th
                dsigdh = Monotonic_Envelope_sens (epsc,dfcdh,deps0dh,depsudh,dfudh,dEcdh,depsdh);    
			}
            
        }
		return dsigdh;
		
}

double SmoothPSConcrete::Transition_p_sens(double epsc,double epsp,double Delta,double Eur,double depsdh,double depspdh,double dDeltadh,double dEurdh){
      //  %Et1 = 0; Et2 = Eur; s1 = 0;
        
        double e1 = epsp + Delta;
        double de1dh = depspdh + dDeltadh;
        double e2 = epsp - Delta;
        double de2dh = depspdh - dDeltadh;
        double s2 = - Delta*Eur;
        double ds2dh = - Delta*dEurdh - dDeltadh*Eur;
        double a = (Eur*(e2-e1)-2*s2)/pow((e2-e1),3); 
        double b = (-Eur*(e2-e1)+3*s2)/pow((e2-e1),2); 
        double dadh = (dEurdh*(e2-e1)+Eur*(de2dh-de1dh)-2*ds2dh)/pow((e2-e1),3) - 3*((Eur*(e2-e1)-2*s2)/pow((e2-e1),4))*(de2dh-de1dh);
        double dbdh = (-dEurdh*(e2-e1)-Eur*(de2dh-de1dh)+3*ds2dh)/pow((e2-e1),2) - 2*((-Eur*(e2-e1)+3*s2)/pow((e2-e1),3))*(de2dh-de1dh);
        double Et  = 3*a*(epsc-e1)*(epsc-e1)+ 2*b*(epsc-e1);
        
        double dsigdh = dadh*pow((epsc-e1),3) + dbdh*(epsc-e1)*(epsc-e1) + Et*(depsdh-de1dh);

		return dsigdh;
}
     
double SmoothPSConcrete::Monotonic_Envelope_Et_sens (double epsilon,double depsdh,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh){

	//    % Modified Popovics-Saenz concrete stress-strain relationship: monotonic envelope
    
      
    
    double K = Ec*eps0/fc;
    double dKdh = dEcdh*eps0/fc + Ec*deps0dh/fc - Ec*eps0/pow(fc,2)*dfcdh;
    double dEtdh =0.0;

    if (epsilon > -eps0)  {      // % ascending branch in compression (Popovics)
        double D = K-1;
        double r = K/D;
        double dDdh = dKdh;
        double Eta = epsilon/(-eps0);
        double drdh = -dKdh/pow(D,2);
        double dEtadh = -depsdh/eps0 + epsilon/pow(eps0,2)*deps0dh;
        double t1=(-1-D*pow(Eta,r) + D*pow(Eta,r)*r);
        double t2=pow((1+D*pow(Eta,r)),2);
        double dt1=pow(Eta,r)*(drdh*log(Eta)+r/Eta*dEtadh);
        double dt2=2*(1+D*pow(Eta,r))*pow(Eta,r)*(dDdh+drdh*log(Eta)+r/Eta*dEtadh);

        dEtdh=deps0dh/pow(eps0,2)*fc*K*t1/t2 - K/eps0*dfcdh*t1/t2 - fc/eps0*dKdh*t1/t2 - fc*K/eps0*dt1/t2 + fc*K/eps0*t1/pow(t2,2)*dt2;
	}
    else if ((epsilon <= -eps0) && (epsilon>-epsmax)){        // descending branch in compression (Saenz)
        double Eta = epsilon/(-eps0);
        double dEtadh = -depsdh/eps0 + epsilon/pow(eps0,2)*deps0dh;
        double Ksig = fc/fcu;
        double dKsig = dfcdh/fcu -fc/pow(fcu,2)*dfudh;
        double Keps = epsu/eps0;
        double dKeps = depsudh/eps0 - epsu/pow(eps0,2)*deps0dh;
        double C = K * (Ksig -1)/pow((Keps - 1),2) - 1/Keps;
        double dCdh = dKdh*(Ksig -1)/pow((Keps - 1),2) + dKsig*K/pow((Keps - 1),2) - 2*K*(Ksig-1)*dKeps/pow((Keps - 1),3) +dKeps/pow(Keps,2);
        double A = C+K-2;
        double dAdh = dCdh + dKdh;
        double B = 1-2*C;
        double dBdh = -2*dCdh;
        double t1 = (-1+B*Eta*Eta+2*C*Eta*Eta*Eta);
        double dt1 = dBdh*Eta*Eta+2*B*Eta*dEtadh+2*dCdh*Eta*Eta*Eta+6*C*Eta*Eta*dEtadh;
        double t2 = pow((1+A*Eta+B*Eta*Eta+C*Eta*Eta*Eta),2);
        double dt2 = 2*(1+A*Eta+B*Eta*Eta+C*Eta*Eta*Eta)*(dAdh*Eta+A*dEtadh+dBdh*Eta*Eta+2*B*Eta*dEtadh+dCdh*Eta*Eta*Eta+3*C*Eta*Eta*dEtadh);
        
        dEtdh=deps0dh/pow(eps0,2)*fc*K*t1/t2 - K/eps0*dfcdh*t1/t2 - fc/eps0*dKdh*t1/t2 - fc*K/eps0*dt1/t2 + fc*K/eps0*t1/pow(t2,2)*dt2;
    }
	else {
		dEtdh= 0.0;
	}

	return dEtdh;
}

        


int  SmoothPSConcrete::compute_epsmax(double * epsmax,double * sigmax){
	double K = Ec*eps0/fc;
	double Ksig = fc/fcu;
	double Keps = epsu/eps0;
	double C = K * (Ksig -1)/pow((Keps - 1),2) - 1/Keps;
	double A = C+K-2;
	double B = 1-2*C;

	double a2=.5*B/C;
	double a0=-.5/C;
	double Q = -a2*a2/9;
	double R = -(27*a0+2*a2*a2*a2)/54;
	double D = pow(Q,3)+R*R;

	double S,T;

	double z=1000;

	if (D<0){
		double a = sqrt(R*R+fabs(D));
		double sita = atan(sqrt(fabs(D))/R);
		S = 2*pow(a,1./3)*cos(1./3*sita);
		T = 2*pow(a,1./3)*sin(1./3*sita);

		double z1[3];
		z1[0]=-a2/3+S;
		z1[1]=-a2/3-.5*(S)+.5*pow(3,.5)*(T);
		z1[2]=-a2/3-.5*(S)-.5*pow(3,.5)*(T);

		z=1.0;
		
		int maxNum=0, minNum=0,medNum=0;
		
		for (int i=0;i<3;i++){
			if (z1[i]>z1[maxNum]) maxNum=i;
			if (z1[i]<z1[minNum]) minNum =i;
			
		}
		
		for (int i=0;i<3;i++) { if ((i!=maxNum) &(i!=minNum)) medNum=i;}
		


		if (z1[maxNum]<1+1.e-14){ 
			opserr<<"wrong parameter in SmoothPSConcrete::compute_epsmax!"<<endln;
			//exit(-1);//Michele December 4 2006: test debug
			z=1.0;
		}
    
		else if (fabs(z1[medNum]-1)<1.e-14) z=z1[maxNum];
		else z = z1[medNum];
	};

	(*epsmax)=z*eps0;
	double emax=(*epsmax)/eps0;
	(*sigmax) = -fc*K*emax/(1+A*emax+B*emax*emax+C*pow(emax,3));
	return 0;
}




double SmoothPSConcrete::getInitialTangentSensitivity(int gradNumber)	{
		// For now, assume that this is only called for initial stiffness 
		if (parameterID == 5) {
			return 1.0; 
		}
		else {
			return 0.0;
		}
	}
