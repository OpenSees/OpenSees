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
                                                                        

// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2=6th order Adams-Bashforth-Moulton, 3=modified Rosenbrock Triple)
// $RelTol:	Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol:	Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-6)
// $MaxHalf: Maximum number of sub-step iterations within an integration step, h=dt*(0.5)^MaxHalf (default value 15)


#include <math.h>
#include <elementAPI.h>
#include <APDVFD.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <Parameter.h>
#define PI acos(-1)
static int numAPDVFD = 0;

void *
OPS_APDVFD(void)

{
  if (numAPDVFD == 0) {
    numAPDVFD++;
    opserr << "APDVFD Model by BUCEA\n";
  }

  opserr << "Due to known issues and unreliable results, this material has been" << endln;
  opserr << "temporarily removed from the compiled versions of OpenSees (Tcl and Py)" << endln;
  opserr << "The material source code remains available. Compile at your own risk." << endln;
  return 0;
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[21];
  int numData = 1;
        // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  APDVFD tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 16 && numData != 17 && numData != 21) {
    opserr << "Invalid #args, want: uniaxialMaterial APDVFD " << iData[0] << "K?G1?G2?Alpha?L?LC?DP?DG?N1?N2?DO1?DO2?DC?S?HP?HC?<LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
       
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args want: uniaxialMaterial APDVFD " << iData[0] << "K?G1?G2?Alpha?L?LC?DP?DG?N1?N2?DO1?DO2?DC?S?HP?HC?<LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    
    return 0;   
  }
  
  if (numData == 16) {
    // Default variables
    dData[16] = 0.0;
    dData[17] = 1;
    dData[18] = 0.000001;
    dData[19] = 0.0000000001;
    dData[20] = 15;
  }
  
  if (numData == 17) {
    // Default variables
    dData[17] = 1;
    dData[18] = 0.000001;
    dData[19] = 0.0000000001;
    dData[20] = 15;
  }
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new APDVFD(iData[0], 
                           dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20]);
  
                                 
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type APDVFD Material\n";
    return 0;
  }
  
  return theMaterial;
}

APDVFD::APDVFD(int tag, double k, double g1, double g2, double a, double l, double lc, double dp, double dg, double n1, double n2, double do1, double do2, double dc, double s, double hp, double hc, double lgap, double nm, double reltol, double abstol, double maxhalf)
:UniaxialMaterial(tag,MAT_TAG_APDVFD), K(k), G1(g1), G2(g2), Alpha(a), L(l), LC(lc), DP(dp), DG(dg), N1(n1), N2(n2), DO1(do1), DO2(do2), DC(dc), S(s), HP(hp), HC(hc), LGap(lgap), NM(nm), RelTol(reltol), AbsTol(abstol), MaxHalf(maxhalf)
{
    if (Alpha < 0.0) {
      opserr << "APDVFD::APDVFD -- Alpha < 0.0, setting to 1.0\n";
      Alpha = 1.0;
    }
    
        
    //initialize variables
    this->revertToStart();
}


APDVFD::APDVFD()
:UniaxialMaterial(0,MAT_TAG_APDVFD),
 K(0.0), G1(0.0), G2(0.0), Alpha(0.0), L(0.0), LC(0.0), DP(0.0), DG(0.0), N1(0.0), N2(0.0), DO1(0.0), DO2(0.0), DC(0.0), S(0.0), HP(0.0), HC(0.0), LGap(0.0),  NM(0.0), RelTol(0.0), AbsTol(0.0), MaxHalf(0.0)
{
        this->revertToStart();
}

APDVFD::~APDVFD()
{
  // does nothing
}

int
APDVFD::setTrialStrain(double strain, double strainRate)
{
    //all variables to the last commit
    double C;
    double Alpha;
   


    this->revertToLastCommit();

 
    // Determine the strain rate and acceleration 

    double Vel, fd0, acc, vel1, vel0;
    if (fabs(strainRate) == 0.0) { //static analysis
        Vel = 0.0;
        acc = 0.0;

    }
    else {
        Vel = strainRate;
        acc = (Vel - TVel) / ops_Dt;
    }

    double smin = pow(0.5, MaxHalf);
    double s = 1.0;
    double stot = 0.0;
    double it = 0.0;
    fd0 = Tstress;

    double h, yt, eps, error;
    vel0 = TVel;  // Velocity of the previous step.


    while (it < 1.0) { //iteration
        h = s * ops_Dt; // Time step 
        vel1 = vel0 + acc * h; // Velocity at the time step h


        // Selection of Numerical Method to solve the ODE
        if (NM == 1.0) {
            DormandPrince(vel0, vel1, fd0, h, yt, eps, error);
        }
        if (NM == 2.0) {
            ABM6(vel0, vel1, fd0, h, yt, eps, error);
        }
        if (NM == 3.0) {
            ROS(vel0, vel1, fd0, h, yt, eps, error);
        }

        // Error check: Adaptive Step Size
        if ((eps <= RelTol) || (s == smin) || (fabs(error) <= AbsTol)) {
            vel0 = vel1;
            fd0 = yt;
            stot = stot + s;
        }
        else {
            if (s > smin) {
                s = 0.5 * s; // step gets smaller -now try this step again.
            }
            else {
                s = smin;
            }
        }

        if (stot == 1.0) { // The total internal stepsize reached dt
            it = 1.0;
        }
    }

    // Effect of gap start 

    if (LGap > 0.) {

        double dStrain = (strain - Tstrain);

        if ((fd0 > 0) && (Tstress < 0)) {  //from negtive to positive
            Tpugr = Tstrain + dStrain * fabs(fd0) / fabs(fd0 - Tstress);  // aproximate displacement for gap initiation
            Tnugr = 0.;

            if (fabs(strain - Tpugr) < LGap) {
                fd0 = 0.;
            }
        }

        if ((fd0 < 0) && (Tstress > 0)) {  //from positive to negative

            Tnugr = Tstrain + dStrain * fabs(fd0) / fabs(fd0 - Tstress);  // aproximate displacement for gap initiation
            Tpugr = 0.;

            if (fabs(strain - Tnugr) < LGap) {
                fd0 = 0.;
            }
        }

        // After gap inititon

        if ((fabs(Tpugr) > 0.) && (Tstress == 0)) {   //from negtive to positive

            if ((strain > Tpugr) && ((strain - Tpugr) < LGap)) {
                fd0 = 0.;
            }
        }



        if ((fabs(Tnugr) > 0.) && (Tstress == 0)) {   //from positive to negative

            if ((strain < Tnugr) && ((strain - Tnugr) > -LGap)) {
                fd0 = 0.;
            }
        }

    }
    // Effect of gap end 


    Tstress = fd0; // Stress 
    TVel = Vel;
    Tstrain = strain;

    return 0;
}

double APDVFD::getStress(void)
{
  return  Tstress;
}

double APDVFD::getTangent(void)
{
  return 0;
}

double APDVFD::getInitialTangent(void)
{
  return 0;
}

double APDVFD::getDampTangent(void)
{
  
  return 0; 
}


double 
APDVFD::getStrain(void)
{
    return Tstrain;
}

double 
APDVFD::getStrainRate(void)
{
  return TVel;
}

int 
APDVFD::commitState(void)
{
  //commit trial  variables
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  CVel = TVel;
  Cpugr = Tpugr;
  Cnugr = Tnugr; 
  
  return 0;
}

int 
APDVFD::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  TVel = CVel;
  Tpugr = Cpugr;
  Tnugr = Cnugr;
  
  return 0;
}

int 
APDVFD::revertToStart(void)
{
  // Initialize state variables
  Tstrain=0.0;
  Tstress=0.0;
  Ttangent = 0.0;
  TVel = 0.0;
  Tpugr = 0.0;
  Tnugr = 0.0;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Ctangent = 0.0;
  CVel = 0.0;
  Cpugr = 0.0;
  Cnugr = 0.0;
  
  return 0;
}

double
APDVFD::sgn(double dVariable){
  if (dVariable<0.0){
    return -1.0;
  }else{
    return 1.0;
  }
}

int 
APDVFD::DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error){
  double k1, k2, k3, k4, k5, k6, k7;
 //int f(int vel0, int y0);
  k1 = f(vel0, y0) * h;

  k2 = f((vel1 - vel0)*(1./5.) + vel0, y0 + (1./5.)*k1) * h;

  k3 = f((vel1 - vel0)*(3./10.) + vel0, y0 + (3./40.)*k1 + (9./40.)*k2) * h;

  k4 = f((vel1 - vel0)*(4./5.) + vel0, y0 + (44./45.)*k1 + (-56./15.)*k2 + (32./9.)*k3) * h;

  k5 = f((vel1 - vel0)*(8./9.) + vel0, y0 + (19372.0/6561.0)*k1 + (-25360.0/2187.0)*k2 + (64448.0/6561.0)*k3 + (-212.0/729.0)*k4) * h;

  k6 = f((vel1 - vel0)*(1.) + vel0, y0 + (9017.0/3168.0)*k1 + (-355.0/33.0)*k2 + (46732.0/5247.0)*k3 + (49.0/176.0)*k4 + (-5103.0/18656.0)*k5) * h;

  yt = y0 + (35./384.)*k1 + (500./1113.)*k3 + (125./192.)*k4 + (-2187./6784.)*k5 + (11./84.)*k6;

  k7 = f((vel1 - vel0)*(1.) + vel0, yt) * h;

  error = (71./57600.)*k1 + (-71./16695.)*k3 + (71./1920.)*k4 + (-17253./339200.)*k5 + (22./525.)*k6 + (-1./40.)*k7;

  eps = fabs(error/ yt);


 return 0;
}

int 
APDVFD::ABM6(double vel0, double vel1, double y0, double h, double& y6, double& eps, double& error){
 double f0, f1, f2, f3, f4, f5, f6, y11, y2, y3, y4, y5, yp6;
 //int f(int vel0, int y0);
    h = h/6.0;

    f0 = f((vel1 - vel0)*(0./6.) + vel0, y0);

    y11 = y0 + h*f0; //predictor

    f1 = f((vel1 - vel0)*(1./6.) + vel0, y11);

    y11 = y0 + h*f1; //corrector

    y2 = y11 + 0.5*h*(3.*f1 - 1.*f0); //predictor

    f2 = f((vel1 - vel0)*(2./6.) + vel0, y2);

    y2 = y11 + 0.5*h*(f2 + f1);  //corrector

    y3 = y2 + h/12.*(23.*f2 - 16.*f1 + 5.*f0); //predictor

    f3 = f((vel1 - vel0)*(3./6.) + vel0, y3);

    y3 = y2 + h/12.*(5.*f3 + 8.*f2 - 1.*f1); //corrector

    y4 = y3 + h/24.*(55.*f3 - 59.*f2 + 37.*f1 - 9.*f0); //predictor

    f4 = f((vel1 - vel0)*(4./6.) + vel0, y4);

    y4 = y3 + h/24.*(9.*f4 + 19.*f3 -5.*f2 + f1); // corrector

    y5 = y4 + h/720.*(1901.*f4 - 2774.*f3 + 2616.*f2 - 1274.*f1 + 251.*f0); // predictor

    f5 = f((vel1 - vel0)*(5./6.) + vel0, y5);

    y5 = y4 + h/720.*(251.*f5 + 646.*f4 -264.*f3 + 106.*f2 -19.*f1); // corrector

    yp6 = y5 + h/1440.*(4277.*f5 -7923.*f4 +9982.*f3 -7298.*f2 + 2877.*f1 -475.*f0); // predictor

    f6 = f((vel1 - vel0)*(6./6.) + vel0, yp6);

    y6 = y5 + (h/1440.)*(475.*f6 +1427.*f5  -798.*f4 + 482.*f3 -173.*f2 + 27.*f1); // corrector

	error = (yp6-y6);

    eps = fabs(error/y6);


 return 0;
}


//void itoa(int n, char* s); 
int
APDVFD::ROS(double vel0, double vel1, double y0, double h, double& y2, double& eps, double& error){
  double k1, k2, k3, J, T, d, e32, W, f0, f1, f2, y3;
  double Alpha;
  double C;
  //int f(int vel0, int y0);
     J = -K / (Alpha*C);
    T = K;
    d = 1. / (2. + sqrt(2.));
    e32 = 6. + sqrt(2.);
    W = 1. - h*d*J;
    f0 = f(vel0, y0);
    k1 = (f0 + h*d*T)/W;
    f1 = f((vel1 - vel0)*(0.5) + vel0, y0 + (0.5)*k1*h);
    k2 = (f1 - k1)/W + k1;
    y2  = y0 + h*k2;
    f2 = f(vel1, y2);
    k3 = 1./W *(f2 - e32*(k2 - f1) - 2.*(k1 - f0) + h*d*T);
    error = h/6.*(k1 -2.*k2 + k3);
    y3 = y2 + error;
    eps = fabs(error/(y3));

 return 0;
}


double
APDVFD::f(double v, double fd) {


    void setTrialStrain(double strain);
    double strain{};

    double i = G1 * L * PI * (pow(DP, 2) - pow(DG, 2));
    double ii = 2 * (3 * Alpha + 1) * (pow(DP, 2) - pow(DG, 2));
    double iii = N1 * Alpha * pow(DO1, (3 * Alpha + 1) / Alpha);
    double C1 = i * pow(ii / iii, Alpha)/1000;

   // double j = pow(PI * (pow(DP, 2) - pow(DG, 2)) / 4, Alpha2 + 1);
    //double jj = ((N2 * Alpha2 * PI * pow(DO2 / 2, (3 * Alpha2 + 1) / Alpha2)) / (3 * Alpha2 + 1)) * pow(1 / 2 * G2 * L, 1 / Alpha2);
   // double jjj = (Alpha2 * PI * pow(HP, (2 * Alpha2 + 1) / Alpha2) * (DC + DP) / (4 * (2 * Alpha2 + 1))) * pow(1 / 2 * G2 * (fabs(Tstrain) - S), 1 / Alpha2);
   // double C2 = (j / pow(jj + jjj, Alpha2))/1000 + C1;
   
    double j = PI * (pow(DP, 2) - pow(DG, 2)) / 4;
    double jj = pow(j, Alpha + 1);
    double mm = N2 * Alpha * PI * pow(DO2 / 2, (3 * Alpha + 1) / Alpha);
    double mmm = mm / (3 * Alpha + 1);
    double mmmm = mmm * pow(1 / (2 * G2 * L), 1 / Alpha);
    double nn = Alpha * PI * pow(HC, (2 * Alpha) / Alpha) * (DP + DC);
    double ttt = Alpha * PI * pow(HP, (2 * Alpha) / Alpha) * (DP + DC);
    double nnn = 4 * (2 * Alpha + 1);
    double nnnn = pow(1 / (2 * G2 * (fabs(Tstrain) - S)), 1 / Alpha);
    double nnnnn = pow(1 / (2 * G2 * L), 1 / Alpha);
    double gg1 = (nn / nnn) * nnnn;
    double gg2 = (ttt / nnn) * nnnn;
    double gg3 = (ttt / nnn) * nnnnn;

    double C4 = (jj / pow(mmmm + gg3, Alpha) / 1000)+C1;
    double C2 = (jj / pow(mmmm + gg1, Alpha)/1000) + C1;
    double C3 = (jj / pow(mmmm + gg2, Alpha)/1000) + C1;
    //double jjjj = (Alpha2 * PI * pow(HP, (2 * Alpha2 + 1) / Alpha2) * (DC + DP) / (4 * (2 * Alpha2 + 1))) * pow(1 / 2 * G2 * L, 1 / Alpha2);
    //double C3 = (j / pow(jj + jjjj, Alpha2))/1000 + C1;
    //double C3 = C1 + C1;
    //return (v - sgn(fd) * pow(fabs(fd) / C, 1.0 / Alpha)) * K;
//}
  if (fabs(Tstrain) <= S){
      return (v - sgn(fd) * pow(fabs(fd) / C1, 1.0 / Alpha)) * K;
  }
  else if (fabs(Tstrain) > S && fabs(Tstrain) < S + LC) {
      return (v - sgn(fd) * pow(fabs(fd) / C2, 1.0 / Alpha)) * K;
  }
  else if (fabs(Tstrain) >= S + LC && fabs(Tstrain) < S + L) {
      return (v - sgn(fd) * pow(fabs(fd) / C3, 1.0 / Alpha)) * K;
  }
  else if (fabs(Tstrain) >= S + L) {
      return (v - sgn(fd) * pow(fabs(fd) / C4, 1.0 / Alpha)) * K;
  }
}

//K(k), G1(g1), G2(g2), Alpha1(a1), Alpha2(a2), L(l), DP(dp), DG(dg), N1(n1), N2(n2), DO1(do1), DO2(do2), DC(dc), S(s), H(h), LGap(lgap), NM(nm), RelTol(reltol), AbsTol(abstol), MaxHalf(maxhalf)

UniaxialMaterial *
APDVFD::getCopy(void)
{
    APDVFD *theCopy = new APDVFD(this->getTag(), K, G1, G2, Alpha, L, LC, DP, DG, N1, N2, DO1, DO2, DC, S, HP, HC, LGap, NM, RelTol, AbsTol, MaxHalf);
    // Converged state variables
        theCopy->Cstrain = Cstrain;
        theCopy->Cstress = Cstress;
        theCopy->Ctangent = Ctangent;
        theCopy->CVel = CVel;
		theCopy->Cpugr = Cpugr;
		theCopy->Cnugr = Cnugr; 
        
        // Trial state variables
		theCopy->Tstrain = Tstrain;
		theCopy->Tstress = Tstress; 
        theCopy->Ttangent = Ttangent;
        theCopy->TVel = TVel;
		theCopy->Tpugr = Tpugr;
		theCopy->Tnugr = Tnugr;        
    return theCopy;
}

int 
APDVFD::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(28);
  data(0) = this->getTag();

  // Material properties
  data(1) = K;
  data(2) = G1;
  data(3) = G2;
  data(4) = Alpha;
  data(5) = L;
  data(6) = LC;
  data(7) = DP;
  data(8) = DG;
  data(9) = N1;
  data(10) = N2;
  data(11) = DO1;
  data(12) = DO2;
  data(13) = DC;
  data(14) = S;
  data(15) = HP;
  data(16) = HC;
  data(17) = LGap;
  data(18) = NM;
  data(19) = RelTol;
  data(20) = AbsTol;
  data(21) = MaxHalf;
  
  // State variables from last converged state
  data(22) = Cstrain;
  data(23) = Cstress;
  data(24) = Ctangent;
  data(25) = CVel;
  data(26) = Cpugr;
  data(27) = Cnugr;

        
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "APDVFD::sendSelf() - failed to send data\n";

  return res;
}

int 
APDVFD::recvSelf(int cTag, Channel &theChannel,
                               FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(28);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "APDVFD::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
        
        // Material properties
        K = data(1);
        G1 = data(2);
        G2 = data(3);
		Alpha = data(4);
        L = data(5);
        LC = data(6);
        DP = data(7);
        DG = data(8);
        N1 = data(9);
        N2 = data(10);
        DO1 = data(11);
        DO2 = data(12);
        DC = data(13);
        S = data(14);
        HP = data(15);
        HC = data(16);
		LGap = data(17);
	    NM = data(18);
        RelTol = data(19);
		AbsTol = data(20);
		MaxHalf = data(21);
        
        // State variables from last converged state  data(6) = L;

        Cstrain = data(22);
        Cstress = data(23);
        Ctangent = data(24);
        CVel = data(25);
		Cpugr = data(26);
		Cnugr = data(27);
          
  }
    
  return res;
}

int
APDVFD::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(K);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"eta") == 0) {
    param.setValue(G1);
    return param.addObject(4, this);
  }
  else if (strcmp(argv[0],"eta") == 0) {
      param.setValue(G2);
      return param.addObject(5, this);
  }
  return -1;
}

int 
APDVFD::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    K = info.theDouble;
    return 0;
  case 4:
    G1 = info.theDouble;
    return 0;
  case 5:
      G2 = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

void 
APDVFD::Print(OPS_Stream &s, int flag)
{
    s << "APDVFD tag: " << this->getTag() << endln;
    s << "  K: " << K << endln; 
    s << "  G1: " << G1 << endln;
    s << "  G2: " << G2 << endln;
    s << "  Alpha: " << Alpha << endln;
    s << "  L: " << L << endln;
    s << "  LC: " << LC << endln;
    s << "  DP: " << DP << endln;
    s << "  DG: " << DG << endln;
    s << "  N1: " << N1 << endln;
    s << "  N2: " << N2 << endln;
    s << "  DO1: " << DO1 << endln;
    s << "  DO2: " << DO2 << endln;
    s << "  DC: " << DC << endln;
    s << "  S: " << S << endln;
    s << "  HP: " << HP << endln;
    s << "  HC: " << HC << endln;
	s << "  LGap: " << LGap << endln; 
	s << "  NM: " << NM << endln; 
    s << "  RelTol: " << RelTol << endln;
	s << "  AbsTol: " << AbsTol << endln;
    s << "  MaxHalf: " << MaxHalf << endln;
        
}
