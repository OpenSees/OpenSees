/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright for this material model is by Quan Gu  Xiamen Unv.   **  
**                                                                    **
** Commercial use of this program without express permission of the   **  
** Authors is strictly prohibited.                                    **
**                                                                    **
** Developed by:                                                      **
**   Quan Gu (quangu@xmu.edu.cn)                                      **  
**                                                                    **     
** reference: Hofstetter, G., J.C. Simo and R.L. Taylor, \A modified  **
** cap-model: Closest point solution algorithms\, Computers &         **
** Structures, 46 (1993),203-214.                                     **
**                                                                    **
** ****************************************************************** */

// stress and strain all defined as: elongation is negative; compression is positive (+)


//#include <math.h>
#include <CapPlasticity.h>
#include <Information.h>

#include <Channel.h>  
#include <FEM_ObjectBroker.h>

#include <T2Vector.h>
#include <MaterialResponse.h>

#include <elementAPI.h>

#include <fstream>            // Quan Gu   2013 March   HK
using std::ofstream;          // Quan Gu   2013 March   HK
using std::ios;               // Quan Gu   2013 March   HK
  

Vector CapPlasticity::tempVector(6);
Matrix CapPlasticity::tempMatrix(6,6);  

void * OPS_ADD_RUNTIME_VPV(OPS_CapPlasticity) {
  int tag;
  int ndm =3;
  double rho = 0.0;
  double G = 1.0e10;
  double K = 1.1e10;
  double X = 1.1032e8;
  double D = 4.6412e-10;
  double W = 0.42;
  double R = 4.43;
  double lambda = 7.9979e6;
  double theta = 0.11;
  double beta = 6.3816e-8;
  double alpha = 2.6614e7;
  double T = -2.0684e6; 
  double tol = 1.0e-10;

  int numArgs = OPS_GetNumRemainingInputArgs();

  int iData[2];
  double dData[10];

  int numData = 2;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial CapPlasticisty \n";
    return 0;
  }  
  tag = iData[0]; ndm = iData[1];


  numData = 3;  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid double values: nDMaterial CapPlasticity " << tag << endln;
    return 0;
  }  
  rho = dData[0]; G = dData[1]; K = dData[2];

  if (numArgs == 10) {
    numData = 10;      
    if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid double values: nDMaterial CapPlasticity " << tag << endln;
      return 0;
    }  
    X = dData[0]; 
    D = dData[1];
    W = dData[2];
    R = dData[3];
    lambda = dData[4];
    theta = dData[5];
    beta = dData[6];
    alpha = dData[7];
    T = dData[8];
    tol = dData[9];
  } //end if
  
  NDMaterial *theMaterial = new CapPlasticity(  tag,
						G,
						K,
						rho,
						X,
						D,
						W,
						R,
						lambda,
						theta,
						beta,
						alpha,
						T, 
						ndm,
						tol) ;
  
  return theMaterial;
}



CapPlasticity::CapPlasticity( int    pTag,
			      double pG,
			      double pK,
			      double pRho,
			      double pX,
			      double pD,
			      double pW,
			      double pR,
			      double pLambda,
			      double pTheta,
			      double pBeta,
			      double pAlpha,
			      double pT,
			      int pNdm, 
			      double pTol_k
			      )
: NDMaterial(pTag,ND_TAG_CapPlasticity),
  CStrain(6),
  CPlastStrain(6),
  CStress(6),
  strain(6), 
  plastStrain(6), 
  stress(6),
  stressDev(6),
  theTangent(6,6)
{
  shearModulus = pG;
  bulkModulus = pK;
  rho = pRho;
  X = pX;
  D = pD;
  W = pW;
  R = pR;
  lambda = pLambda;
  theta = pTheta;
  beta = pBeta;
  alpha = pAlpha;
  T = pT;
  // -------------2012,1,15 --------
  if (T>0) T = -T;
  
  ndm = pNdm;
  tol_k = pTol_k; 
  stressI1 =0.0;
  
  flag =1;
  // --
  revertToStart();
  
  debug =0;
  SHVs =0;
  parameterID =0;
  
  // -- theMode 
  theMode =-10;
  
};

CapPlasticity::CapPlasticity( const CapPlasticity & a)
  : NDMaterial(a.getTag(),ND_TAG_CapPlasticity),CStrain(6),CPlastStrain(6),CStress(6),
    strain(6),plastStrain(6), stress(6), stressDev(6), theTangent(6,6)
{
  shearModulus = a.shearModulus;
  bulkModulus = a.bulkModulus;
  rho = a.rho;
  X = a.X;
  D = a.D;
  W = a.W;
  R = a.R;
  lambda = a.lambda;
  theta = a.theta;
  beta = a.beta;
  alpha = a.alpha;
  T = a.T;
  // -------------2012,1,15 --------
  if (T>0) T = -T;
  
  ndm = a.ndm;
  tol_k =a.tol_k; 
  stressI1 =0.0;
  
  flag =1;
  // --
  revertToStart();
  
  SHVs =0;
  parameterID =0;
	
};

CapPlasticity::~CapPlasticity( ) {
  return;
  
};

double CapPlasticity::getRho(void) {
  return rho;
};

int CapPlasticity::setTrialStrain(const Vector &pStrain) {  // strain from element is eng. strain!
  
  //  static Vector temp(6);
  if (ndm==3 && pStrain.Size()==6) 
    strain = -1.0*pStrain;
  else if (ndm==2 && pStrain.Size()==3) {
    strain[0] = -1.0*pStrain[0];
    strain[1] = -1.0*pStrain[1];
    strain[2] = 0.0;
    strain[3] = -1.0*pStrain[2];
    strain[4] = 0.0;
    strain[5] = 0.0;
  }
  else {
    opserr << "Fatal:CapPlasticity:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << pStrain.Size() << endln;
    //exit(-1);
	opserr<< "Warning: errors in CapPlasticity::setTrialStrain"<<endln;
  }
  
  // ----- change to real strain instead of eng. strain
  
  for ( int i = 3; i<6; i++) 
    strain[i] /=2.0; 
  
  
  
  return 0;
  
};

int CapPlasticity::setTrialStrain(const Vector &pStrain, const Vector &r) {
  
  return setTrialStrain(pStrain);
  
};

int CapPlasticity::setTrialStrainIncr(const Vector &pStrainRate) {
  
  // ----- change to real strain instead of eng. strain
  // ---- since all strain in material is the true strain, not eng.strain. 
  
  for (int i=0; i<3;i++) {
    tempVector(i) =  pStrainRate(i);
    tempVector(i+3) =  pStrainRate(i+3)/2.0;
  }
  
  if (ndm==3 && pStrainRate.Size()==6) 
    strain = CStrain-pStrainRate;
  
  else if (ndm==2 && pStrainRate.Size()==3) {
    strain[0] = CStrain[0]-pStrainRate[0];
    strain[1] = CStrain[1]-pStrainRate[1];
    strain[2] = 0.0;
    strain[3] = CStrain[2]-pStrainRate[2];
    strain[4] = 0.0;
    strain[5] = 0.0;
  }
  else {
    opserr << "Fatal:CapPlasticity:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << pStrainRate.Size() << endln;
	opserr<< "Warning: errors in CapPlasticity::setTrialStrainIncr"<<endln;	
   // exit(-1);
  }
  
  
  
  
  return 0;
  
  
};

int CapPlasticity::setTrialStrainIncr(const Vector &pStrainRate, const Vector &r) {
  return setTrialStrainIncr(pStrainRate);
};



const Vector & CapPlasticity::getStrain(void) {
  
  if (ndm==3){
    tempVector = -1.0*strain;
    return tempVector;
    
  }
  
  else {
    static Vector workV(3);//, temp6(6);
    workV[0] = -1.0*strain[0];
    workV[1] = -1.0*strain[1];
    workV[2] = -1.0*strain[3];
    return workV;
  }
  
  
  
};



int CapPlasticity::commitState(void)  {
  CStrain = strain;
  CStress = stress;
  CPlastStrain = plastStrain;
  CHardening_k =hardening_k;
  
  
  //	debug =1;
  
  return 0;
  
};



int CapPlasticity::revertToLastCommit(void)  {
  return 0;
};

int CapPlasticity::revertToStart(void)  {
  
  CStrain.Zero();
  CPlastStrain.Zero();
  CStress.Zero();
  strain.Zero();
  plastStrain.Zero();
  stress.Zero();
  CHardening_k = Newton_k(tol_k, 0);
  return 0;
};

NDMaterial * CapPlasticity::getCopy(void)  {
  
  CapPlasticity * copy = new CapPlasticity(*this);
  return copy;
  
};

NDMaterial * CapPlasticity::getCopy(const char *code)  {
  if (strcmp(code,this->getType()) == 0) {
    CapPlasticity * copy = new CapPlasticity(*this);
    return copy;
  }
  else
    return 0;
};

const char * CapPlasticity::getType(void) const  {
  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
};

int CapPlasticity::getOrder(void) const  {
  return (ndm == 2) ? 3 : 6;
};

int CapPlasticity::sendSelf(int commitTag, Channel &theChannel)  {return 0;};

int CapPlasticity::recvSelf(int commitTag, Channel &theChannel,  
			    FEM_ObjectBroker &theBroker )  {return 0;};

  
Response*
CapPlasticity::setResponse (const char **argv, int argc, OPS_Stream &output) {

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
    return new MaterialResponse(this, 1, stress);
  
  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
    return new MaterialResponse(this, 2, strain);
  
  else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"Tangent") == 0)
    return new MaterialResponse(this, 3, theTangent);
  
  
  else if (strcmp(argv[0],"plasticStrain") == 0 || strcmp(argv[0],"plasticStrains") == 0)
    return new MaterialResponse(this, 4, plastStrain);
  
  else if (strcmp(argv[0],"k") == 0)
    return new MaterialResponse(this, 5, hardening_k);
  
  else if (strcmp(argv[0],"stress_and_k") == 0 ){
    static Vector dummy(7); 
    return new MaterialResponse(this, 6, dummy);
  }
  
  else
    return NDMaterial::setResponse(argv, argc, output);
  
  
};

int CapPlasticity::getResponse (int responseID, Information &matInfo)  {
  
  switch (responseID) {
  case -1:
    return -1;
  case 1:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) =stress;
    return 0;
    
  case 2:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = strain;
    return 0;
    
  case 3:
    if (matInfo.theMatrix != 0)
      *(matInfo.theMatrix) = theTangent;
    return 0;
    
  case 4:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = plastStrain;
    return 0;
    
  case 5:
    matInfo.setDouble(this->hardening_k);
    return 0;
    
  case 6:
    static Vector dummy(7); 
    for (int i=0; i<6; i++)
      dummy(i) = stress(i);
    dummy(6) = this->hardening_k;
    *(matInfo.theVector) = dummy;
    return 0;
  }
  
  return NDMaterial::getResponse(responseID, matInfo);
  
};

void CapPlasticity::Print(OPS_Stream &s, int flag)  {return;};

// --------------

double CapPlasticity::CapBoundL(double k){
  if (k>0) 
    return k; 
  return 0;
};

// --------------
double CapPlasticity::CapBoundX(double k){  // function of X
  return k+R*failureEnvelop(k);
}; 

// --------------
double CapPlasticity::failureEnvelop(double I){
  
  return alpha - lambda*exp(-beta*I)+theta*I; 
};

// --------------
double CapPlasticity::CapSurface(double normS,double I, double k){  // Fc
  
  double res = 0;
  
  //	if (I< CapBoundX(k))
  res = pow(normS*normS+(I-CapBoundL(k))*(I-CapBoundL(k))/R/R,0.5);
  return res;
};

// --------------
double CapPlasticity::failureEnvelopDeriv(double I){
  return lambda*beta*exp(-beta*I)+theta;
};

// --------------
double CapPlasticity::hardeningParameter_H(double k2, double k1){

  return W*(exp(-D*CapBoundX(k1))-exp(-D*CapBoundX(k2)));
}; 

int CapPlasticity::findMode(double normS, double I1, double k){
  
  
  int mode =-1;
  
  if ((I1 <= T) && (normS <= failureEnvelop(T)))
    mode =1;
  else if ((I1 <= T) && (normS >= failureEnvelop(T)) && (normS <= failureEnvelop(T)+2.0*shearModulus/(9.0*bulkModulus)*(T-I1)/failureEnvelopDeriv(T)))
    mode =2;
  else if ((I1>=k) && (failureEnvelop(k)<=CapSurface( normS, I1, k))) {
    
    mode =3;
    
  }
  else if ((I1<=k) && (normS >= failureEnvelop(k) + 2.0*shearModulus/(9.0*bulkModulus)*(k-I1)/failureEnvelopDeriv(k)))
    mode = 4;
  else if (((normS <= failureEnvelop(I1)) && (I1 >= T) && (I1 < k )) ||(( I1>= k) && (failureEnvelop(k)>=CapSurface( normS, I1, k))))
    mode = 6;
  else if ((normS <= failureEnvelop(k) + 2.0*shearModulus/(9.0*bulkModulus)*(k-I1)/failureEnvelopDeriv(k)) && (normS >= failureEnvelop(T)+2.0*shearModulus/(9.0*bulkModulus)*(T-I1)/failureEnvelopDeriv(T))){
    
    mode = 5;
  }
  return mode;
};

// --------------
double CapPlasticity::Newton_k(double tol, int mode /*, double normS, double I1, double k */){
  
  
  double solution;
  int maxIter =200;
  
  if ( mode == 0){
    // very initial step, get kn for the next step 
    double k = 0; 
    double f = CapBoundX(k) - X;
    int i =1;
    while (( i<= maxIter) && (fabs(f) > tol)){
      double dfdk = 1+ R*failureEnvelopDeriv(k);
      k=k-f/dfdk;
      f = CapBoundX(k) - X;
      i++;
    }//while
    
    if (fabs(f) > tol) {
      opserr<< "Fatal : Newton algorithm does not converge, in CapPlasticity, mode =0! \n";
     // exit(-1);
    }
    solution = k;
    
  }  // mode ==0
  
  else if ( mode ==1 || mode ==2 || mode ==5 ) {
    
    double k = CHardening_k; 
    double f = deltPlastStrainI1 - hardeningParameter_H( k, CHardening_k);
    int i =1;
    while (( i<= maxIter) && (fabs(f) > tol)){
      double dfdk = -W*D*(1+R*failureEnvelopDeriv(k))*exp(-D*CapBoundX(k));
      k=k-f/dfdk;
      f = deltPlastStrainI1 - hardeningParameter_H( k, CHardening_k);
      i++;
    }//while
    
    if (fabs(f) > tol) {
      opserr<< " Newton algorithm does not converge, in CapPlasticity, mode = "<<" "<<"mode"<<endln;
      //opserr<< "Warning: errors in CapPlasticity::setTrialStrainIncr"<<endln;	
    }
    solution = k;
    
  }  // mode = 1,2,5-2
  
 if (solution <0 ) {
   opserr<<"Warning: CapPlasticity:: Newton_k, solution <0! mode is " << mode 
	 << "! k should be adjusted to CHardening_k! " << endln;    // --- April 2013. 
   solution = CHardening_k;         // change to stress I1 according to the definition of I1_p, 2013-6-14 
 
 }
	
 return solution;

}; 

double CapPlasticity::Newton_I1(double tol, int mode, double normS, double I1_trial){
  
  double solution;
  int maxIter =200;
  double relative_tol = tol*fabs(I1_trial);
  if (relative_tol<tol)  relative_tol =tol; 
  
  
  double relative_tol_2 = relative_tol;
  if ( relative_tol_2 >1.0e-7) 
    relative_tol_2 = 1.0e-7;
  
  
  if ( mode == 5){   // get I1, that is the first step of mode 5
    double I1 = I1_trial; 
    int i =1;
    double deltGammar = (normS-failureEnvelop(I1))/(2.0 * shearModulus);
    double f = 9.0 * bulkModulus * failureEnvelopDeriv(I1) * deltGammar + I1_trial-I1;
    
    
    while (( i<= maxIter) && (fabs(f) >= tol)){
      double dfdI1 = 9.0*bulkModulus*(-deltGammar*lambda*beta*beta*exp(-beta*I1)-1.0/(2.0*shearModulus)*failureEnvelopDeriv(I1)*failureEnvelopDeriv(I1))-1.0;
      I1=I1-f/dfdI1;
      i++;
      deltGammar = (normS-failureEnvelop(I1))/(2.0 * shearModulus);
      f = 9.0 * bulkModulus * failureEnvelopDeriv(I1) * deltGammar + I1_trial-I1;
    }//while
    
    
    if (fabs(f) > relative_tol) {
      opserr<< "mode =5. Newton algorithm does not converge, in CapPlasticity, Newton_I1 mode =5! ";
      //exit(-1);
    }
    solution = I1;
    
  }  // mode ==5
  
  else if ( mode == 3) {
    
    
    int i=1;
    double k = CHardening_k; 
    double deltGammar =0;
    double f =0;
    double I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(k, CHardening_k);
    if (k>I1 + relative_tol_2)  
      flag =0;
    else { // flag =1
      if (fabs( k - I1)<relative_tol_2)
	deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
      else
	deltGammar = R*R*hardeningParameter_H(k, CHardening_k)*failureEnvelop(k)/(3.0*(I1-k));
      
      double tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(k));
      double tmp2 = (I1_trial - k)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(k)));
      f = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(k);
    } // else flag =1
    
    double dHdk;
    double dIdk;
    double temp1;
    double ddeltGammardk;
    double temp2;
    double temp3;
    double temp4;
    double temp5;
    double dfdk;
    
    ///////////////////////////////////////
    
    double relative_tol = tol* failureEnvelop(CHardening_k); 
    
    
    while (( i<= maxIter) && (flag==1) && (fabs(f) > tol)){
      dHdk=W*D*(1+R*failureEnvelopDeriv(k))*exp(-D*CapBoundX(k));
      dIdk=-3.0*bulkModulus*dHdk;
      temp1=(dHdk*failureEnvelop(k)+failureEnvelopDeriv(k)*hardeningParameter_H(k, CHardening_k))*(I1-k);
      
      if (fabs(k - I1) < relative_tol_2)
	ddeltGammardk=0.0;
      else
	ddeltGammardk=R*R*(temp1-(dIdk-1)*hardeningParameter_H(k, CHardening_k)*failureEnvelop(k))/(3*(I1-k)*(I1-k));
      
      temp2=1.0+2.0*shearModulus*deltGammar/failureEnvelop(k);
      temp3=R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(k));
      temp4=2.0*shearModulus*(ddeltGammardk*failureEnvelop(k)-deltGammar*failureEnvelopDeriv(k))/(failureEnvelop(k)*failureEnvelop(k));
      temp5=-R-9.0*bulkModulus*deltGammar/(R*failureEnvelop(k))-9.0*bulkModulus*(I1_trial-k)*(ddeltGammardk*failureEnvelop(k)-deltGammar*failureEnvelopDeriv(k))/R/failureEnvelop(k)/failureEnvelop(k);
      dfdk=(-1.0*normS*normS*temp4/(temp2*temp2*temp2)+(I1_trial-k)*temp5/(temp3*temp3*temp3))/ pow(normS*normS/(temp2*temp2)+(I1_trial-k)*(I1_trial-k)/(temp3*temp3),0.5)-failureEnvelopDeriv(k);
      k = k-f/dfdk;
      i++;
      I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(k, CHardening_k);
      if (k>I1 + relative_tol_2)  
	flag =0;
      else { // flag =1
	if (fabs(k - I1) < relative_tol_2)
	  deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
	else
	  deltGammar = R*R*hardeningParameter_H(k, CHardening_k)*failureEnvelop(k)/(3.0*(I1-k));
	
	double tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(k));
	double tmp2 = (I1_trial - k)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(k)));
	f = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(k);
	
      } // else flag =1
      
      
    } // end while 	 
   	 if (k<0) {
		 opserr<<"Warning:  Newton_I1: mode =3. get k<0; adjusted to CHardening_k!!"<<endln; 
         solution = CHardening_k;
	 } 
	 
    solution = k;
  }  // if mode ==3
  
  return solution;
  
}; 

double CapPlasticity::Bisection(double tol, double normS, double I1_trial ){
  
  double x_low = CHardening_k;
  int maxIter = 200; 
  
  
  // -- to get the lower bound using I1_n+1 = k_n+1 , using newton algorithm
  
  int i=1;
  double dHdk =0;
  double dfdk =0;
  
  double k = CHardening_k;
  double f = I1_trial -3.0*bulkModulus*hardeningParameter_H(k, CHardening_k)-k;
  double relative_tol = tol*k;
  
  while ( i <= maxIter && fabs(f) > tol) {
    dHdk = W*D*(1+R*failureEnvelopDeriv(k))*exp(-D*CapBoundX(k));
    dfdk = -3.0*bulkModulus*dHdk-1.0;
    k = k-f/dfdk;
    i++;
    f = I1_trial -3.0*bulkModulus*hardeningParameter_H(k, CHardening_k)-k;
    
  } //while
  
  if (fabs(f) > relative_tol) {
    
    opserr<< "Warning: Newton can not converge in CapPlasticity::Bisection"<<endln;
    //exit(-1);
  }
  double x_up = k;
  
  // two bounds need have opposite sign
  
  // --- x_low 
  
  double I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(x_low, CHardening_k);
  double deltGammar =0;
  if (x_low == I1) {
    double debug1 = failureEnvelop(CHardening_k); 
    deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
  }
  else {
    double debug1 = hardeningParameter_H(x_low, CHardening_k);
    double debug2 = failureEnvelop(x_low);
    
    deltGammar = R*R*hardeningParameter_H(x_low, CHardening_k)*failureEnvelop(x_low)/(3.0*(I1-x_low));
  }
  
  double tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(x_low));
  double tmp2 = (I1_trial - x_low)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(x_low)));
  double f_low = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(x_low);
  
  // ---- x_up
  
  I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(x_up, CHardening_k);
  if (x_low == I1)
    deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
  else
    deltGammar = R*R*hardeningParameter_H(x_up, CHardening_k)*failureEnvelop(x_up)/(3.0*(I1-x_up));
  
  tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(x_low));
  tmp2 = (I1_trial - x_low)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(x_low)));
  double f_up = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(x_up);
  
	double incr = (x_up-x_low)*0.05;
	while ( (f_low*f_up>0) && (x_low < x_up)) {
		
		x_up = x_up -incr;
    
    // --- x_low 
    
    I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(x_low, CHardening_k);
    if (x_low == I1)
      deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
    else
      deltGammar = R*R*hardeningParameter_H(x_low, CHardening_k)*failureEnvelop(x_low)/(3.0*(I1-x_low));
    
    tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(x_low));
    tmp2 = (I1_trial - x_low)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(x_low)));
    f_low = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(x_low);
    
  }// while
  
  if (f_low*f_up>0 ) {
    opserr <<"Warning2: Bisection can not converge in  CapPlasticity::Bisection! "<<endln;
    //exit(-1);
    
  }
  
  // ----- bisection method --
  
  k = (x_low+x_up)/2.0;
  
  I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(k, CHardening_k);
  if (k == I1)
    deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
  else
    deltGammar = R*R*hardeningParameter_H(k, CHardening_k)*failureEnvelop(k)/(3.0*(I1-k));
  
  tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(k));
  tmp2 = (I1_trial - k)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(k)));
  f = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(k);
  
  maxIter = 500;
  i =0;
  
  relative_tol = tol*failureEnvelop(k);
  
  while (fabs(f)>tol && i<maxIter) {
    if (f*f_low<0) {
      x_up = k;
      f_up = f;
    }
    else {
      x_low = k;
      f_low =f;
    }
    
    k = (x_low+x_up)/2.0;
    
    I1 = I1_trial-3.0*bulkModulus*hardeningParameter_H(k, CHardening_k);
    if (k == I1)
      deltGammar = (normS-failureEnvelop(CHardening_k))/(2.0 * shearModulus);
    else
      deltGammar = R*R*hardeningParameter_H(k, CHardening_k)*failureEnvelop(k)/(3.0*(I1-k));
    
    tmp1 = normS/(1+2.0*shearModulus*deltGammar/failureEnvelop(k));
    tmp2 = (I1_trial - k)/(R+9.0*bulkModulus*deltGammar/(R*failureEnvelop(k)));
    f = pow(tmp1*tmp1+tmp2*tmp2,0.5)-failureEnvelop(k);
    
    i++;
  } // while
  
  if (fabs(f) > relative_tol) {
    opserr << "Warning3:No convergence in CapPlasticity::Bisection\n";
    //exit(-1); 
  }
  
  flag =1; 
 if (k <0 ) {
	opserr << "Fatal: CapPlasticity:: Bisection, k <0! mode is 3 ! k is adjusted to CHardening_k !!!!!\n";    // --- April 2013. 
	k = CHardening_k;

 }
  return k; 
  
}; // 



 


const Matrix & CapPlasticity::getTangent(void) {
	
  
  if (ndm==3) 
    return theTangent;
  else {
    static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = theTangent(0,3);
    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = theTangent(1,3);
    workM(2,0) = theTangent(3,0);
    workM(2,1) = theTangent(3,1);
    workM(2,2) = theTangent(3,3);
    return workM;
  }


};

const Matrix & CapPlasticity::getInitialTangent(void) {
  return getTangent();
};

const Vector & CapPlasticity::getStress(void) {
  
  
  double CPlastStrainI1 = CPlastStrain(0)+CPlastStrain(1)+CPlastStrain(2);
  
  Vector CPlastStrainDev = CPlastStrain;
  Vector unitVector2(6);
  for (int   i=0; i<3; i++) {
    unitVector2(i)=1.0;
    unitVector2(i+3)=0.0;
  }
  CPlastStrainDev.addVector(1.0, unitVector2, -CPlastStrainI1/3.0);
  
  double strainI1 = strain(0)+strain(1)+strain(2); 
  
  Vector strainDev = strain;
  Vector deltPlastStrainDev(6); 
  
  strainDev.addVector(1.0, unitVector2, -1.0*strainI1/3.0);
  
  Vector TStressDev = 2.0 * shearModulus * (strainDev - CPlastStrainDev);
  
  double TStressI1 = 3.0*bulkModulus * (strainI1 - CPlastStrainI1);
  
  
  
  // --- iteration to get k0 ----------
  
  //	double k = Newton_k(tol_k, 0);
  
  double normTS = pow( TStressDev && TStressDev, 0.5);  // cannot use norm() 
  double normS =0.0;
  
  int mode = findMode(normTS, TStressI1, CHardening_k);  // used everywhere in this code!!
  
  double deltGammar1 =0.0;
  double deltGammar2 =0.0;
  double deltGammar3 =0.0;
  
  
  if (debug ==1) {
    opserr<<"mode is "<<mode <<endln;
    opserr<<"strain is:" <<strain; 
    opserr<<"stress is:" <<stress<<endln;
  }
  
  if (mode ==1) {
    deltGammar3 = (T-TStressI1)/(9.0*bulkModulus);
    this->stressI1 = T;
    stressDev = TStressDev;
    deltPlastStrainDev.Zero();
    deltPlastStrainI1 = -3.0*deltGammar3;
    
    hardening_k = Newton_k(tol_k, mode); //
    
  } // mode ==1
  
  else if (mode ==2){
    deltGammar1 = (normTS-failureEnvelop(T))/(2.0*shearModulus);
    deltGammar3 = (T-TStressI1)/(9.0*bulkModulus) - deltGammar1*failureEnvelopDeriv(T); 
    this->stressI1 = T;
    normS =failureEnvelop(T);
    //normTS = pow(TStressDev && TStressDev, 0.5);
    
    stressDev.addVector(0.0, TStressDev, normS/ normTS);
    deltPlastStrainDev.addVector( 0.0, TStressDev, deltGammar1/normTS);
    deltPlastStrainI1 = -3.0 * (deltGammar1 * failureEnvelopDeriv(T)+deltGammar3);
    hardening_k = Newton_k(tol_k, mode); //
    
    
  } //mode ==2
  
  else if (mode ==3) {
    
    deltGammar2 =0;
    //double normTS = pow( TStressDev && TStressDev ,0.5);
    
    double relative_tol = tol_k * CHardening_k;
    if ( relative_tol > 1.0e-7) relative_tol = 1.0e-7;
    
    if (fabs(TStressI1 - CHardening_k) < relative_tol){
      deltGammar2 = (normTS-failureEnvelop(CHardening_k))/(2.0*shearModulus);
      hardening_k = TStressI1;
    }
    else {
      //double normTS = pow(TStressDev && TStressDev, 0.5);
      hardening_k = Newton_I1(tol_k, mode, normTS, TStressI1);
      if (flag ==0 )
	hardening_k = Bisection(tol_k, normTS, TStressI1); 
    } //if (TStressI1 == CHardening_k)
    
    this->stressI1 = TStressI1-3.0*bulkModulus*hardeningParameter_H(hardening_k, CHardening_k);
    deltGammar2 =R*R*hardeningParameter_H(hardening_k, CHardening_k)*failureEnvelop(hardening_k)/(3.0*(stressI1-hardening_k));
    normS=normTS/(1+2.0*shearModulus*deltGammar2/failureEnvelop(hardening_k));
    stressDev.addVector(0.0, TStressDev, normS/ normTS);
    this->stressI1=hardening_k+(TStressI1-hardening_k)/(1.0+9.0*bulkModulus*deltGammar2/(R*R*failureEnvelop(hardening_k)));
    deltPlastStrainDev.addVector(0.0, stressDev, deltGammar2/CapSurface(normS, stressI1, hardening_k));	
    deltPlastStrainI1 = 3.0*deltGammar2*(stressI1-hardening_k)/(R*R*CapSurface(normS, stressI1, hardening_k));
    
  } // mode ==3
  
  else if (mode ==4){
    deltGammar1 = (CHardening_k - TStressI1)/(9.0*bulkModulus*failureEnvelopDeriv(CHardening_k));
    deltGammar2 = (normTS-failureEnvelop(CHardening_k))/(2.0*shearModulus)-deltGammar1;
    this->stressI1 = CHardening_k;
    normS = failureEnvelop(CHardening_k);
    stressDev.addVector(0.0, TStressDev, normS/ normTS);
    deltPlastStrainDev.addVector(0.0, stressDev, (deltGammar1+deltGammar2)/normS);
    deltPlastStrainI1 = -3.0*deltGammar1 * failureEnvelopDeriv(CHardening_k);
    hardening_k = CHardening_k;
    
  } //mode ==4
  
  else if (mode ==6) {
    normS = normTS;
    stressDev = TStressDev;
    stressI1 = TStressI1;
    deltPlastStrainDev.Zero();
    deltPlastStrainI1 = 0.0;
    hardening_k = CHardening_k;
    
    
  } // mode ==6
  
  else if (mode ==5){
    
    
    double tol_I1 = tol_k;
    stressI1 = Newton_I1( tol_I1, mode, normTS, TStressI1);
    
    deltGammar1 = (normTS-failureEnvelop(stressI1))/(2.0*shearModulus);
    normS = normTS - 2.0*shearModulus * deltGammar1;
    stressDev.addVector(0.0, TStressDev, normS/ normTS);
    deltPlastStrainDev.addVector(0.0, stressDev, deltGammar1/normS);
    deltPlastStrainI1 = -3.0*deltGammar1*failureEnvelopDeriv(stressI1);
    
    hardening_k = Newton_k(tol_k, mode); //
    
    
  } //mode ==5
  
  // --- compute strain and stress --
  
  double plastStrainI1 = CPlastStrainI1 + deltPlastStrainI1;
  Vector plastStrainDev = CPlastStrainDev + deltPlastStrainDev;
  
  plastStrain.addVector(0.0,  plastStrainDev, 1.0);
  plastStrain.addVector(1.0, unitVector2, plastStrainI1/3.0);
  
  stress.addVector(0.0, stressDev,1.0);
  stress.addVector(1.0, unitVector2, stressI1/3.0);
  
  
  computeConsistentTangent(deltGammar1, deltGammar2, deltGammar3, mode); 
  
  
  theMode = mode;  //  ÎªÁË±£´æmodeµ½³ÉÔ±º¯Êý£¬ÎªÁËgetTangent()µ÷ÊÔÓÃ¡£
  
  
  if (ndm==3){
    tempVector.addVector(0.0, stress,-1.0);
    return tempVector;
  }
  
  else {
    static Vector workV(3);//, temp6(6);
    workV[0] = -1.0*stress[0];
    workV[1] = -1.0*stress[1];
    workV[2] = -1.0*stress[3];
    return workV;
  }
  
  
};

// ---------------- for consistent tangent modulus ----
// ---------------- for consistent tangent modulus ----
// ---------------- for consistent tangent modulus ----
// ---------------- for consistent tangent modulus ----


Matrix & CapPlasticity::dF2dSigma ( int mode){ // the returned matrix has been consistent with notation (i.e., last 3 columns multiplied by 2) 

  tempMatrix.Zero();
  
  Matrix I_dev(6,6);
  I_dev.Zero();
  
  for (int i=0; i<6; i++)  	
    I_dev(i,i) =1.0;
  
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      I_dev(i,j) -= 1.0/3.0;
  
  Vector I2(6);   // unit vector order 2
  I2.Zero();
  
  for (int i=0; i<3; i++)  
    I2(i) = 1.0;
  
  
  
  if (mode ==5) {  // f1
    
    Vector stressDev = stress;
    double I1 = stress(0)+stress(1)+stress(2);
    
    for (int i=0; i<3; i++)
      stressDev(i) -= 1.0/3*I1;
    double normS = pow( stressDev && stressDev, 0.5); // stressDev.Norm(); can not use norm()
    
    Vector N = stressDev;
    N /= normS;
    
    tempMatrix.addMatrix(0.0, I_dev, 1.0/normS);
    
    double tmp = lambda*beta*beta*exp(-beta*I1); 
    
    for (int i=0; i<6; i++) {
      for ( int j=0; j<3; j++)
	tempMatrix(i,j) += -1.0/normS*N(i)*N(j)+tmp*I2(i)*I2(j);
      for (int j=3; j<6; j++)
	tempMatrix(i,j) += (-1.0/normS*N(i)*N(j)+tmp*I2(i)*I2(j))*2.0;  // To be consistent with the transformation between 4th order tensor and matrix
    }
  }
  
  else if (mode ==3) {  // f2
    
    Vector tmp (6);
    tmp = dFdSigma(3);
    
    double Fe = failureEnvelop(hardening_k);
    
    tempMatrix.addMatrix(0.0, I_dev, 1.0/Fe);
    for (int i=0; i<6; i++){
      for ( int j=0; j<3; j++)
	tempMatrix(i,j) += -1.0/Fe*tmp(i)*tmp(j)+1.0/Fe/R/R*I2(i)*I2(j);
      for (int j=3; j<6; j++)
	tempMatrix(i,j) += -1.0/Fe*tmp(i)*tmp(j)*2.0+1.0/Fe/R/R*I2(i)*I2(j)*2.0; // To be consistent with the transformation between 4th order tensor and matrix
    }
  }
  else if (mode ==1) { }
  
  else {
    opserr<< "warning: CapPlasticity::dF2dSigma() should not be called! mode is "<< mode<<endln;
    
    
  }

/*// ---- debug only ----------
   opserr<<" CapPlasticity::dF2dSigma(): " <<endln;
   for (  i =0; i<6; i++) {
		for (int j=0; j<6; j++)
			opserr<<tempMatrix(i,j)<<"   ";
		opserr<<endln;
   }
//*/

	return tempMatrix;
};

Vector & CapPlasticity::dFdSigma (int mode){
			
  Vector stressDev = stress;
  double I1 = stress(0)+stress(1)+stress(2);
  
  for (int i=0; i<3; i++)
    stressDev(i) -= 1.0/3*I1;
  
  double normS = pow( stressDev && stressDev, 0.5); //stressDev.Norm();
  
  Vector I2(6);   // unit vector order 2
  I2.Zero();
  
  for (int i=0; i<3; i++)  
    I2(i) = 1.0;
  
  
  if( mode ==5) {
    tempVector.addVector(0.0, stressDev, 1.0/normS);
    tempVector.addVector(1.0, I2, -1.0*failureEnvelopDeriv(I1));
  }
  else if (mode ==3) {
    double Fc = CapSurface(normS, I1, hardening_k); 
    tempVector.addVector(0.0, stressDev, 1.0/Fc);
    tempVector.addVector(1.0, I2, (I1-CapBoundL(hardening_k))/Fc/R/R);
  }
  else if (mode == 1){
    tempVector.addVector(0.0, I2, -1.0);
  }
  else {
    opserr<< "warning: CapPlasticity::dFdSigma() should not be called! mode is "<< mode<<endln;
    
  }
  return tempVector;
};

double CapPlasticity::dFdk (int OrderOfDerivative){
  
  double result =0.0;
  
  Vector stressDev = stress;
  double I1 = stress(0)+stress(1)+stress(2);
  
  for (int i=0; i<3; i++)
    stressDev(i) -= 1.0/3*I1;
  
  double normS = pow( stressDev && stressDev, 0.5); //stressDev.Norm();
  
  if (OrderOfDerivative ==1){
    result = -(I1 - CapBoundL(hardening_k))/R/R/CapSurface(normS, I1, hardening_k); 
    result -= failureEnvelopDeriv(hardening_k);
    
  }
  else if (OrderOfDerivative ==2){
    
    result = normS*normS/R/R/failureEnvelop(hardening_k);
    result += lambda*beta*beta*exp(-beta*hardening_k);
    
  }
  
  return result;
};

Vector CapPlasticity::dF2dSigmadk ( void){

  Vector I2(6);   // unit vector order 2
  I2.Zero();
  
  for (int i=0; i<3; i++)  
    I2(i) = 1.0;
  
  
  Vector stressDev = stress;
  double I1 = stress(0)+stress(1)+stress(2);
  
  for (int i=0; i<3; i++)
    stressDev(i) -= 1.0/3*I1;
  
  double normS = pow( stressDev && stressDev, 0.5); //stressDev.Norm();
  
  double denominator = R*R*pow(failureEnvelop(hardening_k),3) ;
  
  tempVector.addVector(0.0, stressDev, I1-hardening_k);
  tempVector.addVector(1.0, I2, -1.0*normS*normS);
  
  tempVector  /= denominator;
  
  
  return tempVector;
};

double CapPlasticity::dHdk(double hardening_k){

   double result = 0.0;
   
   result =  D * W *exp(-D*CapBoundX(hardening_k))*(1.0+R* failureEnvelopDeriv(hardening_k) );

	return result;
};


double 
CapPlasticity::tripleTensorProduct (Vector &A, Matrix &B, Vector &C){ // result = A:B:C

	if ((A.Size() !=6) ||(C.Size() !=6) || (B.noCols() !=6) || (B.noRows() !=6)) {
		opserr<<"Fatal: CapPlasticity::tripleTensorProduce() size does not match! "<<endln;
		exit(-1);
	}

	static Vector tmp(6);
	double result = 0.0;

	tmp.addMatrixVector(0.0, B, C, 1.0);
	
	for (int i=0; i<3; i++)
		result += A(i)*tmp(i);
	for (int i=3; i<6; i++)
		result += A(i)*tmp(i)*2.0;

	return result;

}; 



double 
CapPlasticity::dFdIdk(void){

	double result =0;

	if (hardening_k >=0) {

	
		Vector stressDev = stress;
		double I1 = stress(0)+stress(1)+stress(2);
			
		for (int i=0; i<3; i++)
			 stressDev(i) -= 1.0/3*I1;

		double normS = pow( stressDev && stressDev, 0.5); //stressDev.Norm();

		double Fc = CapSurface(normS, I1, hardening_k); 

		result = -R*R*Fc*Fc+(I1-CapBoundL(hardening_k))*(I1-CapBoundL(hardening_k));
		result /= pow(R,4)*pow(Fc,3);


	}

	return result;




};


int 
CapPlasticity::computeConsistentTangent(double deltaGammar1, double deltaGammar2, double deltaGammar3, int mode){
 
    Matrix elasticTangent(6,6);
    Matrix invElastMatrix(6,6);

    Matrix Zig(6,6);
   Zig.Zero();

   elasticTangent.Zero();
   for(int i=0; i<3; i++)
	   for (int j=0; j<3; j++)
		   elasticTangent(i,j) = bulkModulus-2.0/3.0*shearModulus;

   for(int i=0; i<6; i++)
		   elasticTangent(i,i) += 2.0*shearModulus;


       // invElastMatrix = inverse(elasticTangent);
	   invElastMatrix.Zero();
	   elasticTangent.Invert(invElastMatrix);   

// ---------------------------------
   if ((mode ==6) || ((deltaGammar1==0)&&(deltaGammar2==0)&&(deltaGammar3==0))){

	theTangent = elasticTangent;
   
   }
// ----------------------------------

   else if (mode ==3){

	
	   tempMatrix = dF2dSigma(3);

	   Zig = invElastMatrix;
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar2);
       // tempMatrix = inverse(Zig);	   
	   Zig.Invert(tempMatrix);
	   Zig = tempMatrix;



	    Matrix a(2,2);

	    Vector thedFdSigma(6);
	   thedFdSigma = dFdSigma(mode);
	   
	   a(0,0) = tripleTensorProduct (thedFdSigma,Zig,thedFdSigma);

	    Vector thedF2dSigmadk(6);
	   thedF2dSigmadk= dF2dSigmadk();

	   a(0,1) = tripleTensorProduct (thedFdSigma,Zig,thedF2dSigmadk)-1.0/deltaGammar2*dFdk(1);

// ---
	 static Vector stressDev(6);
		stressDev = stress;
	 double I1 = stress(0)+stress(1)+stress(2);
	    
	 for (int i=0; i<3; i++)
		 stressDev(i) -= 1.0/3*I1;

	 double normS = pow( stressDev && stressDev, 0.5); //stressDev.Norm();

 	 double Fc = CapSurface(normS, I1, hardening_k); 
     double dFdI1 =  (I1-CapBoundL(hardening_k))/Fc/R/R;

	 a(1,0) = tripleTensorProduct (thedF2dSigmadk,Zig,thedFdSigma)+1.0/deltaGammar2*dFdI1;

// ---

	 a(1,1) = tripleTensorProduct(thedF2dSigmadk, Zig, thedF2dSigmadk) + 1.0/deltaGammar2*dFdIdk()-1.0/3.0/deltaGammar2/deltaGammar2*dHdk(hardening_k);

	 static Matrix invA(2,2);

     // invA = inverse(a);	   
	 a.Invert(invA);


	
	 static Vector N0(6);
	 static Vector N1(6);

	 N0.addMatrixVector(0.0, Zig, thedFdSigma,1.0); 
	 N1.addMatrixVector(0.0, Zig, thedF2dSigmadk,1.0); 

	 theTangent.Zero();
	 

	 for ( int m=0; m<6; m++){	 
	   for ( int n=0; n<3; n++)				 
	     theTangent(m,n) +=invA(0,0)*N0(m)*N0(n);
	   for (int  n=3; n<6; n++)				 
	     theTangent(m,n) +=invA(0,0)*N0(m)*N0(n)*2.0;
	 }

	 for (int m=0; m<6; m++){	 
	   for ( int n=0; n<3; n++)				 
	     theTangent(m,n) +=invA(0,1)*N0(m)*N1(n);
	   for (int n=3; n<6; n++)				 
			 theTangent(m,n) +=invA(0,1)*N0(m)*N1(n)*2.0;
	 }
	 for (int m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) +=invA(1,0)*N1(m)*N0(n);
		 for (int  n=3; n<6; n++)				 
			 theTangent(m,n) +=invA(1,0)*N1(m)*N0(n)*2.0;
	 }
	 for (int  m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) +=invA(1,1)*N1(m)*N1(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) +=invA(1,1)*N1(m)*N1(n)*2.0;
	 }



	theTangent.addMatrix(-1.0, Zig, 1.0);    
   
   }
//----------------------------------
   else if (mode ==1){

	   tempMatrix.Zero();

	   Zig = invElastMatrix;
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar3);
       // tempMatrix = inverse(Zig);	   
	   Zig.Invert(tempMatrix);
	   Zig = tempMatrix;

	   static Vector N3(6);	
	   static Vector thedFdSigma(6);

	   thedFdSigma = dFdSigma(mode);
	   
	   N3.addMatrixVector(0.0, Zig, thedFdSigma,1.0); 

	   double g33;   
	   
	   g33 = tripleTensorProduct (thedFdSigma,Zig,thedFdSigma);

	   theTangent.addMatrix(0.0, Zig, 1.0);    
	 
	   for ( int m=0; m<6; m++){	    
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=1.0/g33*N3(m)*N3(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=1.0/g33*N3(m)*N3(n)*2.0;
	   }

   
   }
   else if(mode ==5){

	   tempMatrix = dF2dSigma(mode);

	   Zig = invElastMatrix;
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar1);
       // tempMatrix = inverse(Zig);	   
	   Zig.Invert(tempMatrix);
	   Zig = tempMatrix;

	   static Vector N1(6);	
	   static Vector thedFdSigma(6);

	   thedFdSigma = dFdSigma(mode);
	   
	   N1.addMatrixVector(0.0, Zig, thedFdSigma,1.0); 

	   double g11;   
	   
	   g11 = tripleTensorProduct (thedFdSigma,Zig,thedFdSigma);

	   theTangent.addMatrix(0.0, Zig, 1.0);    
	 
	   for ( int m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=1.0/g11*N1(m)*N1(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=1.0/g11*N1(m)*N1(n)*2.0;
	   }
   
   }
//----------------------------------
   else if (mode ==2){

	   static Vector thedFdSigma1(6);
	   static Vector thedFdSigma3(6);

	   thedFdSigma1 = dFdSigma(5);     // dF1/dSigma  --- mode 5
	   thedFdSigma3 = dFdSigma(1);     // dF3/dSigma  --- mode 1

	   tempMatrix = dF2dSigma(5);  // mode=5

	   Zig = invElastMatrix;
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar1);

	   tempMatrix = dF2dSigma(1);  // mode=1
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar3);

       // tempMatrix = inverse(Zig);	   
	   Zig.Invert(tempMatrix);
	   Zig = tempMatrix;

   double g11, g13, g31, g33, invG11,invG13, invG31, invG33; 


   g11 = tripleTensorProduct (thedFdSigma1,Zig,thedFdSigma1);
   g13 = tripleTensorProduct (thedFdSigma1,Zig,thedFdSigma3);
   g31 = tripleTensorProduct (thedFdSigma3,Zig,thedFdSigma1);
   g33 = tripleTensorProduct (thedFdSigma3,Zig,thedFdSigma3);


   double det = g11*g33-g13*g31;

   invG11 =  g33/det;
   invG13 = -g13/det;
   invG31 = -g31/det;
   invG33 =  g11/det;



   static Vector N1(6);   
   static Vector N3(6);   

   N1.addMatrixVector(0.0, Zig, thedFdSigma1,1.0); 
   N3.addMatrixVector(0.0, Zig, thedFdSigma3,1.0); 

   theTangent = Zig;
	 

	 for ( int m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG11*N1(m)*N1(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=invG11*N1(m)*N1(n)*2.0;
	 }

	 for (int  m=0; m<6; m++){	 
		 for (int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG13*N1(m)*N3(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=invG13*N1(m)*N3(n)*2.0;
	 }
	 for (int  m=0; m<6; m++){	 
		 for (int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG31*N3(m)*N1(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=invG31*N3(m)*N1(n)*2.0;
	 }
	 for (int  m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG33*N3(m)*N3(n);
		 for (int  n=3; n<6; n++)				 
			 theTangent(m,n) -=invG33*N3(m)*N3(n)*2.0;
	 }

  
   }

//----------------------------------
   else if (mode ==4){

	   static Vector thedFdSigma1(6);
	   static Vector thedFdSigma2(6);

	   thedFdSigma1 = dFdSigma(5);     // dF1/dSigma  --- mode 5
	   thedFdSigma2 = dFdSigma(3);     // dF2/dSigma  --- mode 3

	   tempMatrix = dF2dSigma(5);  // mode=5

	   Zig = invElastMatrix;
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar1);

	   tempMatrix = dF2dSigma(3);  // mode=3
	   Zig.addMatrix(1.0, tempMatrix, deltaGammar2);

       // tempMatrix = inverse(Zig);	   
	   Zig.Invert(tempMatrix);
	   Zig = tempMatrix;

	   double g11, g12, g21, g22, invG11,invG12, invG21, invG22; 


	   g11 = tripleTensorProduct (thedFdSigma1,Zig,thedFdSigma1);
	   g12 = tripleTensorProduct (thedFdSigma1,Zig,thedFdSigma2);
	   g21 = tripleTensorProduct (thedFdSigma2,Zig,thedFdSigma1);
	   g22 = tripleTensorProduct (thedFdSigma2,Zig,thedFdSigma2);


	   double det = g11*g22-g12*g21;

	   invG11 =  g22/det;
	   invG12 = -g12/det;
	   invG21 = -g21/det;
	   invG22 =  g11/det;



	   static Vector N1(6);   
	   static Vector N2(6);   

	   N1.addMatrixVector(0.0, Zig, thedFdSigma1,1.0); 
	   N2.addMatrixVector(0.0, Zig, thedFdSigma2,1.0); 

	   theTangent = Zig;
	 

	 for ( int m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG11*N1(m)*N1(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=invG11*N1(m)*N1(n)*2.0;
	 }

	 for (int  m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) -=invG12*N1(m)*N2(n);
		 for (int  n=3; n<6; n++)				 
			 theTangent(m,n) -=invG12*N1(m)*N2(n)*2.0;
	 }
	 for (int m=0; m<6; m++){	 
		 for (int  n=0; n<3; n++)				 
			 theTangent(m,n) -=invG21*N2(m)*N1(n);
		 for (int  n=3; n<6; n++)				 
			 theTangent(m,n) -=invG21*N2(m)*N1(n)*2.0;
	 }
	 for (int  m=0; m<6; m++){	 
		 for (int n=0; n<3; n++)		  		 
			 theTangent(m,n) -=invG22*N2(m)*N2(n);
		 for (int n=3; n<6; n++)				 
			 theTangent(m,n) -=invG22*N2(m)*N2(n)*2.0;
	 }	   

   } //mode =4

   
   for (int i=0; i<6; i++)
	  for (int j=3; j<6; j++)
		  theTangent (i,j) /=2.0;

   return 0; 

}

///////////////////////// add sensitivity ///////////////////////////////////
int 
CapPlasticity::setParameter (const char **argv, 
											  int argc, Information &info)
{	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"G") == 0) {
		info.theType = DoubleType;
		return 1;
	}

	if (strcmp(argv[0],"K") == 0) {
		info.theType = DoubleType;
		return 2;
	}

	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		return 3;
	}

	if (strcmp(argv[0],"X") == 0) {
		info.theType = DoubleType;
		return 4;
	}

	if (strcmp(argv[0],"D") == 0) {
		info.theType = DoubleType;
		return 5;
	}

	if (strcmp(argv[0],"W") == 0) {
		info.theType = DoubleType;
		return 6;
	}

	if (strcmp(argv[0],"R") == 0) {
		info.theType = DoubleType;
		return 7;
	}

	if (strcmp(argv[0],"lambda") == 0) {
		info.theType = DoubleType;
		return 8;
	}
	if (strcmp(argv[0],"theta") == 0) {
		info.theType = DoubleType;
		return 9;
	}

	if (strcmp(argv[0],"beta") == 0) {
		info.theType = DoubleType;
		return 10;
	}

	if (strcmp(argv[0],"alpha") == 0) {
		info.theType = DoubleType;
		return 11;
	}

	if (strcmp(argv[0],"T") == 0) {
		info.theType = DoubleType;
		return 12;
	}
	
	else
		opserr << "WARNING: Could not set parameter in CapPlasticity. " << endln;
                
	return -1;
}
	
int 
CapPlasticity::updateParameter(int passedParameterID, 
												 Information &info){

	switch (passedParameterID) {
	case -1:
		return -1;

	case 1:
		this->shearModulus= info.theDouble; // 
		break;

	case 2:
		this->bulkModulus  = info.theDouble;
		break;

	case 3:
		this->rho= info.theDouble; // 
		break;

	case 4:
		this->X  = info.theDouble;
		break;


	case 5:
		this->D= info.theDouble; // 
		break;

	case 6:
		this->W  = info.theDouble;
		break;

	case 7:
		this->R= info.theDouble; // 
		break;

	case 8:
		this->lambda  = info.theDouble;
		break;

	case 9:
		this->theta= info.theDouble; // 
		break;

	case 10:
		this->beta  = info.theDouble;
		break;

	case 11:
		this->alpha= info.theDouble; // 
		break;

	case 12:
		this->T  = info.theDouble;
		
		// -------------2012,1,15 --------
	    if (T>0) 
			this->T = -T;

		break;


	default:
		return -1;
	}
    
 
	return 0;
	}





int CapPlasticity::activateParameter(int passedParameterID){
	parameterID = passedParameterID;
	return 0;						
	}


const Vector & 
CapPlasticity::getStressSensitivity  (int gradNumber, bool conditional){

   opserr<<"Fatal: DDM for Cap Model is not implemented yet!"<<endln;
    exit(-1);
}
