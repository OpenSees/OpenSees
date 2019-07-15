/*/ 
Written by: Quan Gu_1,  Yongdou Liu_1, Wei Guo_23, Weiquan Li_1, Zhiwu Yu_23, Lizhong Jiang_23 and Hanyun Liu_2 
(1.School of Architecture and Civil Engineering, Xiamen University, 361005, China;
 2.School of Civil Engineering, Central South University, 410075, China;
 3.National Engineering Laboratory for High-speed Railway Construction, 410075, China)

 Reference: Quan Gu, Yongdou Liu, et al. A Practical Wheel-Rail Interaction Element for 
 Modelling Vehicle-Track Systems[J]. International Journal of Structural Stability and Dynamics

 Created: 09/2017
 revised : on 10/28/2018 by Yongdou Liu
 Copyright by the writers. */

#include "WheelRail.h"
#include <CrdTransf.h>
#include <Information.h>
#include <ElementResponse.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

#include <cmath>

#include <elementAPI.h>

//--------------------
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <ElementalLoad.h>
#include <Parameter.h>
#include <iostream>
using namespace std;

Vector WheelRail::contactData(7);
Vector WheelRail::localActiveForce(5);
Vector WheelRail::activeData(7);

WheelRail::WheelRail(int pTag, double pDeltT, double pVel, double pInitLocation, int pNd1, 
		double pRWheel,double pI,double pE,double pA,CrdTransf *pCoordTransf,int pnLoad,
		Vector * pNodeList,	Vector * pDeltaYList,Vector * pDeltaYLocationList)
 :Element(pTag,ELE_TAG_WheelRail),P(0),
   theTangent(0),connectedExternalNodes(0),activeDof(5),
   rearRailNode(2),frontRailNode(2),railDisp(3), shapFun1(2),shapFun2(4)
{
//-----------members in the construtor list----------------------- 
	deltT = pDeltT;
	vel = pVel;   
	initLocation = pInitLocation;
	wheelNodeNum = pNd1;   // 轮节点编号

	rollingRadiusWheel=pRWheel;
	I=pI;
	E=pE;
	A=pA;//     useless and should be deleted.

	theCoordTransf= pCoordTransf;
	//theCoordTransf = theCoordTransf.getCopy2d();
	nLoad=pnLoad;

	if (pNodeList !=0){
		theNodeList = new Vector(*pNodeList);
	}
	if ((pDeltaYList !=0)&&((pDeltaYLocationList !=0))){
		theDeltaYList = new Vector(*pDeltaYList);
		theDeltaYLocationList=new Vector(*pDeltaYLocationList);
	}

	numRailNodeList=pNodeList->Size();
	theNumOfDeltaYList=(*theDeltaYList).Size();

	this->connectedExternalNodes.resize(numRailNodeList+1);
	connectedExternalNodes(0) = pNd1;
	for(int i=1;i<numRailNodeList+1;i++)
		connectedExternalNodes(i) =int( (*theNodeList)(i-1) );

	this->P=new Vector((numRailNodeList+1)*3);
	P->Zero();

	this->theTangent = new Matrix( (numRailNodeList+1)*3,(numRailNodeList+1)*3);
	theTangent->Zero();

	currentLocation = initLocation;
	this->getDeltaY();

	Fhz = 0.0;
	G=4.57*1.0e-8*pow(rollingRadiusWheel,-0.149);
	
	deltaU = 0;
	uF = 0;

	loadStep  = 1;
}

WheelRail::~WheelRail()
{
	/*/
	if (theNumOfNodeList !=0) {
		delete theNodeList;
		theNodeList = 0;
	}

	if (theNumOfDeltaYList !=0) {
		delete theDeltaYList;
		theDeltaYList = 0;
		delete theDeltaYLocationList;
		theDeltaYLocationList = 0;
	}//*/
}

int
WheelRail::getNumExternalNodes(void) const
{
    return (this->numRailNodeList+1);
}

const ID &
WheelRail::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
WheelRail::getNodePtrs(void)
{
  return theNodes;
}

int
WheelRail::getNumDOF(void) 
{
    return (this->numRailNodeList+1)*3;
}

void	
WheelRail::setDomain(Domain *thePassedDomain)
{
	theDomain = thePassedDomain;
	this->DomainComponent::setDomain(theDomain);//将信息存入domain

	theNodes = new Node *[numRailNodeList+1];
	for(int i=0;i<numRailNodeList+1;i++){
		theNodes[i]=theDomain->getNode(connectedExternalNodes(i));
	}

	activeBeamIndex=0;
	frontRailNode = 	theNodes[activeBeamIndex+2]->getCrds();
	while((activeBeamIndex<numRailNodeList-2)&&(currentLocation>frontRailNode(0)+1.e-14)){
		activeBeamIndex++;
		frontRailNode = 	theNodes[activeBeamIndex+2]->getCrds();
	}
	rearRailNode = 	theNodes[activeBeamIndex+1]->getCrds();

	this->getDeltaY();
	
	this->getShapeFuns();

	this->getActiveDof();

}

int
WheelRail::commitState()
{
	loadStep++;

	if(loadStep>nLoad){
		currentLocation = currentLocation + deltT*vel;
		this->getDeltaY();

		frontRailNode = 	theNodes[activeBeamIndex+2]->getCrds();
		while((activeBeamIndex<numRailNodeList-2)&&(currentLocation>frontRailNode(0)+1.e-14)){
			activeBeamIndex++;
			frontRailNode = 	theNodes[activeBeamIndex+2]->getCrds();
		}
		rearRailNode = 	theNodes[activeBeamIndex+1]->getCrds();

		if(activeBeamIndex>numRailNodeList-2){
			opserr<<"the location of the wheel is "<<currentLocation<<
				" which is larger than the front element node frontRailNode "<<frontRailNode(0)<<endln;
			exit(-1);
		}
		
		this->getShapeFuns();

		this->getActiveDof();

	}

	int retVal=this->Element::commitState();
	if (retVal<0) {
		opserr << "WheelRail::commitState() - failed in base class\n";
		return retVal;
	}
	
  return retVal;
}

int
WheelRail::revertToLastCommit()
{
  return 0;
}

int
WheelRail::revertToStart()
{
  return 0;
}

int
WheelRail::update(void)
{
	//-------------------------------------------get element resisting force-----------------------------
	const Vector &disp1 = theNodes[0]->getTrialDisp();//轮节点的位移
	const Vector &disp2 = theNodes[activeBeamIndex+1]->getTrialDisp();//轮下梁节点位移
	const Vector &disp3 = theNodes[activeBeamIndex+2]->getTrialDisp();//轮下梁节点位移

	railDisp.Zero();
	railDisp(0)=disp2(0)*shapFun1(0)+disp3(0)*shapFun1(1);
	railDisp(1)=disp2(1)*shapFun2(0)+disp2(2)*shapFun2(1)+disp3(1)*shapFun2(2)+disp3(2)*shapFun2(3);

	Vector limits(2);
	double deltaU0=railDisp(1)-disp1(1)+theDeltaY;

	Fhz=0;
	uF=0;
	deltaU=deltaU0;

	if(deltaU0 >0){
		limits(0)=0;
		limits(1)=pow(deltaU0/G,1.5);
		this->NewtonBisection( limits, disp1(1) );
		//Fhz=FalsePostionAlgorithm(limits, uWheel); //eather one is ok
		uF=Fhz*pow(b,3)*pow(a,3)/3/E/I/theEleLength/theEleLength/theEleLength;
		deltaU=railDisp(1)-uF- disp1(1)+ theDeltaY;
	}

	//---求由于赫兹力产生的等效节点荷载，参考《结构力学1》（第二版）281页，
	P->Zero();
	(*P)(activeDof(0))=-Fhz;
	for(int i=0;i<4;i++)
		(*P)(activeDof(i+1))=Fhz*shapFun2(i);

	//-------------------------get element stiffness---------------------------
	theTangent->Zero();
	if(Fhz>0){

		Vector dRdFH(5),dFHdU(5);

		dRdFH(0) = -1.0;
		for(int i=0;i<4;i++)
			dRdFH(i+1)=shapFun2(i);

		double dDeltadFH = 2*G*pow(Fhz,-1/3.0)/3.0;//注明：暂时放在getUfDeltaU（）更新
		double beamFlexibility = pow(a*b,3.0)/3/E/I/theEleLength/theEleLength/theEleLength;

		dFHdU(0) = -1/(beamFlexibility+dDeltadFH);
		for(int i=0;i<4;i++)
			dFHdU(i+1)=shapFun2(i)/(beamFlexibility+dDeltadFH);

		for(int i = 0;i<5;i++)
			for(int j = 0;j<5;j++)
				(*theTangent)(activeDof(i),activeDof(j)) = dRdFH(i)*dFHdU(j);

	}else{

		if(loadStep<=nLoad){
			*theTangent = this->getInitialStiff();
		}else{
			/*/
			opserr<<"The wheel "<<wheelNodeNum<<" separates with "<<
				activeBeamIndex+1<<"th rail beam element at location of x= "<<currentLocation<<endln;//*/
		}

	}

  return 0;
}


const Matrix & WheelRail::getTangentStiff(void)
{
	return *theTangent;
}


const Matrix &
WheelRail::getInitialStiff(void)   
{
	//An empirical tangent stiffness used when the wheel and rail separate each other in static analysis.
	Matrix KG(9,9),kg(5,5);
	kg.Zero();
	KG(1,1) = 1.64002e+006;
	KG(1,2) = 1639.94;
	KG(1,4) = 0.00305968;
	KG(1,5) = -0.0408973;
	KG(1,7) = -1.64002e+006;

	KG(2,2) = 16419.9;
	KG(2,4) = 0.030635;
	KG(2,5) = -0.409484;
	KG(2,7) = -1.64002e+006;

	KG(4,4) = 0.00114313;
	KG(4,5) = -0.0152797;
	KG(4,7) = -612730;

	KG(5,5) = 0.204237;
	KG(5,7) = 8.19009e+006;

	KG(7,7) =  821077;
	for(int i = 0;i<9;i++){
		for(int j = i;j<9;j++){
			if(i!=j){
				KG(j,i) = KG(i,j);			
			}
		}
	}
	kg(0,0)=KG(7,7);
	kg(0,1)=KG(7,1);
	kg(0,2)=KG(7,2);
	kg(0,3)=KG(7,4);
	kg(0,4)=KG(7,5);

	kg(1,1)=KG(1,1);
	kg(1,2)=KG(1,2);
	kg(1,3)=KG(1,4);
	kg(1,4)=KG(1,5);

	kg(2,2)=KG(2,2);
	kg(2,3)=KG(2,4);
	kg(2,4)=KG(2,5);

	kg(3,3)=KG(4,4);
	kg(3,4)=KG(4,5);

	kg(4,4)=KG(5,5);

	for(int i = 0;i<5;i++){
		for(int j = 0;j<5;j++){
			(*theTangent)(activeDof(i),activeDof(j)) = kg(i,j);
		}
	}
	return *theTangent; 
}
    
void 
WheelRail::zeroLoad(void)
{
  return;  // ok
}

int 
WheelRail::addLoad(const Vector &addP)
{
  return 0;  // ok
}

int 
WheelRail::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;  // ok
}

const Vector & WheelRail::getResistingForce()
{	
	return *P;
}


const Vector &
WheelRail::getResistingForceIncInertia()
{	
    *P = this->getResistingForce();
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    *P += this->getRayleighDampingForces();

  return *P;
}


int
WheelRail::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;   //ok
}

int
WheelRail::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;    //ok
}


int
WheelRail::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
  return 0;   //ok
}


void
WheelRail::Print(OPS_Stream &s, int flag)
{
  return;    //ok
}


Response*
WheelRail::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;
  //
  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","WheelRail");
  output.attr("eleTag",this->getTag());

  for(int i=0;i<numRailNodeList+1;i++){
	sprintf(outputData,"node%d",i);
	output.attr(outputData,connectedExternalNodes[i]);
  }

  /*/ 记录激活量
      激活自由度（5个） 激活单元号（1个） commited位置（1个）//*/
  if (strcmp(argv[0],"activeData") == 0 || strcmp(argv[0],"activeDatas") == 0) {
	for(int i=0;i<5;i++)	
		activeData(i)=activeDof(i);
	activeData(5)=this->activeBeamIndex;
	activeData(6)=currentLocation;

    theResponse =  new ElementResponse(this, 2, activeData);
  
  /*/ 记录有效单元节点力
      P（5个）//*/
  }    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {
	for(int i=0;i<5;i++) 
		localActiveForce(i)=(*P)(activeDof(i));

    theResponse = new ElementResponse(this, 3, localActiveForce);

  /*/ 记录接触量
      嵌入量（1个）uF（1个） Fhz（1个） theDeltaY（1个） uUnderWheel（3个）//*/
  }    else if (strcmp(argv[0],"contactData") == 0 || strcmp(argv[0],"contactDatas") == 0) {

	contactData(0)=this->deltaU;
	contactData(1)=this->uF;
	contactData(2)=this->Fhz;
	contactData(3)=this->theDeltaY;
	for(int i=0;i<3;i++)
		contactData(i+4)= railDisp(i);
    theResponse = new ElementResponse(this, 4, contactData);

  } // ElementOutput//*/

  return theResponse;
}


int 
WheelRail::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // activeData
 /*/ 记录激活量
      激活自由度（5个） 激活单元号（1个） commited位置（1个）//*/
	for(int i=0;i<5;i++)
		activeData(i)=activeDof(i);
	activeData(5)=this->activeBeamIndex;
	activeData(6)=currentLocation;
    return eleInfo.setVector(activeData);
    
  case 3:
	/*/ 记录有效单元节点力
    P（5个）//*/
	for(int i=0;i<5;i++)
		localActiveForce(i)=(*P)(activeDof(i));
	return eleInfo.setVector(localActiveForce);
    
  case 4:
	  /*/ 记录接触量
      嵌入量（1个）uF（1个） Fhz（1个） theDeltaY（1个） uUnderWheel（3个）//*/
	contactData(0)=this->deltaU;
	contactData(1)=this->uF;
	contactData(2)=this->Fhz;
	contactData(3)=this->theDeltaY;
	for(int i=0;i<3;i++)
		contactData(i+4)= railDisp(i);
    return eleInfo.setVector(contactData);

  default:
    return -1;
  }
}


int
WheelRail::setParameter(const char **argv, int argc, Parameter &param)
{
  return 0;   //ok
}
    

int
WheelRail::updateParameter(int parameterID, Information &info)
{
  return -1;   //ok
}

double WheelRail::getResidualOfDeltaU(double pFhz,double uWheel){
	//The residual of wheel-rail inter-penetration obtained using Geometric equation 
	//and using constitutive equation.  
	uF=pFhz*pow(b,3)*pow(a,3)/3/E/I/theEleLength/theEleLength/theEleLength;
	double deltaU1=railDisp(1)-uF-uWheel+theDeltaY;
	double deltaU2=0;	
	if(pFhz>0){
	  deltaU2=pow(pFhz,0.666666666666667)*G;
	}
	return (deltaU1-deltaU2);
}


void WheelRail::NewtonBisection(Vector limits,double uWheel){
	//The Newton-Bisection algorithm used to solve the Fhz that meets the Geometric equation
	//between wheel and rail.	

	//该牛顿二分法1 可行
	int maxIterT=30;double tol=1.0e-5;
	double FHL=limits(0),FHH=limits(1);
	double FHzi=0.5*(FHL+FHH);
	double dDeltaUi=0,dDeltaU=0;
	double R = pow(a*b,3.0)/3/E/I/theEleLength/theEleLength/theEleLength;
	double dUbaldFH=0;
	int i=0;
	bool converge=false;
	while (i<maxIterT&& converge==false) {
		dDeltaUi=getResidualOfDeltaU(FHzi,uWheel);
		double dUdFH=-2/3.0*G*pow(FHzi,-1/3.0)-R;
		Fhz=FHzi-dDeltaUi/dUdFH;

		if(Fhz>fmax(FHH,FHL)||Fhz<fmin(FHH,FHL)){
			Fhz=0.5*(FHL+FHH);
			dDeltaU=getResidualOfDeltaU(Fhz,uWheel);
			if(dDeltaU==0){
				converge=true;
			}else if(dDeltaU*getResidualOfDeltaU(FHH,uWheel)<0){
				FHL=Fhz;
			}else {
				FHH=Fhz;
			}
		}

		if(abs(Fhz-FHzi)<tol && abs(dDeltaU-dDeltaUi)<tol*1.0e-6) 
			converge=true;
		FHzi=Fhz;

	}
	if (i>maxIterT)	opserr<<maxIterT<<"次迭代后失败！";//previous process*/ 
//===========================
	/*
	int maxIterT=30;double tol=1.0e-5;
	double FHL=limits(0),FHH=limits(1);
	double FHzi=0.5*(FHL+FHH);
	double dDeltaUi=0,dDeltaU=0;
	double R = pow(a*b,3.0)/3/E/I/theEleLength/theEleLength/theEleLength;
	double dUbaldFH=0;
	int i=0;
	bool converge=false;
	while (i<maxIterT&& converge==false) {
		dDeltaUi=getdDeltaU(FHzi,uWheel);
		double dUbaldFH=-2/3.0*G*pow(FHzi,-1/3.0)-R;
		Fhz=FHzi-dDeltaUi/dUbaldFH;

		if((Fhz>FHH)||(Fhz<FHL)){
			Fhz=0.5*(FHL+FHH);
			dDeltaU=getdDeltaU(Fhz,uWheel);
			if(dDeltaU>0){
				FHL=Fhz;
			}else if(dDeltaU<0){
				FHH=Fhz;
			}else{
				converge=true;
			}
		}

		if(abs(Fhz-FHzi)<tol && abs(dDeltaU-dDeltaUi)<tol*1.0e-6) 
			converge=true;
		FHzi=Fhz;

	}
	if (i>maxIterT)	opserr<<maxIterT<<"次迭代后失败！";//previous process*/ 
//========================================
	/*
	int maxIterT=15;double tol=1.0e-6;
	double FhzL=limits(0),FhzH=limits(1);
	double FhzLastIter=0.5*(FhzL+FhzH);
	//FhzLastIter=FhzH*(1-a*b/theEleLength/theEleLength);//
	double PhiLastIter=getdDeltaU(FhzLastIter,uWheel),Phi=0;
	double R = pow(a*b,3.0)/3/E/I/theEleLength/theEleLength/theEleLength;
	double dPhi_dFhz=0;
	int i=1;
	bool converge=false;
	while(i<maxIterT&& converge==false) {
		dPhi_dFhz= -2/3.0*G*pow(FhzLastIter,-1/3.0) - R;
		Fhz=FhzLastIter-PhiLastIter/dPhi_dFhz;

		if(Fhz>max(FhzH,FhzL)||Fhz<min(FhzH,FhzL) )
			Fhz=0.5*(FhzL+FhzH);
		Phi=getdDeltaU(Fhz,uWheel);

		if(Phi*getdDeltaU(FhzH,uWheel)<0){
			FhzL=FhzH;
		}			

		FhzH=Fhz;

		if(abs(Fhz-FhzLastIter)/Fhz<tol && abs(Phi)<tol*1.0e-9) {
			converge=true;
		}
		//opserr<<"i= "<<i<<":  "<<Fhz<<"   "<<Phi<< endln;
		FhzLastIter=Fhz;
		PhiLastIter=getdDeltaU(FhzLastIter,uWheel);

		i++;
	}
	if (i>=maxIterT)	opserr<<maxIterT<<"次迭代后失败！";//revised on2018/9/4 by Yongdou Liu/*/

}

double WheelRail::FalsePostionAlgorithm(Vector limits,double uWheel){
	//The false-position algorithm used to solve the Fhz that meets the Geometric equation
	//between wheel and rail.
	int maxIterT=30;
	double tol=1.0e-6;
	double FHz1=limits(0),FHz2=limits(1);
	double Phi=0,Phi1=0, Phi2=0;
	double dPhi_dFhz=0;
	int i=0;
	double fhz=0;//debug
	for (i=0;i<maxIterT;i++) {
		Phi1=getResidualOfDeltaU(FHz1,uWheel);
		Phi2=getResidualOfDeltaU(FHz2,uWheel);
		if(abs(FHz1-FHz2)<tol && abs(Phi1-Phi2)<tol*1.0e-9){
			fhz=FHz1;
			return fhz;
		}
		else{
			dPhi_dFhz=(Phi2-Phi1)/(FHz2-FHz1);
			fhz=FHz2-Phi2/dPhi_dFhz;
		}
		
		Phi=getResidualOfDeltaU(fhz,uWheel);

		if (abs(FHz2-fhz)/fhz<tol || abs(Phi)<tol*1e-9)	return fhz;

		if(Phi*Phi2<0){
			FHz1=FHz2;
			Phi1=Phi2;
		}
		FHz2=fhz;
		Phi2=Phi;
		//opserr<<"i="<<i<<"  FalsePostionAlgorithm:  fhz=  "<<fhz<<" "<<Phi<<endln;
	}
	if (i>=maxIterT)
		opserr<<maxIterT<<" WHEEL RAIL MAX ITER EXCEEDED";
	return 0.;
}

void WheelRail::getDeltaY(){
	//The rail irregularity
	int i=0; 
	while ((i<theNumOfDeltaYList)&&(currentLocation>(*theDeltaYLocationList)(i)+1.e-14))
		i++;
	if ( (i==0)||(i > theNumOfDeltaYList) ){
		theDeltaY = 0;
	}else { // general case
		theDeltaY  = (*theDeltaYList)(i-1)+(currentLocation-(*theDeltaYLocationList)(i-1))
			*((*theDeltaYList)(i)-(*theDeltaYList)(i-1))/((*theDeltaYLocationList)(i) -(*theDeltaYLocationList)(i-1));
	}
}

void WheelRail::getActiveDof(){
	//The relationship between active DOFs and the active rail beam element index( 1th, 2th or nth).
	activeDof(0)=1;
	activeDof(1)=3*(activeBeamIndex+1)+1;
	activeDof(2)=3*(activeBeamIndex+1)+2;
	activeDof(3)=3*(activeBeamIndex+2)+1;
	activeDof(4)=3*(activeBeamIndex+2)+2;
}

void WheelRail::getShapeFuns(){

	theEleLength = sqrt((rearRailNode(0)-frontRailNode(0))*(rearRailNode(0)-frontRailNode(0))+ (rearRailNode(1)-frontRailNode(1))*(rearRailNode(1)-frontRailNode(1)));
	a=currentLocation-rearRailNode(0);
	b=theEleLength-a;

	double ksi=2*(currentLocation-rearRailNode(0))/theEleLength-1;

	shapFun1(0)=(1-ksi)/2;    shapFun1(1)=(1+ksi)/2;

	shapFun2(0)=0.25*(1-ksi)*(1-ksi)*(2+ksi);    shapFun2(1)=0.125*theEleLength*(1-ksi)*(1-ksi)*(1+ksi);
	shapFun2(2)=0.25*(1+ksi)*(1+ksi)*(2-ksi);   shapFun2(3)=-0.125*theEleLength*(1+ksi)*(1+ksi)*(1-ksi);
}
