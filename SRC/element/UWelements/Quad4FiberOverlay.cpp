////////////////////////////////////////////////////////////////////////////
//Written: M. Chiaramonte,  P. Arduino,  P.Mackenzie-Helnwein, U.Washington
//			03.29.2011
//Modified: 

// Decription: This file contains the implementation for the Quad4FiberOverlay

#include <Quad4FiberOverlay.h>

double Quad4FiberOverlay::pts[4][2];
double Quad4FiberOverlay::wts[2];
Matrix Quad4FiberOverlay::FiberK(8,8);
Vector Quad4FiberOverlay::P(8);
static int num_Quad4FiberOverlay = 1;
	 

void *
OPS_Quad4FiberOverlay(void)  
{
if (num_Quad4FiberOverlay == 0) {
	num_Quad4FiberOverlay++;
    OPS_Error("Quad4FiberOverlay element - Written: M.Chiaramonte, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n", 1);
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 9) {
    opserr << "Want: Quad4FiberOverlay ?tag ?nd1 ?nd2 ?nd3 ?nd4 ?matTag ?AreaFiber B1  B2 \n";
    return 0;
  }
    
  int    iData[6];
  double dData[3];
  int matTag = 0;
  int eleTag = 0;
  int numData = 0;
  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Quad4FiberOverlay" << endln;
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Quad4FiberOverlay: invalid matTag for element: " << eleTag << "\n";
    return 0;
  }
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Quad4FiberOverlay " << eleTag << endln;
    return 0;	
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matTag);
  
  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element " << eleTag << endln;
    return 0;

  }

  // Parsing was successful, allocate the material
  theElement = new Quad4FiberOverlay(iData[0], iData[1], iData[2], iData[3], iData[4], *theMaterial, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Quad4FiberOverlay\n";
    return 0;
  }

  return theElement;

  opserr << "Element Passed First Inspection" << endln;
}


//first constructor ////////////////////////////////////////////////////////////////////////////
Quad4FiberOverlay::Quad4FiberOverlay(int tag, int nd1, int nd2, int nd3, int nd4,
UniaxialMaterial &m, double AreaFiber,double B1, double B2)
:Element(tag, ELE_TAG_Quad4FiberOverlay)
 ,theMaterial(0)
 ,externalNodes(SL_NUM_NODE)
 ,g1(SL_NUM_NDF)
 ,dualg1(SL_NUM_NDF)
 ,g2(SL_NUM_NDF)
 ,dualg2(SL_NUM_NDF)
 ,beta1(B1)
 ,beta2(B2)
 ,dNidxAlphai(SL_NUM_NODE,SL_NUM_NDF)
 ,Q1(SL_NUM_NDF) ,Q2(SL_NUM_NDF),Q3(SL_NUM_NDF),Q4(SL_NUM_NDF)
 ,Qfi(SL_NUM_NDF)
 ,Qfj(SL_NUM_NDF)
 ,Vf(SL_NUM_NDF)
 ,A(3)
 ,AA(3)
 ,Bb(SL_NUM_DOF)
 ,Af(AreaFiber)
 ,u(SL_NUM_DOF)
{	
	//determine int pts location based on fiber orientation
	nFi.Zero();
	nFj.Zero();
	A.Zero();	
	AA.Zero();
	nFi[0] = -1.0;
	nFi[1] = (beta1 - 0.5) * 2.0;
	nFj[0] = 1.0;
	nFj[1] = (beta2 - 0.5) * 2.0;
	A = nFj - nFi; 
	A.Normalize();
	AA[0] = A(0)*A(0);
	AA[1] = A(1)*A(1);
	AA[2] = A(1)*A(0);  // this must be 1* A1 A2 since strain uses gamma = 2 eps12
	//set up integration paramaters (2 intgr pts)
	pts[0][0] = nFi(0)+A(0)*(1-0.5773502691896258);   
	pts[0][1] = nFi(1)+A(1)*(1-0.5773502691896258);					    
	pts[1][0] = nFj(0)-A(0)*(1-0.5773502691896258);    
	pts[1][1] = nFj(1)-A(1)*(1-0.5773502691896258);					   

	wts[0] = 1.0;
	wts[1] = 1.0;
	externalNodes(0) = nd1;
	externalNodes(1) = nd2;
	externalNodes(2) = nd3;
	externalNodes(3) = nd4;
	theMaterial = m.getCopy();
	for(int i = 0; i<4;i++){ 
	theNodes[i] = 0;
	}

}
//second constructor ////////////////////////////////////////////////////////////////////////////	
Quad4FiberOverlay::Quad4FiberOverlay()
:Element(0, ELE_TAG_Quad4FiberOverlay)
 ,theMaterial(0)
 ,externalNodes(0)
 ,g1(0)
 ,dualg1(0)
 ,g2(0)
 ,dualg2(0)
 ,AA(0)
 ,Bb(0)
 ,dNidxAlphai(SL_NUM_NDF,SL_NUM_NODE)
{	
	opserr << "WARNING material with tag " << endln;
}	
//destructor //////////////////////////////////////////////////////////////////////////////////		
Quad4FiberOverlay::~Quad4FiberOverlay()
{

	}
// // domain 	///////////////////////////////////////////////////////////////////////////////////
void
Quad4FiberOverlay::setDomain(Domain *theDomain)	
{
	 // Check Domain is not null - invoked when object removed from a domain
     if (theDomain == 0) {
	 theNodes[0] = 0;
	 theNodes[1] = 0;
	 theNodes[2] = 0;
	 theNodes[3] = 0;
	 return;
     }

     int Nd1 = externalNodes(0);
     int Nd2 = externalNodes(1);
     int Nd3 = externalNodes(2);
     int Nd4 = externalNodes(3);
	
     theNodes[0] = theDomain->getNode(Nd1);
     theNodes[1] = theDomain->getNode(Nd2);
     theNodes[2] = theDomain->getNode(Nd3);
     theNodes[3] = theDomain->getNode(Nd4);

	 Q1 = theNodes[0]->getCrds();
	 Q2 = theNodes[1]->getCrds();
	 Q3 = theNodes[2]->getCrds();
	 Q4 = theNodes[3]->getCrds();


	 Qfi = Q1 + (Q4 - Q1)*beta1;
	 Qfj = Q2 + (Q3 - Q2)*beta2;
	 Vf = Qfj - Qfi;
	 Lf = Vf.Norm();

	 this->DomainComponent::setDomain(theDomain);

 }
 // nodes 	///////////////////////////////////////////////////////////////////////////////////
 int
Quad4FiberOverlay::getNumExternalNodes(void) const
 {
     return SL_NUM_NODE;
 }
 const ID &
Quad4FiberOverlay::getExternalNodes(void)
 {
   return externalNodes;
 }
 
 Node **
Quad4FiberOverlay::getNodePtrs(void)
 {
   return theNodes;
 }
 
 int
Quad4FiberOverlay::getNumDOF(void) {
     return SL_NUM_DOF;
 }
// // current state //////////////////////////////////////////////////////////////////////////////
int
Quad4FiberOverlay::commitState()
{
    return theMaterial->commitState();
}
 
int
Quad4FiberOverlay::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}
 
int
Quad4FiberOverlay::revertToStart()
{
    return theMaterial->revertToStart();
}
int
Quad4FiberOverlay::update()
{
   strain = 0;
   strain = this->computeCurrentStrain();
  // set the strain in the materials                                                  
  theMaterial->setTrialStrain(strain);
  this->getTangentStiff();
  return 0;
}

double
Quad4FiberOverlay::computeCurrentStrain()
{  // determine the current strain given trial displacements at nodes   
	u.Zero();
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	u[0] = disp1(0);
	u[1] = disp1(1);
	u[2] = disp2(0);
	u[3] = disp2(1);
	u[4] = disp3(0);
	u[5] = disp3(1);
	u[6] = disp4(0);
	u[7] = disp4(1);
	// Loop over the integration points
    strain = 0;
	for(int ip = 0; ip < 2; ip++) {		        //This llop calculates twice the strain, since it is adding the strain at two GP...
		this->getEltBb(pts[ip][0],pts[ip][1]);
		for (int i = 0; i < SL_NUM_DOF; i++) {
			strain += Bb(i)*u(i);
		}
	}
	strain = strain/2.0; //Since it was calculated twice. this function should be redeveloped
	return strain;
}
// // tangent stiff //////////////////////////////////////////////////////////////////////////////
int
Quad4FiberOverlay::Dual()
{	
	dualg1.Zero();
	dualg2.Zero();
//  g1Wedgeg2 = g1*g2' - g2*g1';
//	dualg1 = (g1Wedgeg2*g2)/(g1' * g1Wedgeg2 * g2);
	dualg1(0) = g2(1)/(  g1(0)*g2(1)-g1(1)*g2(0)  );
	dualg1(1) = -g2(0)/(g1(0)*g2(1)-g1(1)*g2(0));
//	dualg2 = (g1'*g1Wedgeg2)'/(g1' * g1Wedgeg2 * g2);
	dualg2(0) = -g1(1)/(  g1(0)*g2(1)-g1(1)*g2(0)   );
	dualg2(1) = g1(0)/(  g1(0)*g2(1)-g1(1)*g2(0)  );
//	opserr << "Dual g1 " << dualg1 << endln;
//	opserr << "Dual g2 " << dualg2 << endln;
	 return 0;
}

int
Quad4FiberOverlay::UpdateBase(double Xi, double Eta)
{	
	Matrix xAlphai(SL_NUM_NDF,SL_NUM_NODE);
	xAlphai.Zero();
	dNidxAlphai.Zero();
	g1.Zero();
	g2.Zero();
	xAlphai(0,0) = -1;
	xAlphai(0,1) = +1;
	xAlphai(0,2) = +1;
	xAlphai(0,3) = -1;
	xAlphai(1,0) = -1;
	xAlphai(1,1) = -1;
	xAlphai(1,2) = +1;
	xAlphai(1,3) = +1;
	Vector Qi;
	for(int i=0; i<4; i++){
		Qi.Zero(); 
		dNidxAlphai(i,0) = xAlphai(0,i)*(1+Eta*xAlphai(1,i))/4.0;
		dNidxAlphai(i,1) = xAlphai(1,i)*(1+Xi*xAlphai(0,i))/4.0;
		Qi = theNodes[i]->getCrds();
		g1 += dNidxAlphai(i,0)*Qi;
		g2 += dNidxAlphai(i,1)*Qi;
	}
    return 0;
}

int
Quad4FiberOverlay::getEltBb(double Xi, double Eta)
{		
		Matrix B(SL_NUM_NDF,SL_NUM_NODE);
		B.Zero();
		Bb.Zero();
		this->UpdateBase(Xi,Eta);
		this->Dual();
		for(int i=0; i<4; i++) {
				for(int j =0; j<2;j++){
				    B(j,i) = dNidxAlphai(i,0)*dualg1(j) + dNidxAlphai(i,1)*dualg2(j); 
				}
				//B(0,i) = dNidxAlphai(i,0)*dualg1(0) + dNidxAlphai(i,1)*dualg1(1); 
				//B(1,i) = dNidxAlphai(i,0)*dualg2(0) + dNidxAlphai(i,1)*dualg2(1); 
				//opserr << "B = " << B(0,i) << B(1,i) << endln;
		}

		for(int i =0; i<4;i++) {
			Bb[(i+1)*2-2] = AA(0)*B(0,i)+AA(2)*B(1,i);
			Bb[(i+1)*2-1] = AA(1)*B(1,i)+AA(2)*B(0,i);
		}
		return 0;
}

const Matrix&
Quad4FiberOverlay::getTangentStiff()
{
		FiberK.Zero();
		double Ef = theMaterial->getTangent();
		for (int ip = 0; ip<2; ip++) { 
			this->getEltBb(pts[ip][0],pts[ip][1]);
			for (int i = 0; i < SL_NUM_DOF; i++){
				for (int j = 0; j < SL_NUM_DOF; j++){
					FiberK(i,j) += Lf/2.0*Af*Ef*wts[ip]*Bb(i) * Bb(j);
					//opserr << "K = " << "(" <<i<<","<<j<<") ="<<FiberK(i,j) << endln;
				}
			}
		}

	return FiberK;
}

const Matrix &
Quad4FiberOverlay::getInitialStiff(void)
{
return getTangentStiff();
}
// // resisting force ////////////////////////////////////////////////////////////////////////////
const Vector&
Quad4FiberOverlay::getResistingForce()
{
	P.Zero();
	//const Vector &disp1 = theNodes[0]->getTrialDisp();
	//const Vector &disp2 = theNodes[1]->getTrialDisp();
	//const Vector &disp3 = theNodes[2]->getTrialDisp();
	//const Vector &disp4 = theNodes[3]->getTrialDisp();
	//Vector u(8);
	//u(0) = disp1(0);
	//u(1) = disp1(1);
	//u(2) = disp2(0);
	//u(3) = disp2(1);
	//u(4) = disp3(0);
	//u(5) = disp3(1);
	//u(6) = disp4(0);
	//u(7) = disp4(1);
	//this->getTangentStiff();
	//P = FiberK ^ u;
//	opserr << "finished resisting force" << endln;
	for (int ip = 0; ip<2; ip++) { 
		this->getEltBb(pts[ip][0],pts[ip][1]);
		for (int i = 0; i < SL_NUM_DOF; i++){
			P(i) += Lf/2.0*Af*wts[ip]*Bb(i)*theMaterial->getStress();
		}
	}
	return P;
}
// // output /////////////////////////////////////////////////////////////////////////////////////

void
Quad4FiberOverlay::Print(OPS_Stream &s, int flag)
{
	if (flag == 2) {
	    s << "#Quad4Overlay\n";
        int i;
		const int numNodes = 4;
		const int nstress = 1 ;
		for (i=0; i<numNodes; i++) {
			const Vector &nodeCrd = theNodes[i]->getCrds();
			const Vector &nodeDisp = theNodes[i]->getDisp();
			s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
		}
    }
}

Response *
Quad4FiberOverlay::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  output.tag("ElementOutput");
  output.attr("eleType",this->getClassType());
  output.attr("eleTag",this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];
 
  for (int i=0; i<numNodes; i++) {
    sprintf(nodeData,"node%d",i+1);
    output.attr(nodeData,nodes(i));
  }
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
		strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
		const Vector &force = this->getResistingForce();
		int size = force.Size();
		for (int i=0; i<size; i++) {
			sprintf(nodeData,"P%d",i+1);
			output.tag("ResponseType",nodeData);
		}
		theResponse = new ElementResponse(this, 1, this->getResistingForce());
  } else if (strcmp(argv[0],"axialForce") ==0) {
		strain = this->computeCurrentStrain();
		theMaterial->setTrialStrain(strain);
		theResponse = new ElementResponse(this, 2, Af*theMaterial->getStress());
  }
   output.endTag();
  return theResponse;
}
 
 
int
Quad4FiberOverlay::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
  case -1:
    return -1;
  case 1: // global forces                                                                                               
    return eleInfo.setVector(this->getResistingForce());
  case 2:
	strain = this->computeCurrentStrain();
	theMaterial->setTrialStrain(strain);
    return eleInfo.setDouble(Af*theMaterial->getStress());
  default:
    return 0;
  }
}

int
Quad4FiberOverlay::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	int dimension = 2;
    static Vector v1(3);
    static Vector v2(3);

    if (displayMode == 1 || displayMode == 2) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();    
	 
	 for (int i=0; i<dimension; i++) {
	    v1(i) = Qfi(i)+end1Disp(i)*fact;
	    v2(i) = Qfj(i)+end2Disp(i)*fact;    
	}
	
	// compute the strain and axial force in the member
	double strain, force;
	if (Lf == 0.0) {
	    strain = 0.0;
	    force = 0.0;
	} else {
	    strain = this->computeCurrentStrain();
	    theMaterial->setTrialStrain(strain);
	    force = Af*theMaterial->getStress();    
	}
    
	if (displayMode == 2) // use the strain as the drawing measure
	  return theViewer.drawLine(v1, v2, (float)strain, (float)strain);	
	else { // otherwise use the axial force as measure
	  return theViewer.drawLine(v1,v2, (float)force, (float)force);
	}
    } else if (displayMode < 0) {
      int mode = displayMode  *  -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = Qfi(i) + eigen1(i,mode-1)*fact;
	  v2(i) = Qfj(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = Qfi(i);
	  v2(i) = Qfj(i);
	}    
      }
      return theViewer.drawLine(v1, v2, 1.0, 1.0);	
    }
    return 0;
}

// // database & parallel ////////////////////////////////////////////////////////////////////////
int
Quad4FiberOverlay::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
return 0;
}
int
Quad4FiberOverlay::sendSelf(int commitTag, Channel &theChannel)
{
return 0;
}


