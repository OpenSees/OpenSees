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
                                                                        
// Written by: Gonzalo Torrisi, Universidad Nacional de Cuyo
// Input modified 20/19/2022

// we specify what header files we need
#include "BeamGT.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix BeamGT::BeamK(6,6);
Vector BeamGT::BeamR(6);

void* OPS_BeamGT()

{
  // print out a message about who wrote this element & any copyright info wanted
  //if (numMyBeam == 0) {
//	opserr << " \n";
//    opserr << "           Beam with Flexure and Shear Hinges\n";
//    opserr << "   Written by Gonzalo Torrisi UNCuyo Copyright 2016\n";
//	opserr << "                 Only in plane X-Y \n";
//    opserr << "                Use at your Own Peril\n";

  //  numMyBeam++;
  //}

  Element *theBeam = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theBeam = new BeamGT();
    return theBeam;
  }

  if (numRemainingArgs != 10) {
    //opserr << "ERROR - BeamGT not enough args provided, want: element BeamGT tag? Node1? Node2?  matTag? matTag2? matTag3? E? G? A? I? Lp1? Lp2? Lr? fc? \n";
	opserr << "ERROR - BeamGT not enough args provided, want: element BeamGT tag? Node1? Node2?  matTag? matTag2? matTag3? Lp1? Lp2? Lr? fc? \n";
   // numMyBeam++;
  }

  // get the id and end nodes 
  int iData[6];
 double dData[4];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
    int eleTag = iData[0];
 
  numData =1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING error reading element material 1 tag for element " << eleTag << endln;
    return 0;
  }
      numData =1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING error reading element material 2 tag for element " << eleTag << endln;
    return 0;
  }
      numData =1;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING error reading element material 3 tag for element " << eleTag << endln;
    return 0;
  }
   int matID = iData[3];
  int matID2 = iData[4];
  int matID3 = iData[5];
 
  numData =4 ;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
   opserr << "WARNING error reading Elastic properties for element" << eleTag << endln;
   return 0;
  }


  	  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
	  UniaxialMaterial *theMaterial2 =OPS_GetUniaxialMaterial(matID2);
	  UniaxialMaterial *theMaterial3 =OPS_GetUniaxialMaterial(matID3);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element " << eleTag << endln;
    return 0;
  }

    if (theMaterial2 == 0) {
    opserr << "WARNING material with tag " << matID2 << "not found for element " << eleTag << endln;
    return 0;
  }
	   if (theMaterial3 == 0) {
    opserr << "WARNING material with tag " << matID3 << "not found for element " << eleTag << endln;
    return 0;
  }

    // now create the truss and add it to the Domain

  theBeam = new BeamGT(eleTag, iData[1], iData[2], *theMaterial, *theMaterial2, *theMaterial3, dData[0],dData[1],dData[2],dData[3]);

  if (theBeam == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    delete [] theMaterial;
	delete [] theMaterial2;
	delete [] theMaterial3;
	return 0;
  }

  return theBeam;
}


// typical constructor
BeamGT::BeamGT(int tag, 
                   int Nd1, int Nd2,
                   UniaxialMaterial &theMat, UniaxialMaterial &theMat2, UniaxialMaterial &theMat3,
                    double lp1, double lp2, double lr, double fc)
:Element(tag, ELE_TAG_BeamGT),     
 externalNodes(2),
 trans(4,4),LP1(lp1),LP2(lp2),LR(lr),FC(fc), theMaterial(0),theMaterial2(0),theMaterial3(0),Tm(6,6),TTm(6,6),Cdefor(3),Tdefor(3),Cdespla(6),Tdespla(6),Stifloc(6,6),Stif0(6,6),Cesf(3),Tesf(3),RR(3)
 {       
    // allocate memory for numMaterials1d uniaxial material models
  theMaterial = new UniaxialMaterial *[2];
  theMaterial2= 0;
  theMaterial3= 0;


  if ( theMaterial == 0 ) {
    opserr << "FATAL BeamGT::BeamGT - failed to create a 1d  material or direction array\n";
    exit(-1);
  }

  // get a copy of the material and check we obtained a valid copy
 
      theMaterial[0]= theMat.getCopy();
	  theMaterial[1]= theMat.getCopy();
	  theMaterial2= theMat2.getCopy();
	  theMaterial3= theMat3.getCopy();

 for (int i=0; i<2; i++) {
	  if (theMaterial[i] == 0) {
        opserr << "FATAL BeamGT::BeamGT - failed to get a copy of material\n" ;
        exit(-1);
      }
 }
  
   	 if (theMaterial2== 0){
		 opserr << "FATAL BeamGT::BeamGT - failed to get a copy of material2\n";
	 }

   	 if (theMaterial3== 0){
		 opserr << "FATAL BeamGT::BeamGT - failed to get a copy of material3\n";
	 }


  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL BeamGT::BeamGT() - out of memory, could not create an ID of size 2\n";
    exit(-1);
  }

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;   


  theNodes[0] = 0; 
  theNodes[1] = 0;
}

// constructor which should be invoked by an FE_ObjectBroker only
BeamGT::BeamGT()
:Element(0, ELE_TAG_BeamGT),     
 theMaterial(0), 
 theMaterial2(0),
 theMaterial3(0),
 externalNodes(2),
 trans(4,4), LP1(0.0),LP2(0.0), LR(0.0),FC(0.0),Tm(6,6),TTm(6,6),Cdefor(3),Tdefor(3),Cdespla(6),Tdespla(6),Stifloc(6,6),Stif0(6,6),Cesf(3),Tesf(3),RR(3)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
BeamGT::~BeamGT()
{
    // clean up the memory associated with the element, this is
    // memory the Truss2D objects allocates and memory allocated 
    // by other objects that the Truss2D object is responsible for 
    // cleaning up, i.e. the MaterialObject.
	 for (int i=0; i<2; i++) {
	    if (theMaterial[i] != 0)
        delete theMaterial[i];    
	 }
 	delete theMaterial2;
	delete theMaterial3;
}

int
BeamGT::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
BeamGT::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
BeamGT::getNodePtrs(void) 
{
  return theNodes;
}

int
BeamGT::getNumDOF(void) {
    return 6;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.

void
BeamGT::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        return;
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    Node *end1Ptr, *end2Ptr;
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);

    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);  
    

    if (end1Ptr == 0) {
      opserr << "WARNING BeamGT::setDomain() - at Beam " << this->getTag() << " node " <<
        Nd1 << "  does not exist in domain\n";
                                
        return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {        
      opserr << "WARNING BeamGT::setDomain() - at Beam " << this->getTag() << " node " <<
        Nd2 << " does not exist in domain\n";

        return;  // don't go any further - otherwise segemntation fault
    }   

        theNodes[0] = end1Ptr;
		theNodes[1] = end2Ptr;
        
        
        // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();       


    if ((dofNd1 != 3 || dofNd2 != 3)) {
      opserr << "BeamGT::setDomain(): 3 dof required at nodes\n";
      return;
    }   

    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds(); 
     

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1); 
	double  L = sqrt(dx*dx + dy*dy);
            if (L == 0.0) {
      opserr << "WARNING BeamGT::setDomain() - BeamGT " << this->getTag() <<" has zero length\n";
      return;  // don't go any further - otherwise divide by 0 error
    }

        // cosenos de las diagonales
        double coseno=dx/L;
        double seno=dy/L;
		trans(0,0)=L;
		trans(0,1)=coseno;
		trans(0,2)=seno;
		this->revertToStart();
		this->revertToLastCommit();
		
}        

int
BeamGT::commitState()
{
	 int ecode=0;
	 //commit material models
	 	 for (int i=0; i<2; i++) {
		 ecode += theMaterial[i]->commitState();
		 }
	 //commit the base class
	     ecode += theMaterial2->commitState();
	     ecode += theMaterial3->commitState();
	     Cdeltares=Tdeltares;
	for (int i=0; i<3; i++) {
		Cdefor[i]=Tdefor[i];
	}
	for (int i=0; i<6; i++) {
		Cdespla[i]=Tdespla[i];
	}
			for (int i=0; i<3; i++) {
				Cesf[i]=Tesf[i];
			}

				RR[0]=dcur1c;
				RR[1]=dcur2c;
				RR[2]=dgamc;
				RR[3]=daxc;

  return ecode;
}

int
BeamGT::revertToLastCommit()
{
//    return theMaterial->revertToLastCommit();
    int code=0;
    
    // revert state for 1d materials
	for (int i=0; i<2; i++) {
     code += theMaterial[i]->revertToLastCommit();
	}
    code += theMaterial2->revertToLastCommit();
    code += theMaterial3->revertToLastCommit();
	Tdeltares=Cdeltares;
		for (int i=0; i<3; i++) {
		Tdefor[i]=Cdefor[i];
	}
	for (int i=0; i<6; i++) {
		Tdespla[i]=Cdespla[i];
	}
	for (int i=0; i<3; i++) {
		Tesf[i]=Cesf[i];
	}

		dcur1c=RR[0];
		dcur2c=RR[1];
		dgamc=RR[2];
		daxc = RR[3];

    return code;
}

int
BeamGT::revertToStart()
{
  //  return theMaterial->revertToStart();
            int code=0;
    
    // revert state for 1d materials
			for (int i=0; i<2; i++) {
        code += theMaterial[i]->revertToStart();
			}
    code += theMaterial2->revertToStart();
    code += theMaterial3->revertToStart();
	for (int i=0; i<3; i++) {
		Cdefor[i]=0.0;
		Tdefor[i]=0.0;
	}
	for (int i=0; i<6; i++) {
		Cdespla[i]=0.0;
	    Tdespla[i]=0.0;
	}
	 		for (int i=0; i<3; i++) {
				Cesf[i]=0.0;
				Tesf[i]=0.0;
			}

			for (int i=0; i<6; i++){
				for (int j=0; j<6; j++){
					Stifloc(i,j)=0.0;
					Stif0(i,j)=0.0;
				}
			}
			for (int i=0; i<3; i++){
				RR(i)=0.0;
			}

		return code;
}

int
BeamGT::update()
{
  // determine the current strain given trial displacements at nodes
 // double strain = this->computeCurrentStrain();

  // set the strain in the materials
  //theMaterial->setTrialStrain(strain);

        // compute strain and rate; set as current trial for material
 //    strain     = this->computeCurrentStrain(mat );
   // strainRate = this->computeCurrentStrain();

	double str[3];

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();  
	const Vector &Idisp1 =theNodes[0]->getIncrDisp();
	const Vector &Idisp2 =theNodes[1]->getIncrDisp();

     double L=trans(0,0);
     double c1=trans(0,1);
     double	s1=trans(0,2);
	//global displacements
		double dcur1r=0.0;
	   double dcur2r=0.0;
		double ddcur1=0.0;
	    double ddcur2=0.0;
		double ddgam=0.0;
	// 
	// momentos convergidos

     // displacements at the ends of the beam
		double dx1=Idisp1(0);
		double dy1=Idisp1(1);
		double r1=Idisp1(2);
		double dx2=Idisp2(0);
		double dy2=Idisp2(1);
		double r2=Idisp2(2);
		double x1=disp1(0);
		double x2=disp2(0);
		double y1=disp1(1);
		double y2=disp2(1);
		//Displacements at the ends of the beam in local coordinates
		double u1 = dx1*c1+dy1*s1;
		double v1 = -dx1*s1+dy1*c1;
		double u2 = dx2*c1+dy2*s1;
		double v2 = -dx2*s1+dy2*c1;
		double xnorm = 1.0;
		double du1 = x1*c1+y1*s1;
		double du2 = x2*c1+y2*s1;
		 //CURVATURES
		 //incremento de esfuerzos de extremo y curvaturas
		double DM1=0.0;
	    double DM2=0.0;
		double DV=0.0;
		double DN = 0.0;

 		//tengo que traer la matriz de rigidez convergida
		double ei = theMaterial[0]->getInitialTangent();
		double EA = theMaterial3->getInitialTangent();
		double ret = 0;
		
			//matriz de rigidez en locales
			double a1 = theMaterial[0]->getTangent();
			double b1 = theMaterial[0]->getInitialTangent();
			double re1 = a1 / b1;
			double f1 = (1.0 - re1) * LP1 / (re1 * ei);

			double a2 = theMaterial[1]->getTangent();
			double b2 = theMaterial[1]->getInitialTangent();
			double re2 = a2 / b2;
			double f2 = (1.0 - re2) * LP2 / (re2 * ei);

			double a3 = theMaterial2->getTangent();
			double b3 = theMaterial2->getInitialTangent();
			double re3 = a3 / b3;
			double f3 = (1.0 - re3) * LR / (re3 * b3 * L * L);
			// Tangent Flexibility matrix
			double f11 = L / (3.0 * ei) + f1 + FC * 1.2 / (L * b3) + f3;
			double f22 = L / (3.0 * ei) + f2 + FC * 1.2 / (L * b3) + f3;
			double f12 = -L / (6.0 * ei) + FC * 1.2 / (L * b3) + f3;
			double det = f11 * f22 - f12 * f12;
			double fa1 = f22 / det;
			double fb1 = -f12 / det;
			double fc1 = -f12 / det;
			double fd1 = f11 / det;
			double x = 1.0 / L;
			double Am = (fa1 + 2 * fb1 + fd1) * x * x;
			double B = (fa1 + fb1) * x;
			double C = (fb1 + fd1) * x;
			// Momentos convergidos MCE y corte convergido VCE


			DV = Am * (v1 - v2) + B * r1 + C * r2;
			DM1 = (B * (v1 - v2) + fa1 * r1 + fb1 * r2);
			DM2 = C * (v1 - v2) + fb1 * r1 + fd1 * r2;

			// incremento de curvaturas
			double dcur1 = DM1 / theMaterial[0]->getTangent();
			double dcur2 = DM2 / theMaterial[1]->getTangent();
			double dgam = DV / theMaterial2->getTangent();
			double dax = (du2 - du1) / L;

			dcur1c = dcur1c + dcur1 ;
			dcur2c = dcur2c + dcur2 ;
		    dgamc = dgamc + dgam;
			daxc = 0.0;
			daxc = daxc + dax;
			ret = theMaterial3->setTrialStrain(daxc);
			double N = 1.0 * theMaterial3->getStress();

				ret = theMaterial2->setTrialStrain(dgamc, N);
				ret = theMaterial[0]->setTrialStrain(dcur1c, N);
				ret = theMaterial[1]->setTrialStrain(dcur2c, N);

  return ret;
}

const Matrix &
BeamGT::getTangentStiff(void)
{
	   double L=trans(0,0);
        double c1=trans(0,1);
        double	s1=trans(0,2);
		
  	 //Tangent Flexibility matrix 
	 double a1=theMaterial[0]->getTangent();
	 double b1=theMaterial[0]->getInitialTangent();
	 double re1=a1/b1;

	double f1=(1.0-re1)*LP1/(re1*b1);
	 double a2=theMaterial[1]->getTangent();
	 double b2=theMaterial[1]->getInitialTangent();
	 double re2=a2/b2;
	 double f2=(1.0-re2)*LP2/(re2*b2);
	 double a3=theMaterial2->getTangent();
	 double b3=theMaterial2->getInitialTangent();
	 double re3=a3/b3;
	 double f3=(1.0-re3)*LR/(re3*b3*L*L);
	 double EA=theMaterial3->getTangent();

	 if(re3>1.0){
		 re3=1.0;
	 }

	 double f11 = L / (3.0 * b2) + f1 + FC * 1.2 / (L * b3) + f3;
	 double f22 = L / (3.0 * b2) + f2 + FC * 1.2 / (L * b3) + f3;
	 double f12 = -L / (6.0 * b2) + FC * 1.2 / (L * b3) + f3;

	 //TAngent stiffness
	 double det=f11*f22-f12*f12;
	 double fa1=f22/det;
	 double fb1=-f12/det;
	 double fc1=-f12/det;
	 double fd1=f11/det;
	 double x=1.0/L;

	 double Am=(fa1+2*fb1+fd1)*x*x;
	 double B=(fa1+fb1)*x;
	 double C=(fb1+fd1)*x;

	 //Beam tangent stiffness matrix in global coordinates
	 BeamK(0,0)=EA/L*c1*c1+Am*s1*s1;
     BeamK(0,1)=EA/L*c1*s1-Am*c1*s1;
	 BeamK(0,2)=-B*s1;
	 BeamK(0,3)=-EA/L*c1*c1-Am*s1*s1;
	 BeamK(0,4)=-EA/L*c1*s1+Am*s1*c1;
	 BeamK(0,5)=-C*s1;

	 BeamK(1,0)=BeamK(0,1);
	 BeamK(1,1)=EA/L*s1*s1+Am*c1*c1;
	 BeamK(1,2)=B*c1;
	 BeamK(1,3)=-EA/L*s1*c1+Am*s1*c1;
	 BeamK(1,4)=-EA/L*s1*s1-Am*c1*c1;
	 BeamK(1,5)=C*c1;

	 BeamK(2,0)=BeamK(0,2);
	 BeamK(2,1)=BeamK(1,2);
	 BeamK(2,2)=fa1;
	 BeamK(2,3)=B*s1;
	 BeamK(2,4)=-B*c1;
	 BeamK(2,5)=fb1;

	 BeamK(3,0)=-EA/L*c1*c1-Am*s1*s1;
	 BeamK(3,1)=-EA/L*c1*s1+Am*s1*c1;
	 BeamK(3,2)=B*s1;
	 BeamK(3,3)=EA/L*c1*c1+Am*s1*s1;
	 BeamK(3,4)=EA/L*s1*c1-Am*s1*c1;
	 BeamK(3,5)=C*s1;

	 BeamK(4,0)=-EA/L*c1*s1+Am*c1*s1;
	 BeamK(4,1)=-EA/L*s1*s1-Am*c1*c1;
	 BeamK(4,2)=-B*c1;
	 BeamK(4,3)=EA/L*c1*s1-Am*s1*c1;
	 BeamK(4,4)=EA/L*s1*s1+Am*c1*c1;
	 BeamK(4,5)=-C*c1;

	 BeamK(5,0)=-C*s1;
	 BeamK(5,1)=C*c1;
	 BeamK(5,2)=fb1;
	 BeamK(5,3)=C*s1;
	 BeamK(5,4)=-C*c1;
	 BeamK(5,5)=fd1;
    return BeamK;
}

const Matrix &
BeamGT::getInitialStiff(void)
{
 //Calculo la flexibilidad y luego la invierto
	 double L=trans(0,0);
     double c1=trans(0,1);
     double	s1=trans(0,2);
	 double ei= theMaterial[0]->getInitialTangent();
	 double EA = theMaterial3->getInitialTangent();
	 double GA = theMaterial2->getInitialTangent();
	 double f11 = L / (3.0 * ei) + FC * 1.2 / (L * GA);
	 double f22 = L / (3.0 * ei) + FC * 1.2 / (L * GA);
	 double f12 = -L / (6.0 * ei) + FC * 1.2 / (L * GA);


	 double det=f11*f22-f12*f12;
	 double fa1=f22/det;
	 double fb1=-f12/det;
	 double fc1=-f12/det;
	 double fd1=f11/det;
	 double x=1.0/L;

	 double Am=(fa1+2*fb1+fd1)*x*x;
	 double B=(fa1+fb1)*x;
	 double C=(fb1+fd1)*x;

	 //Beam tangent stiffness matrix in global coordinates
	 BeamK(0,0)=EA/L*c1*c1+Am*s1*s1;
     BeamK(0,1)=EA/L*c1*s1-Am*c1*s1;
	 BeamK(0,2)=-B*s1;
	 BeamK(0,3)=-EA/L*c1*c1-Am*s1*s1;
	 BeamK(0,4)=-EA/L*c1*s1+Am*s1*c1;
	 BeamK(0,5)=-C*s1;

	 BeamK(1,0)=BeamK(0,1);
	 BeamK(1,1)=EA/L*s1*s1+Am*c1*c1;
	 BeamK(1,2)=B*c1;
	 BeamK(1,3)=-EA/L*s1*c1+Am*s1*c1;
	 BeamK(1,4)=-EA/L*s1*s1-Am*c1*c1;
	 BeamK(1,5)=C*c1;

	 BeamK(2,0)=BeamK(0,2);
	 BeamK(2,1)=BeamK(1,2);
	 BeamK(2,2)=fa1;
	 BeamK(2,3)=B*s1;
	 BeamK(2,4)=-B*c1;
	 BeamK(2,5)=fb1;

	 BeamK(3,0)=-EA/L*c1*c1-Am*s1*s1;
	 BeamK(3,1)=-EA/L*c1*s1+Am*s1*c1;
	 BeamK(3,2)=B*s1;
	 BeamK(3,3)=EA/L*c1*c1+Am*s1*s1;
	 BeamK(3,4)=EA/L*s1*c1-Am*s1*c1;
	 BeamK(3,5)=C*s1;

	 BeamK(4,0)=-EA/L*c1*s1+Am*c1*s1;
	 BeamK(4,1)=-EA/L*s1*s1-Am*c1*c1;
	 BeamK(4,2)=-B*c1;
	 BeamK(4,3)=EA/L*c1*s1-Am*s1*c1;
	 BeamK(4,4)=EA/L*s1*s1+Am*c1*c1;
	 BeamK(4,5)=-C*c1;

	 BeamK(5,0)=-C*s1;
	 BeamK(5,1)=C*c1;
	 BeamK(5,2)=fb1;
	 BeamK(5,3)=C*s1;
	 BeamK(5,4)=-C*c1;
	 BeamK(5,5)=fd1;

	// 	 	opserr<<" initial stiffness "<< "\n";
    return BeamK;
}

const Vector &
BeamGT::getResistingForce()
{       

		const Vector &disp1 = theNodes[0]->getTrialDisp();
		const Vector &disp2 = theNodes[1]->getTrialDisp();
		double L=trans(0,0);
        double c1=trans(0,1);
        double s1=trans(0,2);

		double dx1=disp1(0);
		double dy1=disp1(1);
		double dr1=disp1(2);
		double dx2=disp2(0);
		double dy2=disp2(1);
		double dr2=disp2(2);
		double du1=dx1*c1+dy1*s1;
		double dv1=-dx1*s1+dy1*c1;
		double du2=dx2*c1+dy2*s1;
		double dv2=-dx2*s1+dy2*c1;
		

	
		double M1=theMaterial[0]->getStress();
		double M2=theMaterial[1]->getStress();
		double V1=theMaterial2->getStress();
		double N=-1.0*theMaterial3->getStress();
		double V2=-V1;
		double N2=-N;
		BeamR(0)=1.0*(N*c1-V1*s1);
		BeamR(1)=1.0*(N*s1+V1*c1);
		BeamR(2)=1.0*M1;
		BeamR(3)=(N2*c1-V2*s1);
		BeamR(4)=(N2*s1+V2*c1);
		BeamR(5)=1.0*M2;
      return BeamR;
}

int
BeamGT::sendSelf(int commitTag, Channel &theChannel)
{
    int res;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    
	int dataTag = this->getDbTag();

    // Truss2D packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments
	
	Vector data(17);
	data(0)=this->getTag();
	data(5)=LP1;
	data(6)=LP2;
	data(7)=LR;
	data(8)=theMaterial[0]->getClassTag();
	data(9)=theMaterial[1]->getClassTag();
	data(10)=theMaterial2->getClassTag();
    data(14)=theMaterial3->getClassTag();
	int matDbTag1=theMaterial[0]->getDbTag();
	int matDbTag2=theMaterial[1]->getDbTag();
	int matDbTag3=theMaterial2->getDbTag();
	int matDbTag4=theMaterial3->getDbTag();
	
	if (matDbTag1==0){
		matDbTag1=theChannel.getDbTag();
		 if (matDbTag1 !=0)
			 theMaterial[0]->setDbTag(matDbTag1);
	}
	data(11)=matDbTag1;

	if (matDbTag2==0){
		matDbTag2=theChannel.getDbTag();
		 if (matDbTag2 !=0)
			 theMaterial[1]->setDbTag(matDbTag2);
	}
	data(12)=matDbTag2;

		if (matDbTag3==0){
		matDbTag3=theChannel.getDbTag();
		 if (matDbTag3 !=0)
			 theMaterial2->setDbTag(matDbTag3);
	}
	data(13)=matDbTag3;



		if (matDbTag4==0){
		matDbTag4=theChannel.getDbTag();
		 if (matDbTag4 !=0)
			 theMaterial3->setDbTag(matDbTag4);
	}
	data(15)=matDbTag4;
	data(16) = FC;

  	res=0;
     res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send Vector\n";
      return -1;
    }         

    // Truss2D then sends the tags of it's two end nodes
     res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send ID\n";
      return -2;
    }

    // finally Truss2D asks it's material object to send itself
    res = theMaterial[0]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send the Material\n";
      return -3;
    }
    res = theMaterial[1]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send the Material\n";
      return -3;
    }
	    res = theMaterial2->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send the Material\n";
      return -3;
    }

	    res = theMaterial3->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BeamGT::sendSelf() - failed to send the Material\n";
      return -3;
    }

    return 0;
}

int
BeamGT::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{


    int res;
    int dataTag = this->getDbTag();

    // Truss2D creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector

    Vector data(17);
    res = theChannel.recvVector(dataTag, commitTag, data);
    res=-1;
	if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive Vector\n";
      return -1;
    }         

    this->setTag((int)data(0));
   

	LP1=data(5);
	LP2=data(6);
	LR=data(7);
	
	    
    // Truss2D now receives the tags of it's two external nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive ID\n";
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    int matClass1 = data(8);
	int matClass2 = data(9);
	int matClass3 = data(10);
    int matDb1 = data(11);
	int matDb2 = data(12);
	int matDb3 = data(13);
	int matClass4 = data(14);
	int matDb4 = data(15);
	FC = data(16);
    theMaterial[0] = theBroker.getNewUniaxialMaterial(matClass1);
    if (theMaterial[0] == 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to create a Material[0]\n";
      return -3;
    }


  theMaterial[1] = theBroker.getNewUniaxialMaterial(matClass2);
    if (theMaterial[1] == 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to create a Material[1]\n";
      return -3;
    }
	  theMaterial2 = theBroker.getNewUniaxialMaterial(matClass3);
    if (theMaterial2 == 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to create a Material2\n";
      return -3;
    }
	  theMaterial3 = theBroker.getNewUniaxialMaterial(matClass4);
    if (theMaterial3 == 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to create a Material3\n";
      return -3;
    }




    // we set the dbTag before we receive the material  - this is important
    theMaterial[0]->setDbTag(matDb1); 
    res = theMaterial[0]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive the Material\n";
      return -3;
    }
    theMaterial[1]->setDbTag(matDb2); 
    res = theMaterial[1]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive the Material\n";
      return -3;
    }

	    theMaterial2->setDbTag(matDb3); 
    res = theMaterial2->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive the Material\n";
      return -3;
    }


	    theMaterial3->setDbTag(matDb4); 
    res = theMaterial2->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING BeamGT::recvSelf() - failed to receive the Material\n";
      return -3;
    }

    return 0;
}

void
BeamGT::Print(OPS_Stream &s, int flag)
{

  s << " " << "\n";
  s << " " << "\n";
  s << "Element: " << this->getTag(); 
  s << " type: BeamGT " << "\n";
  s << " " << "\n";
  s << "+--------------------------------------------------------+"<< "\n";
  s << "|        Beam with Flexure and Shear Hinges              |\n";
  s << "|   Written by Gonzalo Torrisi UNCuyo Copyright 2016     |\n";
  s << "|                 Only in plane X-Y                      |\n";
  s << "|                Use at your Own Peril                   |\n";
  s << "+--------------------------------------------------------+"<<"\n";
  s << "             Nodes: " << "\n";
  s << "Nodo 1  :"<< externalNodes(0)<< "\n";
  s << "Nodo 2  :"<< externalNodes(1)<< "\n";
  s << "         BeamGT Materials: " << "\n";
  s << "Material for Flexure 1 :" << *theMaterial[0]<< "\n";
  s << "Material for Flexure 2 :" << *theMaterial[1]<< "\n";
  s << "Material for Shear     :" << *theMaterial2<< "\n";
  s << "Material for Axial     :" << *theMaterial3<< "\n";
  s << "Nonlinear properties:" << "\n";
  s << "Plastic hinge length end 1 :" << LP1 << "\n";
  s << "Plastic hinge length end 2 :" << LP2 << "\n";
  s << "Length subjected to shear  :" << LR << "\n";
  s << "Factor for shear correction:" << FC << "\n";
  s << " " << "\n";
}

Response*
BeamGT::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;
//			opserr<<"En setResponse"<<"\n";

    output.tag("ElementOutput");
    output.attr("eleType","BeamGT");
    output.attr("eleTag",this->getTag());
    output.attr("node1 ",externalNodes[0]);
    output.attr("node2 ",externalNodes[1]);

    char outputData[10];

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            char outputData[10];
//           int numDOFperNode = numDOF/2;
            for (int i=0; i<4; i++) {
               sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
           for (int j=0; j<4; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(3));



    } else if ((strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) ||
	       (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)) {

       for (int i=0; i<6; i++) {
           sprintf(outputData,"P%d",i+1);
           output.tag("ResponseType",outputData);


       }
       theResponse = new ElementResponse(this, 2, Vector(6));


    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	       strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"basicDeformation") == 0) {

           for (int i=0; i<6; i++) {
               sprintf(outputData,"e%d",i+1);
               output.tag("ResponseType",outputData);
           }
           theResponse = new ElementResponse(this, 3, Vector(4));
	



    } else if (strcmp(argv[0],"basicStiffness") == 0) {

            for (int i=0; i<6; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 13, Matrix(6,6));

    } else if ((strcmp(argv[0],"defoANDforce") == 0) ||
        (strcmp(argv[0],"deformationANDforces") == 0) ||
        (strcmp(argv[0],"deformationsANDforces") == 0)) {

            int i;
            for (i=0; i<4; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            for (i=0; i<4; i++) {
                sprintf(outputData,"P%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 4, Vector(2*4));

    // a material quantity
    } else if (strcmp(argv[0],"material") == 0) {
      if (argc > 2) {
	int matNum = atoi(argv[1]);
	if (matNum >= 1 && matNum <= 2)
	  theResponse =  theMaterial[matNum-1]->setResponse(&argv[2], argc-2, output);
      }
    	  theResponse =  theMaterial2->setResponse(&argv[2], argc-2, output);
	}


    output.endTag();

    return theResponse;
}

int 
BeamGT::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
	
    const Vector  diff  = disp2-disp1;
//				opserr<<"En getResponse"<<"\n";
    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    case 2:
        if (eleInformation.theVector != 0) {
       //     for (int i = 0; i < 2; i++)
       //         (*(eleInformation.theVector))(i) = theMaterial[i]->getStress();
       //            (*(eleInformation.theVector))(2) = theMaterial2->getStress();       		
       //           (*(eleInformation.theVector))(3) = theMaterial3->getStress();     

				  (*(eleInformation.theVector))(0) = -1.0*theMaterial3->getStress();
				  (*(eleInformation.theVector))(1) = theMaterial2->getStress();
				  (*(eleInformation.theVector))(2) = theMaterial[0]->getStress();
				  (*(eleInformation.theVector))(3) = theMaterial3->getStress();
				  (*(eleInformation.theVector))(4) = -1.0*theMaterial2->getStress();
				  (*(eleInformation.theVector))(5) = theMaterial[1]->getStress();

         }
        return 0;

    case 3:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < 2; i++)
             (*(eleInformation.theVector))(i) = theMaterial[i]->getStrain();
			 (*(eleInformation.theVector))(2) = theMaterial2->getStrain();
			 (*(eleInformation.theVector))(3) = theMaterial3->getStrain();
		}
        return 0;

    case 13:
        if (eleInformation.theMatrix != 0) {
            for (int i = 0; i < 2; i++){
	      (*(eleInformation.theMatrix))(i,i) = theMaterial[i]->getTangent();
		}
			(*(eleInformation.theMatrix))(2,2) = theMaterial2->getTangent();
			(*(eleInformation.theMatrix))(3,3) = theMaterial3->getTangent();
		}
        return 0;

    case 4:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < 2; i++) {
                (*(eleInformation.theVector))(i) = theMaterial[i]->getStrain();
                (*(eleInformation.theVector))(i+4) = theMaterial[i]->getStress();    
			}
			     (*(eleInformation.theVector))(2) = theMaterial2->getStrain();
			     (*(eleInformation.theVector))(6) = theMaterial2->getStress();   

			     (*(eleInformation.theVector))(3) = theMaterial3->getStrain();
			     (*(eleInformation.theVector))(7) = theMaterial3->getStress();    

		}
        return 0;    

    default:
        return -1;
    }
}


double
BeamGT::computeCurrentStrain(int mat) const
//BeamGT::computeCurrentStrain(void) const
{
			// NOTE this method will never be called with L == 0.0
		double str[6];
		double strain;
		// determine the strain
		const Vector &disp1 = theNodes[0]->getTrialDisp();
		const Vector &disp2 = theNodes[1]->getTrialDisp();  
     
	    
        double L=trans(0,0);
        double c1=trans(0,1);
        double s1=trans(0,2);
   
	//deformaciones en locales del elemento.
		double dx1=disp1(0);
		double dy1=disp1(1);
		double dr1=disp1(2);
		double dx2=disp2(0);
		double dy2=disp2(1);
		double dr2=disp2(2);

		double du1=dx1*c1+dy1*s1;
		double dv1=-dx1*s1+dy1*c1;
		double du2=dx2*c1+dy2*s1;
		double dv2=-dx2*s1+dy2*s1;

		double ro1=dr1+(dv2-dv1)/L;
		double ro2=dr2+(dv2-dv1)/L;
		double dg=(dv2-dv1)/L;

		str[0]=du1;
		str[1]=dv1;
		str[2]=dr1;
		str[3]=du2;
		str[4]=dv2;
		str[5]=dr2;
		strain=str[mat];
		return strain;
}

int
BeamGT::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
//BeamGT::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	          int code=0;
	// first determine the two end points of the CorotTruss2 based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	

//						
	const Vector &end1Disp = theNodes[0]->getDisp();
	const Vector &end2Disp = theNodes[1]->getDisp();   
    
	static Vector v1(3);
	static Vector v2(3);
	static Vector rgb(3);
	static Vector v3(3);
	static Vector v4(3);


	theNodes[0]->getDisplayCrds(v3, fact, displayMode);
	theNodes[1]->getDisplayCrds(v4, fact, displayMode);

	for (int i = 0; i < 2; i++) {
		v1(i) = end1Crd(i)+end1Disp(i)*fact;
		v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}

	// compute the strain and axial force in the member
	double strain[6], force[6];
		    
        double L=trans(0,0);
        double c1=trans(0,1);
        double	s1=trans(0,2);
   
	//deformaciones en locales del elemento.

	for (int i = 0; i < 1; i++) {

	    strain[i] = this->computeCurrentStrain(i);
	}
		double ro1=strain[2]+(strain[4]-strain[1])/L;
		double ro2=strain[5]+(strain[4]-strain[1])/L;
		double dg=(strain[4]-strain[1])/L;

		double EA=theMaterial3->getInitialTangent();
		force[0]=EA/L*(strain[3]-strain[0]);
		force[1]=theMaterial2->getStress();
		force[2]=theMaterial[0]->getStress();
		force[3]=force[0];
		force[4]=force[1];
		force[5]=theMaterial[1]->getStress();

		//		
	if (displayMode == 1) { //  use the moment as measure
	code = 0;
	code += theViewer.drawLine(v1, v2, (float)force[2], (float)force[5]);

	return code;
	}

	else if (displayMode == 2) // use the strain as the drawing measure
		{
			code=0;
		  code += theViewer.drawLine(v1, v2, (float)strain[2], (float)strain[5]);
		  return code;
}

	else if (displayMode == 3) // use the strain as the drawing measure
	{
		code = 0;
		double ww = 10*(force[2] + force[5]) /L;
		//code += theViewer.drawLine(v1, v2, (float)force[2], (float)force[5],0,0,ww,1);
		code += theViewer.drawLine(v1, v2, (float)force[2], (float)force[5]);
		return code;
	}

	else if(displayMode < 0)
	{
		code = 0;
		code += theViewer.drawLine(v3, v4, 1.0, 1.0, this->getTag(), 0);
		return code;
	}


	else { 
		code=0;
		  code += theViewer.drawLine(v1, v2, 0, 0);

		  return code;
	}



	return 0;
}

 
