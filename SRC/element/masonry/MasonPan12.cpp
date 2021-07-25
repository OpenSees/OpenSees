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

// we specify what header files we need
#include "MasonPan12.h"
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
#include <Parameter.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix MasonPan12::PanelK(36,36);
Vector MasonPan12::PanelR(36);

void* OPS_MasonPan12()

//#define OPS_Export 

//OPS_Export void *
//OPS_MasonPan12element()
{
  // print out a message about who wrote this element & any copyright info wanted
  //if (numMyPanel == 0) {
  // opserr << " \n";
  //  opserr << "                 REFINED MASONRY PANEL\n";
  //  opserr << "   Written by Gonzalo Torrisi UNCuyo Copyright 2015\n";
  //	opserr << "         Model with 6 compression struts-36 dof\n";
  //	opserr << "                 Only in plane X-Y \n";
   // opserr << "                Use at your Own Peril\n";

//    numMyPanel++;
//  }

   Element *thePanel = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    thePanel = new MasonPan12();
    return thePanel;
  }

  if (numRemainingArgs != 18) {
    opserr << "ERROR - Masonry Panel not enough args provided, want: element MasonryPanel tag? Node1? Node2? Node3? Node4?  Node5?  Node6?  Node7?  Node8?  Node9?   Node10?   Node11?   Node12?   matTag? matTag2? thick? wfactor? w1?\n";
//    numMyPanel++;
  }

  // get the id and end nodes 
  int iData[15];
  double dData[3];
  int numData;

  numData = 13;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }
    int eleTag = iData[0];
 
  numData =1;
  if (OPS_GetIntInput(&numData, &iData[13]) != 0) {
    opserr << "WARNING error reading element material 1 tag for element " << eleTag << endln;
    return 0;
  }
      numData =1;
  if (OPS_GetIntInput(&numData, &iData[14]) != 0) {
    opserr << "WARNING error reading element material 2 tag for element " << eleTag << endln;
    return 0;
  }

  int matID = iData[13];
  int matID2 = iData[14];

  numData =3 ;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element areas, thickness and properties for element" << eleTag << endln;
    return 0;
  }

      UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
	  UniaxialMaterial *theMaterial2 =OPS_GetUniaxialMaterial(matID2);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element " << eleTag << endln;
    return 0;
  }

    // now create the truss and add it to the Domain

  thePanel= new MasonPan12(eleTag, iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8], iData[9], iData[10], iData[11], iData[12], *theMaterial, *theMaterial2, dData[0], dData[1], dData[2]);

  if (thePanel == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    delete theMaterial;
	delete theMaterial2;
    return 0;
  }

  return thePanel;
}


// typical constructor
MasonPan12::MasonPan12(int tag, 
                   int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6,
                   int Nd7, int Nd8, int Nd9, int Nd10, int Nd11, int Nd12,
                   UniaxialMaterial &theMat, UniaxialMaterial &theMat2,
                    double thick, double wr,double w1)
:Element(tag, ELE_TAG_MasonPan12),     
 externalNodes(12),
 trans(8,4),rig1(6),rig2(6),rig3(6), TH(thick), W1(w1), WR(wr), theMaterial(0),theMaterial2(0)
 //L1(0.0),L2(0.0),L3(0.0),L4(0.0),L5(0.0),L6(0.0)
{       
  // get a copy of the material object for our own use
 // theMaterial = theMat.getCopy();
  //if (theMaterial == 0) {
  // opserr << "FATAL MassonPan::MasonPan12() - out of memory, could not get a copy of the Material\n";
  // exit(-1);
 // }

    // allocate memory for numMaterials1d uniaxial material models
  theMaterial = new UniaxialMaterial *[6];
  theMaterial2= 0;


  if ( theMaterial == 0 ) {
    opserr << "FATAL MasonPan12::MasonPan12 - failed to create a 1d  material or direction array\n";
    exit(-1);
  }

  // get a copy of the material and check we obtained a valid copy
 //for (int i=0; i<6; i++) {
 
      theMaterial[0] = theMat.getCopy();
	  theMaterial[3] = theMat.getCopy();
      
	  if (theMaterial[0] == 0) {
        opserr << "FATAL MasonPan12::MasonPan12 - failed to get a copy of material\n" ;
        exit(-1);
      }
	  	  if (theMaterial[3] == 0) {
        opserr << "FATAL MasonPan12::MasonPan12 - failed to get a copy of material\n" ;
        exit(-1);
      }
  
     theMaterial[1] = theMat2.getCopy();
     theMaterial[2] = theMat2.getCopy();
     theMaterial[4] = theMat2.getCopy();
     theMaterial[5] = theMat2.getCopy();

 // 	 if (theMaterial2== 0){
//		 opserr << "FATAL MasonPan12::MasonPan12 - failed to get a copy of material2\n";
//	 }

  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 12) {
    opserr << "FATAL MassonPan::MasonPan12() - out of memory, could not create an ID of size 12\n";
    exit(-1);
  }

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;   
  externalNodes(2) = Nd3;   
  externalNodes(3) = Nd4;   
  externalNodes(4) = Nd5;   
  externalNodes(5) = Nd6;   
  externalNodes(6) = Nd7;   
  externalNodes(7) = Nd8;   
  externalNodes(8) = Nd9; 
  externalNodes(9) = Nd10; 
  externalNodes(10) = Nd11; 
  externalNodes(11) = Nd12; 


  theNodes[0] = 0; 
  theNodes[1] = 0;
  theNodes[2] = 0; 
  theNodes[3] = 0;
  theNodes[4] = 0; 
  theNodes[5] = 0;
  theNodes[6] = 0; 
  theNodes[7] = 0;
  theNodes[8] = 0; 
  theNodes[9] = 0;
  theNodes[10] = 0;
  theNodes[11] = 0;
}

// constructor which should be invoked by an FE_ObjectBroker only
MasonPan12::MasonPan12()
:Element(0, ELE_TAG_MasonPan12),     
 theMaterial(0), theMaterial2(0),
 externalNodes(12),
 trans(8,4),rig1(6),rig2(6),rig3(6), TH(0.0),W1(0.0),WR(0.0)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
  theNodes[2] = 0; 
  theNodes[3] = 0;
  theNodes[4] = 0; 
  theNodes[5] = 0;
  theNodes[6] = 0; 
  theNodes[7] = 0;
  theNodes[8] = 0; 
  theNodes[9] = 0;
  theNodes[10] = 0;
  theNodes[11] = 0;
}

//  destructor - provided to clean up any memory
MasonPan12::~MasonPan12()
{
    // clean up the memory associated with the element, this is
    // memory the Truss2D objects allocates and memory allocated 
    // by other objects that the Truss2D object is responsible for 
    // cleaning up, i.e. the MaterialObject.
	for (int i=0; i<6; i++) {
    if (theMaterial[i] != 0)
        delete theMaterial[i];    
}
	delete [] theMaterial;
	delete theMaterial2;
}

int
MasonPan12::getNumExternalNodes(void) const
{
    return 12;
}

const ID &
MasonPan12::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
MasonPan12::getNodePtrs(void) 
{
  return theNodes;
}

int
MasonPan12::getNumDOF(void) {
    return 36;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.

void
MasonPan12::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        return;
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    Node *end1Ptr, *end2Ptr, *end3Ptr, *end4Ptr, *end5Ptr, *end6Ptr, *end7Ptr, *end8Ptr, *end9Ptr, *end10Ptr, *end11Ptr, *end12Ptr;
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
        int Nd3 = externalNodes(2);
        int Nd4 = externalNodes(3);
        int Nd5 = externalNodes(4);
        int Nd6 = externalNodes(5);
        int Nd7 = externalNodes(6);
        int Nd8 = externalNodes(7);
        int Nd9 = externalNodes(8);
        int Nd10= externalNodes(9);
        int Nd11=externalNodes(10);
        int Nd12= externalNodes(11);

    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);  
    end3Ptr = theDomain->getNode(Nd3);  
    end4Ptr = theDomain->getNode(Nd4);  
    end5Ptr = theDomain->getNode(Nd5);  
    end6Ptr = theDomain->getNode(Nd6);  
    end7Ptr = theDomain->getNode(Nd7);  
    end8Ptr = theDomain->getNode(Nd8);  
    end9Ptr = theDomain->getNode(Nd9);  
    end10Ptr = theDomain->getNode(Nd10);        
    end11Ptr = theDomain->getNode(Nd11);        
    end12Ptr = theDomain->getNode(Nd12);        

    if (end1Ptr == 0) {
      opserr << "WARNING MasonPan12::setDomain() - at truss " << this->getTag() << " node " <<
        Nd1 << "  does not exist in domain\n";
                                
        return;  // don't go any further - otherwise segemntation fault
    }
    if (end12Ptr == 0) {        
      opserr << "WARNING MasonPan12::setDomain() - at truss " << this->getTag() << " node " <<
        Nd2 << "  does not exist in domain\n";

        return;  // don't go any further - otherwise segemntation fault
    }   

        theNodes[0] = end1Ptr;
		theNodes[1] = end2Ptr;
		theNodes[2] = end3Ptr;
		theNodes[3] = end4Ptr;
		theNodes[4] = end5Ptr;
		theNodes[5] = end6Ptr;
        theNodes[6] = end7Ptr;
		theNodes[7] = end8Ptr;
		theNodes[8] = end9Ptr;
		theNodes[9] = end10Ptr;
        theNodes[10] = end11Ptr;
		theNodes[11] = end12Ptr;
        
        
        // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();       
    int dofNd3 = end3Ptr->getNumberDOF();
    int dofNd4 = end4Ptr->getNumberDOF();       
    int dofNd5 = end5Ptr->getNumberDOF();
    int dofNd6 = end6Ptr->getNumberDOF();       
    int dofNd7 = end7Ptr->getNumberDOF();
    int dofNd8 = end8Ptr->getNumberDOF();       
    int dofNd9 = end9Ptr->getNumberDOF();
    int dofNd10 = end10Ptr->getNumberDOF();     
    int dofNd11 = end11Ptr->getNumberDOF();
    int dofNd12 = end12Ptr->getNumberDOF();     

    if ((dofNd1 != 3 || dofNd2 != 3)) {
      opserr << "MasonPan12::setDomain(): 2 dof required at nodes\n";
      return;
    }   

    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds(); 
    const Vector &end3Crd = end3Ptr->getCrds();
    const Vector &end4Crd = end4Ptr->getCrds(); 
    const Vector &end5Crd = end5Ptr->getCrds();
    const Vector &end6Crd = end6Ptr->getCrds(); 
    const Vector &end7Crd = end7Ptr->getCrds();
    const Vector &end8Crd = end8Ptr->getCrds(); 
    const Vector &end9Crd = end9Ptr->getCrds();
    const Vector &end10Crd = end10Ptr->getCrds();       
    const Vector &end11Crd = end11Ptr->getCrds();
    const Vector &end12Crd = end12Ptr->getCrds();       

    double dx = end4Crd(0)-end1Crd(0);
    double dy = end10Crd(1)-end1Crd(1); 
         
    if (dx == 0.0) {
      opserr << "WARNING MasonPan12::setDomain() - MasonPan12 " << this->getTag() << 
        " has zero length\n";
      return;  // don't go any further - otherwise divide by 0 error
    }
    if (dy == 0.0) {
      opserr << "WARNING MasonPan12::setDomain() - MasonPan12 " << this->getTag() << 
        " has zero height\n";
      return;  // don't go any further - otherwise divide by 0 error
    }

        // cosenos de las diagonales
        double dx1 =end4Crd(0)-end10Crd(0);
        double dy1= end4Crd(1)-end10Crd(1);
		double L1=sqrt(dx1*dx1+dy1*dy1);
        double dx2 =end3Crd(0)-end11Crd(0);
        double dy2= end3Crd(1)-end11Crd(1);
		double L2=sqrt(dx2*dx2+dy2*dy2);
        double dx3 =end5Crd(0)-end9Crd(0);
        double dy3= end5Crd(1)-end9Crd(1);
		double L3=sqrt(dx3*dx3+dy3*dy3);
        
        double dx4 =end7Crd(0)-end1Crd(0);
		double dy4= end7Crd(1)-end1Crd(1);
		double L4=sqrt(dx4*dx4+dy4*dy4);
        double dx5 =end6Crd(0)-end2Crd(0);
        double dy5= end6Crd(1)-end2Crd(1);
		double L5=sqrt(dx5*dx5+dy5*dy5);
        double dx6 =end8Crd(0)-end12Crd(0);
        double dy6= end8Crd(1)-end12Crd(1);
		double L6=sqrt(dx6*dx6+dy6*dy6);
		double Area1=L1*WR*TH*W1;
		double Area2=L1*WR*TH*(1-W1)/2;
		double Area3=L1*WR*TH*(1-W1)/2;
		double Area4=L1*WR*TH*W1;
		double Area5=L1*WR*TH*(1-W1)/2;
		double Area6=L1*WR*TH*(1-W1)/2;
		double Lpan=end4Crd(0)-end1Crd(0);
		double Apan=Lpan*TH;

		//rigidez del resorte
		 //suma de cosenos l cuadrado
		// double coef=(dx1/L1*dx1/L1)/(L1)+(dx2/L2*dx2/L2)/(L2)+(dx3/L3*dx3/L3)/(L3)+(dx4/L4*dx4/L4)/(L4)+(dx5/L5*dx5/L5)/(L5)+(dx6/L6*dx6/L6)/(L6);
		 // rigidez de las bielas
		 double Emam=0.0;
		  Emam= theMaterial[0]->getInitialTangent();
		 //double Gcor=theMaterial2->getInitialTangent();

		//double coef2=Emam*(Area1+Area2+Area3)*coef;
		//rigidez del resorte
	//	double Rigspring0=WS*coef2;
		//modulo de corte equivalente
	//	double Gspr=Rigspring0*fabs(dy1)/Apan;
		//rigidez relativa
		//double Rigspring=fabs(dy1)/Apan;
	//	double Rigspring=Rigspring0/Emam;

	//	double coef3=1.0-WS;
		trans(0,0)=L1;
		trans(0,1)=dx1/L1;
		trans(0,2)=dy1/L1;
		trans(0,3)=Area1;
		trans(1,0)=L2;
		trans(1,1)=dx2/L2;
		trans(1,2)=dy2/L2;
		trans(1,3)=Area2;
		trans(2,0)=L3;
		trans(2,1)=dx3/L3;
		trans(2,2)=dy3/L3;
		trans(2,3)=Area3;
		trans(3,0)=L4;
		trans(3,1)=dx4/L4;
		trans(3,2)=dy4/L4;
		trans(3,3)=Area4;
		trans(4,0)=L5;
		trans(4,1)=dx5/L5;
		trans(4,2)=dy5/L5;
		trans(4,3)=Area5;
		trans(5,0)=L6;
		trans(5,1)=dx6/L6;
		trans(5,2)=dy6/L6;
		trans(5,3)=Area6;
		trans(6,0)=Apan;
		trans(6,1)=0.0;
		trans(6,2)=0.0;
		trans(6,3)=0.0;
		trans(7,0)=dy;
		for (int im=0; im<6; im++) {
			rig1(im) = trans(im,1)*trans(im,1)*trans(im,3)/trans(im,0);
			rig2(im) = trans(im,1)*trans(im,2)*trans(im,3)/trans(im,0);
		    rig3(im) = trans(im,2)*trans(im,2)*trans(im,3)/trans(im,0);
		}

}        

int
MasonPan12::commitState()
{
	 int ecode=0;
	 //commit material models
	 for (int i=0; i<6; i++)
		 ecode += theMaterial[i]->commitState();
	 //commit the base class
	//     ecode += theMaterial2->commitState();
	 
	 ecode += this->Element::commitState();
	  Cdeltares=Tdeltares;

	 return ecode;
}

int
MasonPan12::revertToLastCommit()
{
//    return theMaterial->revertToLastCommit();
    int code=0;
    
    // revert state for 1d materials
    for (int i=0; i<6; i++)
        code += theMaterial[i]->revertToLastCommit();
  //  code += theMaterial2->revertToLastCommit();
	Tdeltares=Cdeltares;
    return code;
}

int
MasonPan12::revertToStart()
{
  //  return theMaterial->revertToStart();
            int code=0;
    
    // revert state for 1d materials
    for (int i=0; i<6; i++)
        code += theMaterial[i]->revertToStart();
   // code += theMaterial2->revertToStart();

		return code;
}

int
MasonPan12::update()
{
  // determine the current strain given trial displacements at nodes
 // double strain = this->computeCurrentStrain();

  // set the strain in the materials
  //theMaterial->setTrialStrain(strain);

        // compute strain and rate; set as current trial for material
 //////////////////////////////////////       strain     = this->computeCurrentStrain(mat );
   // strainRate = this->computeCurrentStrain();
			double str[7];

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();  
    const Vector &disp3 = theNodes[2]->getTrialDisp();  
    const Vector &disp4 = theNodes[3]->getTrialDisp();  
    const Vector &disp5 = theNodes[4]->getTrialDisp();  
    const Vector &disp6 = theNodes[5]->getTrialDisp();  
    const Vector &disp7 = theNodes[6]->getTrialDisp();  
    const Vector &disp8 = theNodes[7]->getTrialDisp();  
    const Vector &disp9 = theNodes[8]->getTrialDisp();  
    const Vector &disp10 = theNodes[9]->getTrialDisp(); 
    const Vector &disp11 = theNodes[10]->getTrialDisp();        
    const Vector &disp12 = theNodes[11]->getTrialDisp();        


        double c1=trans(0,1);
        double s1=trans(0,2);
      double	L1=trans(0,0);
	  double A1=trans(0,3);
        double du1=((disp4(0)-disp10(0))*c1+(disp4(1)-disp10(1))*s1)/L1;
        c1=trans(1,1);
       double s2=trans(1,2);
		double L2=trans(1,0);
		 double A2=trans(1,3);
        double du2=((disp3(0)-disp11(0))*c1+(disp3(1)-disp11(1))*s2)/L2;
        c1=trans(2,1);
      double  s3=trans(2,2);
		 double L3=trans(2,0);
		  double A3=trans(2,3);
        double du3=((disp5(0)-disp9(0))*c1+(disp5(1)-disp9(1))*s3)/L3;
        c1=trans(3,1);
     double   s4=trans(3,2);
		double L4=trans(3,0);
		 double A4=trans(3,3);
        double du4=((disp7(0)-disp1(0))*c1+(disp7(1)-disp1(1))*s4)/L4;
        c1=trans(4,1);
     double   s5=trans(4,2);
	double	L5=trans(4,0);
	 double A5=trans(4,3);
        double du5=((disp6(0)-disp2(0))*c1+(disp6(1)-disp2(1))*s5)/L5;
        c1=trans(5,1);
    double    s6=trans(5,2);
	double	L6=trans(5,0);
	 double A6=trans(5,3);
        double du6=((disp8(0)-disp12(0))*c1+(disp8(1)-disp12(1))*s6)/L6;
   //resorte
		
		double deltares=((disp7(0)+disp10(0))-(disp1(0)+disp4(0)))/trans(7,0);

          str[0]=du1;
		  str[1]=du2;
		  str[2]=du3;
		  str[3]=du4;
		  str[4]=du5;
		  str[5]=du6;
		  str[6]=deltares;
		
		  int idir;
		if (deltares>=0) {
		   idir=1;
		  }
		  else
		  {
		   idir=2;
		 }
		  Tdeltares=deltares;
	
		  int ret = 0;
	for (int mat=0; mat<6; mat++) {
                 
				 double strain=0;
				 strain=str[mat];
        ret+= theMaterial[mat]->setTrialStrain(strain);
    }

	//double Axial = 0.0;
	//double FAxial = 0.0;
	//   double force1 = A1*theMaterial[0]->getStress();  
	//   double force2 = A2*theMaterial[1]->getStress(); 
	//   double force3 = A3*theMaterial[2]->getStress(); 
	//   double force4 = A4*theMaterial[3]->getStress(); 
	//   double force5 = A5*theMaterial[4]->getStress(); 
	//   double force6 = A6*theMaterial[5]->getStress(); 
	//   FAxial=force1*abs(s1)+force2*abs(s2)+force3*abs(s3)+force4*abs(s4)+force5*abs(s5)+force6*abs(s6);
    //   Axial=FAxial/trans(6,0);
	 //  Axial=FAxial/(TH*WR*L1);
	//   if (WS == 0.0)
	//	   Axial=0.0;

	 // double FFAXIAL=FAxial;
	 //  ret+=theMaterial2->setTrialStrain(deltares,Axial);

    return ret;


//  return 0;
}

const Matrix &
MasonPan12::getTangentStiff(void)
{

     double Et;

//diagonal 1
	     Et=theMaterial[0]->getTangent();
       	PanelK(9,9)=rig1(0)*Et;
        PanelK(9,10)=rig2(0)*Et;
        PanelK(10,9)=rig2(0)*Et;
        PanelK(10,10)=rig3(0)*Et;
        PanelK(27,27)=rig1(0)*Et;
        PanelK(27,28)=rig2(0)*Et;
        PanelK(28,27)=rig2(0)*Et;
        PanelK(28,28)=rig3(0)*Et;
        PanelK(9,27)=-rig1(0)*Et;
        PanelK(9,28)=-rig2(0)*Et;
        PanelK(10,27)=-rig2(0)*Et;
        PanelK(10,28)=-rig3(0)*Et;
        PanelK(27,9)=-rig1(0)*Et;
        PanelK(27,10)=-rig2(0)*Et;
        PanelK(28,9)=-rig2(0)*Et;
        PanelK(28,10)=-rig3(0)*Et;
   
//diagonal 2
		 Et=theMaterial[1]->getTangent();
        PanelK(6,6)=rig1(1)*Et;
        PanelK(6,7)=rig2(1)*Et;
        PanelK(7,6)=rig2(1)*Et;
        PanelK(7,7)=rig3(1)*Et;
        PanelK(30,30)=rig1(1)*Et;
        PanelK(30,31)=rig2(1)*Et;
        PanelK(31,30)=rig2(1)*Et;
        PanelK(31,31)=rig3(1)*Et;
        PanelK(6,30)=-rig1(1)*Et;
        PanelK(6,31)=-rig2(1)*Et;
        PanelK(7,30)=-rig2(1)*Et;
        PanelK(7,31)=-rig3(1)*Et;
        PanelK(30,6)=-rig1(1)*Et;
        PanelK(30,7)=-rig2(1)*Et;
        PanelK(31,6)=-rig2(1)*Et;
        PanelK(31,7)=-rig3(1)*Et;
    
//diagonal 3
		 Et=theMaterial[2]->getTangent();
        PanelK(12,12)=rig1(2)*Et;
        PanelK(12,13)=rig2(2)*Et;
        PanelK(13,12)=rig2(2)*Et;
        PanelK(13,13)=rig3(2)*Et;
        PanelK(24,24)=rig1(2)*Et;
        PanelK(24,25)=rig2(2)*Et;
        PanelK(25,24)=rig2(2)*Et;
        PanelK(25,25)=rig3(2)*Et;
        PanelK(12,24)=-rig1(2)*Et;
        PanelK(12,25)=-rig2(2)*Et;
        PanelK(13,24)=-rig2(2)*Et;
        PanelK(13,25)=-rig3(2)*Et;
        PanelK(24,12)=-rig1(2)*Et;
        PanelK(24,13)=-rig2(2)*Et;
        PanelK(25,12)=-rig2(2)*Et;
        PanelK(25,13)=-rig3(2)*Et;

//diagonal 4
		 Et=theMaterial[3]->getTangent();
        PanelK(0,0)=rig1(3)*Et;
        PanelK(0,1)=rig2(3)*Et;
        PanelK(1,0)=rig2(3)*Et;
        PanelK(1,1)=rig3(3)*Et;
        PanelK(18,18)=rig1(3)*Et;
        PanelK(18,19)=rig2(3)*Et;
        PanelK(19,18)=rig2(3)*Et;
        PanelK(19,19)=rig3(3)*Et;
        PanelK(0,18)=-rig1(3)*Et;
        PanelK(0,19)=-rig2(3)*Et;
        PanelK(1,18)=-rig2(3)*Et;
        PanelK(1,19)=-rig3(3)*Et;
        PanelK(18,0)=-rig1(3)*Et;
        PanelK(18,1)=-rig2(3)*Et;
        PanelK(19,0)=-rig2(3)*Et;
        PanelK(19,1)=-rig3(3)*Et;
    
//diagonal 5
		 Et=theMaterial[4]->getTangent();
        PanelK(3,3)=rig1(4)*Et;
        PanelK(3,4)=rig2(4)*Et;
        PanelK(4,3)=rig2(4)*Et;
        PanelK(4,4)=rig3(4)*Et;
        PanelK(15,15)=rig1(4)*Et;
        PanelK(15,16)=rig2(4)*Et;
        PanelK(16,15)=rig2(4)*Et;
        PanelK(16,16)=rig3(4)*Et;
        PanelK(3,15)=-rig1(4)*Et;
        PanelK(3,16)=-rig2(4)*Et;
        PanelK(4,15)=-rig2(4)*Et;
        PanelK(4,16)=-rig3(4)*Et;
        PanelK(15,3)=-rig1(4)*Et;
        PanelK(15,4)=-rig2(4)*Et;
        PanelK(16,3)=-rig2(4)*Et;
        PanelK(16,4)=-rig3(4)*Et;

//diagonal 6
		Et=theMaterial[5]->getTangent();
        PanelK(21,21)=rig1(5)*Et;
        PanelK(21,22)=rig2(5)*Et;
        PanelK(22,21)=rig2(5)*Et;
        PanelK(22,22)=rig3(5)*Et;
        PanelK(33,33)=rig1(5)*Et;
        PanelK(33,34)=rig2(5)*Et;
        PanelK(34,33)=rig2(5)*Et;
        PanelK(34,34)=rig3(5)*Et;
        PanelK(21,33)=-rig1(5)*Et;
        PanelK(21,34)=-rig2(5)*Et;
        PanelK(22,33)=-rig2(5)*Et;
        PanelK(22,34)=-rig3(5)*Et;
        PanelK(33,21)=-rig1(5)*Et;
        PanelK(33,22)=-rig2(5)*Et;
        PanelK(34,21)=-rig2(5)*Et;
        PanelK(34,22)=-rig3(5)*Et;
    
    
   //     double Et2= theMaterial2->getInitialTangent();
  //  	 Et= theMaterial2->getTangent();
//		double Rigsp=trans(6,3)*Et;
//		if (WS==0)
//			Rigsp=0;

//	const Vector &disp1 = theNodes[0]->getTrialDisp();
//    const Vector &disp4 = theNodes[3]->getTrialDisp();  
//    const Vector &disp7 = theNodes[6]->getTrialDisp();  
//    const Vector &disp10 = theNodes[9]->getTrialDisp(); 

//	double deltares=((disp7(0)+disp10(0))-(disp1(0)+disp4(0)))/trans(7,0);
	
//	int idir;
//		if (deltares-Cdeltares>=0) {
//		   idir=1;
//		  }
//		  else
//		  {
//		   idir=2;
//		 }
//					
//		if (idir==1)
//		{
//		
//		PanelK(27,27)+=Rigsp;
//		PanelK(9,9)+=Rigsp;
//		PanelK(27,9)-=Rigsp;
//		PanelK(9,27)-=Rigsp;
//		}
//		else
//		{
//		PanelK(0,0)+=Rigsp;
//		PanelK(18,18)+=Rigsp;
//		PanelK(0,18)-=Rigsp;
//		PanelK(18,0)-=Rigsp;
 //         }		
    // return the matrix
    return PanelK;
}

const Matrix &
MasonPan12::getInitialStiff(void)
{
	double E;
   E = theMaterial[0]->getInitialTangent();
//diagonal 1
       	PanelK(9,9)=rig1(0)*E;
        PanelK(9,10)=rig2(0)*E;
        PanelK(10,9)=rig2(0)*E;
        PanelK(10,10)=rig3(0)*E;
        PanelK(27,27)=rig1(0)*E;
        PanelK(27,28)=rig2(0)*E;
        PanelK(28,27)=rig2(0)*E;
        PanelK(28,28)=rig3(0)*E;
        PanelK(9,27)=-rig1(0)*E;
        PanelK(9,28)=-rig2(0)*E;
        PanelK(10,27)=-rig2(0)*E;
        PanelK(10,28)=-rig3(0)*E;
        PanelK(27,9)=-rig1(0)*E;
        PanelK(27,10)=-rig2(0)*E;
        PanelK(28,9)=-rig2(0)*E;
        PanelK(28,10)=-rig3(0)*E;
   E = theMaterial[1]->getInitialTangent();
//diagonal 2
        PanelK(6,6)=rig1(1)*E;
        PanelK(6,7)=rig2(1)*E;
        PanelK(7,6)=rig2(1)*E;
        PanelK(7,7)=rig3(1)*E;
        PanelK(30,30)=rig1(1)*E;
        PanelK(30,31)=rig2(1)*E;
        PanelK(31,30)=rig2(1)*E;
        PanelK(31,31)=rig3(1)*E;
        PanelK(6,30)=-rig1(1)*E;
        PanelK(6,31)=-rig2(1)*E;
        PanelK(7,30)=-rig2(1)*E;
        PanelK(7,31)=-rig3(1)*E;
        PanelK(30,6)=-rig1(1)*E;
        PanelK(30,7)=-rig2(1)*E;
        PanelK(31,6)=-rig2(1)*E;
        PanelK(31,7)=-rig3(1)*E;
    E = theMaterial[2]->getInitialTangent();
//diagonal 3
       PanelK(12,12)=rig1(2)*E;
        PanelK(12,13)=rig2(2)*E;
        PanelK(13,12)=rig2(2)*E;
        PanelK(13,13)=rig3(2)*E;
        PanelK(24,24)=rig1(2)*E;
        PanelK(24,25)=rig2(2)*E;
        PanelK(25,24)=rig2(2)*E;
        PanelK(25,25)=rig3(2)*E;
        PanelK(12,24)=-rig1(2)*E;
        PanelK(12,25)=-rig2(2)*E;
        PanelK(13,24)=-rig2(2)*E;
        PanelK(13,25)=-rig3(2)*E;
        PanelK(24,12)=-rig1(2)*E;
        PanelK(24,13)=-rig2(2)*E;
        PanelK(25,12)=-rig2(2)*E;
        PanelK(25,13)=-rig3(2)*E;
    E = theMaterial[3]->getInitialTangent();
//diagonal 4
        PanelK(0,0)=rig1(3)*E;
        PanelK(0,1)=rig2(3)*E;
        PanelK(1,0)=rig2(3)*E;
        PanelK(1,1)=rig3(3)*E;
        PanelK(18,18)=rig1(3)*E;
        PanelK(18,19)=rig2(3)*E;
        PanelK(19,18)=rig2(3)*E;
        PanelK(19,19)=rig3(3)*E;
        PanelK(0,18)=-rig1(3)*E;
        PanelK(0,19)=-rig2(3)*E;
        PanelK(1,18)=-rig2(3)*E;
        PanelK(1,19)=-rig3(3)*E;
        PanelK(18,0)=-rig1(3)*E;
        PanelK(18,1)=-rig2(3)*E;
        PanelK(19,0)=-rig2(3)*E;
        PanelK(19,1)=-rig3(3)*E;
    E = theMaterial[4]->getInitialTangent();
//diagonal 5
        PanelK(3,3)=rig1(4)*E;
        PanelK(3,4)=rig2(4)*E;
        PanelK(4,3)=rig2(4)*E;
        PanelK(4,4)=rig3(4)*E;
        PanelK(15,15)=rig1(4)*E;
        PanelK(15,16)=rig2(4)*E;
        PanelK(16,15)=rig2(4)*E;
        PanelK(16,16)=rig3(4)*E;
        PanelK(3,15)=-rig1(4)*E;
        PanelK(3,16)=-rig2(4)*E;
        PanelK(4,15)=-rig2(4)*E;
        PanelK(4,16)=-rig3(4)*E;
        PanelK(15,3)=-rig1(4)*E;
        PanelK(15,4)=-rig2(4)*E;
        PanelK(16,3)=-rig2(4)*E;
        PanelK(16,4)=-rig3(4)*E;
    E = theMaterial[5]->getInitialTangent();
//diagonal 6
        PanelK(21,21)=rig1(5)*E;
        PanelK(21,22)=rig2(5)*E;
        PanelK(22,21)=rig2(5)*E;
        PanelK(22,22)=rig3(5)*E;
        PanelK(33,33)=rig1(5)*E;
        PanelK(33,34)=rig2(5)*E;
        PanelK(34,33)=rig2(5)*E;
        PanelK(34,34)=rig3(5)*E;
        PanelK(21,33)=-rig1(5)*E;
        PanelK(21,34)=-rig2(5)*E;
        PanelK(22,33)=-rig2(5)*E;
        PanelK(22,34)=-rig3(5)*E;
        PanelK(33,21)=-rig1(5)*E;
        PanelK(33,22)=-rig2(5)*E;
        PanelK(34,21)=-rig2(5)*E;
        PanelK(34,22)=-rig3(5)*E;

//		double Et2= theMaterial2->getInitialTangent();
//		double Rigsp=trans(6,3)*Et2;
//				if (WS==0)
//			Rigsp=0;
//		
//		PanelK(0,0)+=Rigsp/2;
//		PanelK(18,18)+=Rigsp/2;
//		PanelK(9,9)+=Rigsp/2;
//		PanelK(27,27)+=Rigsp/2;
//		PanelK(0,18)-=Rigsp/2;
//		PanelK(18,0)-=Rigsp/2;
//		PanelK(9,27)-=Rigsp/2;
//		PanelK(27,9)-=Rigsp/2;
    // return the matrix
    return PanelK;
}

const Vector &
MasonPan12::getResistingForce()
{       

       double Area1=trans(0,3);
      double  c1=trans(0,1);
      double  s1=trans(0,2);
        double force = Area1*theMaterial[0]->getStress();
	PanelR(9)=force*c1;
    PanelR(10)=force*s1;
    PanelR(27)=-force*c1;
    PanelR(28)=-force*s1;

         double Area2=trans(1,3);
        c1=trans(1,1);
        s1=trans(1,2);
      force = Area2*theMaterial[1]->getStress();
	PanelR(6)=force*c1;
    PanelR(7)=force*s1;
    PanelR(30)=-force*c1;
    PanelR(31)=-force*s1;
       double  Area3=trans(2,3);
        c1=trans(2,1);
        s1=trans(2,2);
      force = Area3*theMaterial[2]->getStress();
    PanelR(12)=force*c1;
    PanelR(13)=force*s1;
    PanelR(24)=-force*c1;
    PanelR(25)=-force*s1;
       double  Area4=trans(3,3);
        c1=trans(3,1);
        s1=trans(3,2);
     force = Area4*theMaterial[3]->getStress();

    PanelR(18)=force*c1;
    PanelR(19)=force*s1;
	PanelR(0)=-force*c1;
    PanelR(1)=-force*s1;
      double   Area5=trans(4,3);
        c1=trans(4,1);
        s1=trans(4,2);
     force = Area5*theMaterial[4]->getStress();
    PanelR(3)=-force*c1;
    PanelR(4)=-force*s1;
    PanelR(15)=force*c1;
    PanelR(16)=force*s1;

       double  Area6=trans(5,3);
        c1=trans(5,1);
        s1=trans(5,2);
     force = Area6*theMaterial[5]->getStress();
    PanelR(21)=force*c1;
    PanelR(22)=force*s1;
    PanelR(33)=-force*c1;
    PanelR(34)=-force*s1;

//	double shear;
//	shear=trans(6,0)*theMaterial2->getStress();
//			if (WS==0)
//			shear=0;
//
//	const Vector &disp1 = theNodes[0]->getTrialDisp();
//   const Vector &disp4 = theNodes[3]->getTrialDisp();  
//    const Vector &disp7 = theNodes[6]->getTrialDisp();  
 //   const Vector &disp10 = theNodes[9]->getTrialDisp(); 
//
//	double deltares=((disp7(0)+disp10(0))-(disp1(0)+disp4(0)))/trans(7,0);
//
//	int idir;
//		if ((deltares-Cdeltares)>=0) {
//		  idir=1;
//		  }
//		  else
//		  {
//		  idir=2;
//		 }
//		
//		if ((deltares-Cdeltares)>=0) {
//		 PanelR(9) -= shear;
//		 PanelR(27) += shear;
//		 }
//		 else
//		 {
//		 PanelR(0) += shear;
//		 PanelR(18) -=  shear;
//		 }
    return PanelR;
}

int
MasonPan12::sendSelf(int commitTag, Channel &theChannel)
{
    int res;
	

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    
	//int dataTag = this->getDbTag();

    // Truss2D packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments
	
    //Vector data(5);
    //data(0) = this->getTag();
    //data(1) = A;
    //data(2) = theMaterial->getClassTag();
    //int matDbTag = theMaterial->getDbTag();
    
	// NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    
	//if (matDbTag == 0) {
    //    matDbTag = theChannel.getDbTag();
    //    if (matDbTag != 0)
    //         theMaterial->setDbTag(matDbTag);
   // }
    //data(3) = matDbTag;
	res=0;
    // res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING MasonPan12::sendSelf() - failed to send Vector\n";
      return -1;
    }         

    // Truss2D then sends the tags of it's two end nodes
    // res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING MasonPan12::sendSelf() - failed to send ID\n";
      return -2;
    }

    // finally Truss2D asks it's material object to send itself
    //res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING MasonPan12::sendSelf() - failed to send the Material\n";
      return -3;
    }

    return 0;
}

int
MasonPan12::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res;


   // int dataTag = this->getDbTag();

    // Truss2D creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector

    //Vector data(5);
    //res = theChannel.recvVector(dataTag, commitTag, data);
    res=-1;
	if (res < 0) {
      opserr << "WARNING MasonPan12::recvSelf() - failed to receive Vector\n";
      return -1;
    }         

   // this->setTag((int)data(0));
   // A = data(1);
    
    // Truss2D now receives the tags of it's two external nodes
    //res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING MasonPan12::recvSelf() - failed to receive ID\n";
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    //int matClass = data(2);
    //int matDb = data(3);

    //theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING MasonPan12::recvSelf() - failed to create a Material\n";
      return -3;
    }

    // we set the dbTag before we receive the material  - this is important
    //theMaterial->setDbTag(matDb); 
    //res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING MasonPan12::recvSelf() - failed to receive the Material\n";
      return -3;
    }

    return 0;
}

void
MasonPan12::Print(OPS_Stream &s, int flag)
{

  s << " " << "\n";
  s << " " << "\n";
  s << "Element: " << this->getTag(); 
  s << " type: MasonPan12 " << "\n";
  s << " " << "\n";
  s << "+--------------------------------------------------------+"<< "\n";
  s << "|                 REFINED MASONRY PANEL                  |\n";
  s << "|   Written by Gonzalo Torrisi UNCuyo Copyright 2015     |\n";
  s << "|       Model with 6 compression struts-36 dof           |\n";
  s << "|                 Only in plane X-Y                      |\n";
  s << "|                Use at your Own Peril                   |\n";
  s << "+--------------------------------------------------------+"<<"\n";
  s << "             Nodes: " << "\n";
  s << "Nodo 1  :"<< externalNodes(0)<< "\n";
  s << "Nodo 2  :"<< externalNodes(1)<< "\n";
  s << "Nodo 3  :"<< externalNodes(2)<< "\n";
  s << "Nodo 4  :"<< externalNodes(3)<< "\n";
  s << "Nodo 5  :"<< externalNodes(4)<< "\n";
  s << "Nodo 6  :"<< externalNodes(5)<< "\n";
  s << "Nodo 7  :"<< externalNodes(6)<< "\n";
  s << "Nodo 8  :"<< externalNodes(7)<< "\n";
  s << "Nodo 9  :"<< externalNodes(8)<< "\n";
  s << "Nodo 10 :"<< externalNodes(9)<< "\n";
  s << "Nodo 11 :"<< externalNodes(10)<< "\n";
  s << "Nodo 12 :"<< externalNodes(11)<< "\n";
  s << "        MasonPan12 Factors: " << "\n";
  s << "Panel Thickness                     :"<< TH<< "\n";  
  s << "Factor wd (total strut width)       :"<< WR<< "\n";
  s << "Factor w1 (percent to strut 1)      :"<< W1<< "\n";
//  s << "Factor Gs (percent to shear spring) :"<< WS<< "\n";
  s << "           MasonPan12 Areas: " << "\n";
  s << "Area 1-4 :"<< trans(0,3)<< "\n";
  s << "Area 2-5 :"<< trans(1,3)<< "\n";
  s << "Area 3-6 :"<< trans(2,3)<< "\n";
  s << "Area panel :" << trans(6,0)<< "\n";
   s << "         MasonPan12 Materials: " << "\n";
  s << "Material for central struts :" << *theMaterial[0]<< "\n";
  s << "Material for lateral struts :" << *theMaterial[1]<< "\n";
//  s << "Material for spring :" << *theMaterial2<< "\n";
  s << " " << "\n";
}

Response*
MasonPan12::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Masonpan");
    output.attr("eleTag",this->getTag());
    output.attr("node1 ",externalNodes[0]);
    output.attr("node2 ",externalNodes[1]);
    output.attr("node3 ",externalNodes[2]);
    output.attr("node4 ",externalNodes[3]);
    output.attr("node5 ",externalNodes[4]);
    output.attr("node6 ",externalNodes[5]);
    output.attr("node7 ",externalNodes[6]);
    output.attr("node8 ",externalNodes[7]);
    output.attr("node9 ",externalNodes[8]);
    output.attr("node10",externalNodes[9]);
    output.attr("node11",externalNodes[10]);
    output.attr("node12",externalNodes[11]);	
    char outputData[10];

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            char outputData[10];
 //           int numDOFperNode = numDOF/2;
            for (int i=0; i<3; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<3; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(36));

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
            theResponse = new ElementResponse(this, 3, Vector(6));

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
            for (i=0; i<6; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            for (i=0; i<6; i++) {
                sprintf(outputData,"P%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 4, Vector(2*6));

    // a material quantity
    } else if (strcmp(argv[0],"material") == 0) {
      if (argc > 2) {
	int matNum = atoi(argv[1]);
	if (matNum >= 1 && matNum <= 6)
	  theResponse =  theMaterial[matNum-1]->setResponse(&argv[2], argc-2, output);
		  theResponse =  theMaterial2->setResponse(&argv[2], argc-2, output);
      }
    
	}


    output.endTag();

    return theResponse;
}

int 
MasonPan12::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector& disp3 = theNodes[2]->getTrialDisp();	
    const Vector& disp4 = theNodes[3]->getTrialDisp();	
    const Vector& disp5 = theNodes[4]->getTrialDisp();	
    const Vector& disp6 = theNodes[5]->getTrialDisp();	
    const Vector& disp7 = theNodes[6]->getTrialDisp();
    const Vector& disp8 = theNodes[7]->getTrialDisp();
    const Vector& disp9 = theNodes[8]->getTrialDisp();
    const Vector& disp10 = theNodes[9]->getTrialDisp();
    const Vector& disp11 = theNodes[10]->getTrialDisp();
    const Vector& disp12 = theNodes[11]->getTrialDisp();
	
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    case 2:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < 6; i++)
                (*(eleInformation.theVector))(i) = trans(i,3)*theMaterial[i]->getStress();
  //               (*(eleInformation.theVector))(6) = trans(6,0)*theMaterial2->getStress();              			 
        }
        return 0;

    case 3:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < 6; i++)
                (*(eleInformation.theVector))(i) = theMaterial[i]->getStrain();
  //            (*(eleInformation.theVector))(6) = theMaterial2->getStrain();
        }
        return 0;

    case 13:
        if (eleInformation.theMatrix != 0) {
            for (int i = 0; i < 6; i++)
	      (*(eleInformation.theMatrix))(i,i) = theMaterial[i]->getTangent();
        }
        return 0;

    case 4:
        if (eleInformation.theVector != 0) {
            for (int i = 0; i < 6; i++) {
                (*(eleInformation.theVector))(i) = theMaterial[i]->getStrain();
                (*(eleInformation.theVector))(i+6) = trans(i,3)*theMaterial[i]->getStress();
			}
	//		     (*(eleInformation.theVector))(6) = theMaterial2->getStrain();
  //              (*(eleInformation.theVector))(13) = trans(6,0)*theMaterial2->getStress();
		}
        return 0;      

    default:
        return -1;
    }
}


double
MasonPan12::computeCurrentStrain(int mat) const
//MasonPan12::computeCurrentStrain(void) const
{
			// NOTE this method will never be called with L == 0.0
		double str[6];
		double strain;
		// determine the strain
		const Vector &disp1 = theNodes[0]->getTrialDisp();
		const Vector &disp2 = theNodes[1]->getTrialDisp();  
		const Vector &disp3 = theNodes[2]->getTrialDisp();  
		const Vector &disp4 = theNodes[3]->getTrialDisp();  
		const Vector &disp5 = theNodes[4]->getTrialDisp();  
		const Vector &disp6 = theNodes[5]->getTrialDisp();  
		const Vector &disp7 = theNodes[6]->getTrialDisp();  
		const Vector &disp8 = theNodes[7]->getTrialDisp();  
		const Vector &disp9 = theNodes[8]->getTrialDisp();  
		const Vector &disp10 = theNodes[9]->getTrialDisp(); 
		const Vector &disp11 = theNodes[10]->getTrialDisp();        
		const Vector &disp12 = theNodes[11]->getTrialDisp();        

       double  c1=trans(0,1);
       double s1=trans(0,2);
	double	L1=trans(0,0);
        double du1=((disp4(0)-disp10(0))*c1+(disp4(1)-disp10(1))*s1)/L1;
        c1=trans(1,1);
        s1=trans(1,2);
	double	L2=trans(1,0);
        double du2=((disp3(0)-disp11(0))*c1+(disp3(1)-disp11(1))*s1)/L2;
        c1=trans(2,1);
        s1=trans(2,2);
	double	L3=trans(2,0);
        double du3=((disp5(0)-disp9(0))*c1+(disp5(1)-disp9(1))*s1)/L3;
        c1=trans(3,1);
        s1=trans(3,2);
	double	L4=trans(3,0);
        double du4=((disp7(0)-disp1(0))*c1+(disp7(1)-disp1(1))*s1)/L4;
        c1=trans(4,1);
        s1=trans(4,2);
	double	L5=trans(4,0);
        double du5=((disp6(0)-disp2(0))*c1+(disp6(1)-disp2(1))*s1)/L5;
        c1=trans(5,1);
        s1=trans(5,2);
	double	L6=trans(5,0);
        double du6=((disp8(0)-disp12(0))*c1+(disp8(1)-disp12(1))*s1)/L6;
 
           str[0]=du1;
		   str[1]=du2;
		   str[2]=du3;
		   str[3]=du4;
		   str[4]=du5;
		   str[5]=du6;
  		 strain=str[mat];
		return strain;
}

int
//MasonPan12::displaySelf(Renderer &theViewer, int displayMode, float fact)
MasonPan12::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	  
	int code=0;
	// first determine the two end points of the CorotTruss2 based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	
	const Vector &end3Crd = theNodes[2]->getCrds();	
	const Vector &end4Crd = theNodes[3]->getCrds();	
	const Vector &end5Crd = theNodes[4]->getCrds();	
	const Vector &end6Crd = theNodes[5]->getCrds();	
	const Vector &end7Crd = theNodes[6]->getCrds();	
	const Vector &end8Crd = theNodes[7]->getCrds();	
	const Vector &end9Crd = theNodes[8]->getCrds();	
	const Vector &end10Crd = theNodes[9]->getCrds();	
	const Vector &end11Crd = theNodes[10]->getCrds();	
	const Vector &end12Crd = theNodes[11]->getCrds();		
	
	
	const Vector &end1Disp = theNodes[0]->getDisp();
	const Vector &end2Disp = theNodes[1]->getDisp();   
	const Vector &end3Disp = theNodes[2]->getDisp();
	const Vector &end4Disp = theNodes[3]->getDisp();
	const Vector &end5Disp = theNodes[4]->getDisp();
	const Vector &end6Disp = theNodes[5]->getDisp();
	const Vector &end7Disp = theNodes[6]->getDisp();
	const Vector &end8Disp = theNodes[7]->getDisp();
	const Vector &end9Disp = theNodes[8]->getDisp();
	const Vector &end10Disp = theNodes[9]->getDisp();
	const Vector &end11Disp = theNodes[10]->getDisp();	
	const Vector &end12Disp = theNodes[11]->getDisp();	
    
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	static Vector v9(3);
	static Vector v10(3);
	static Vector v11(3);
	static Vector v12(3);	
	static Vector v13(3);
	static Vector v14(3);

	static Vector v1a(3);
	static Vector v2a(3);
	static Vector v3a(3);
	static Vector v4a(3);
	static Vector v5a(3);
	static Vector v6a(3);
	static Vector v7a(3);
	static Vector v8a(3);
	static Vector v9a(3);
	static Vector v10a(3);
	static Vector v11a(3);
	static Vector v12a(3);




	static Vector rgb(3);

	theNodes[3]->getDisplayCrds(v1a, fact, displayMode);
	theNodes[9]->getDisplayCrds(v2a, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3a, fact, displayMode);
	theNodes[10]->getDisplayCrds(v4a, fact, displayMode);
	theNodes[4]->getDisplayCrds(v5a, fact, displayMode);
	theNodes[8]->getDisplayCrds(v6a, fact, displayMode);
	theNodes[6]->getDisplayCrds(v7a, fact, displayMode);
	theNodes[0]->getDisplayCrds(v8a, fact, displayMode);
	theNodes[5]->getDisplayCrds(v9a, fact, displayMode);
	theNodes[1]->getDisplayCrds(v10a, fact, displayMode);
	theNodes[7]->getDisplayCrds(v11a, fact, displayMode);
	theNodes[11]->getDisplayCrds(v12a, fact, displayMode);


	for (int i = 0; i < 2; i++) {
		v1(i) = end4Crd(i)+end4Disp(i)*fact;
		v2(i) = end10Crd(i)+end10Disp(i)*fact;    
		v3(i) = end3Crd(i)+end3Disp(i)*fact;  
		v4(i) = end11Crd(i)+end11Disp(i)*fact;  	
		v5(i) = end5Crd(i)+end5Disp(i)*fact;  
		v6(i) = end9Crd(i)+end9Disp(i)*fact;  
		
		v7(i) = end7Crd(i)+end7Disp(i)*fact;
		v8(i) = end1Crd(i)+end1Disp(i)*fact;    
		v9(i) = end6Crd(i)+end6Disp(i)*fact;  
		v10(i) = end2Crd(i)+end2Disp(i)*fact;  	
		v11(i) = end8Crd(i)+end8Disp(i)*fact;  
		v12(i) = end12Crd(i)+end12Disp(i)*fact;  	

//		v13(0)=(v10(0)+v6(0))/2;
//		v14(0)=(v3(0)+v11(0))/2;
//
//		v13(1)=(v8(1)+v2(1))/2;
//		v14(1)=(v1(1)+v7(1))/2;
//		v13(2)=0.0;
//		v14(2)=0.0;
	}

	// compute the strain and axial force in the member
	double strain[7], force[7];

	for (int i = 0; i < 6; i++) {

	    strain[i] = this->computeCurrentStrain(i);
	    theMaterial[i]->setTrialStrain(strain[i]);
	    force[i] = theMaterial[i]->getStress();    

	}
//	    strain[6]=(end10Disp(0)+end7Disp(0))-(end1Disp(0)+end4Disp(0));
//		theMaterial2->setrialStrain(strain[6],FFAXIAL);
//		force[6]=theMaterial2->getStress();
	

//				code=0;
//		  code +=theViewer.drawLine(v1,v2,1.0,1.0);
//		  code +=theViewer.drawLine(v3,v4,1.0,1.0);
//	      code +=theViewer.drawLine(v5,v6,1.0,1.0);
//		  code +=theViewer.drawLine(v7,v8,1.0,1.0);
//		  code +=theViewer.drawLine(v9,v10,1.0,1.0);
//	      code +=theViewer.drawLine(v11,v12,1.0,1.0);
//	      code +=theViewer.drawLine(v13,v14,(float)strain[6],(float)strain[65]);
//		  return code;



	if (displayMode == 2) // use the strain as the drawing measure
		{
			code=0;
		  code +=theViewer.drawLine(v1,v2,(float)strain[0],(float)strain[0]);
		  code +=theViewer.drawLine(v3,v4,(float)strain[1],(float)strain[1]);
	      code +=theViewer.drawLine(v5,v6,(float)strain[2],(float)strain[2]);
		  code +=theViewer.drawLine(v7,v8,(float)strain[3],(float)strain[3]);
		  code +=theViewer.drawLine(v9,v10,(float)strain[4],(float)strain[4]);
	      code +=theViewer.drawLine(v11,v12,(float)strain[5],(float)strain[5]);

		  return code;
}

	else if (displayMode < 0)
	{
		code = 0;
		code += theViewer.drawLine(v1a, v2a, 1.0, 1.0, this->getTag(), 0);
		code += theViewer.drawLine(v3a, v4a, 1.0, 1.0, this->getTag(), 0);
		code += theViewer.drawLine(v5a, v6a, 1.0, 1.0, this->getTag(), 0);
		code += theViewer.drawLine(v7a, v8a, 1.0, 1.0, this->getTag(), 0);
		code += theViewer.drawLine(v9a, v10a, 1.0, 1.0, this->getTag(), 0);
		code += theViewer.drawLine(v11a, v12a, 1.0, 1.0, this->getTag(), 0);
		return code;

	}

	

	else { // otherwise use the axial force as measure
		code=0;
		  //code +=theViewer.drawLine(v1,v2,(float)force[0],(float)force[0],0,0,2,0);
		  code += theViewer.drawLine(v1, v2, (float)force[0], (float)force[0]);
		  code +=theViewer.drawLine(v3,v4,(float)force[1],(float)force[1]);
	      code +=theViewer.drawLine(v5,v6,(float)force[2],(float)force[2]);
		  //code +=theViewer.drawLine(v7,v8,(float)force[3],(float)force[3],0,0,2,0);
		  code += theViewer.drawLine(v7, v8, (float)force[3], (float)force[3]);
		  code +=theViewer.drawLine(v9,v10,(float)force[4],(float)force[4]);
	      code +=theViewer.drawLine(v11,v12,(float)force[5],(float)force[5]);


		  return code;
	}

	return 0;
}

 
