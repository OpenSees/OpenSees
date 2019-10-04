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
                                                                        
// $Revision: 1.10 $
// $Date: 2014/07/01 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shellDKGQ/shellDKGQ.h,v $

// Written: Shuhao Zhang & Xinzheng Lu
//
// Three node flat shell element with membrane and drill DOF
// Ref: Plate Bending Part - DKT, thin plate element
//      Membrane Part - GT9, a membrane element with drilling DOF
//

/// Jose Abell added support for self weight 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellDKGT.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <DataFileStream.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#include <ElementalLoad.h>

#define min(a,b) ( (a)<(b) ? (a):(b) )

static int numShellDKGT = 0;

void *
OPS_ShellDKGT(void)         
{
  if (numShellDKGT == 0) {
//    opserr << "Using ShellDKGT - Developed by: Shuhao Zhang and Xinzheng Lu\n";
    numShellDKGT++;
  }

  Element *theElement = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 5) {
    opserr << "Want: element ShellDKGT $tag $iNode $jNoe $kNode $secTag";   
    return 0;	
  }
  
  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKGT \n";
    return 0;
  }

  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(iData[4]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellDKGT " << iData[0] << "section " << iData[4] << " not found\n";
    return 0;
  }
  

    double b_data[3] = {0, 0, 0};

    int num_remaining_args = OPS_GetNumRemainingInputArgs();
    if (num_remaining_args > 3) {
        num_remaining_args = 3;
    }
    if (num_remaining_args > 0) {
        if (OPS_GetDoubleInput(&num_remaining_args, b_data) < 0) {
            opserr << "WARNING: invalid double b_data\n";
            return 0;
        }
    }

  theElement = new ShellDKGT(iData[0], iData[1], iData[2], iData[3],
                               *theSection, b_data[0], b_data[1], b_data[2]);

  return theElement;
}


//static data
Matrix  ShellDKGT::stiff(18,18) ;                   
Vector  ShellDKGT::resid(18) ;
Matrix  ShellDKGT::mass(18,18) ;

//quadrature data

const double  ShellDKGT::one_over_three =1.0/3.0 ;
const double  ShellDKGT::one_over_five = 1.0/5.0 ;
const double ShellDKGT:: three_over_five= 3.0/5.0 ;
const double ShellDKGT:: wg1= -9.0/16.0;
const double ShellDKGT:: wg2= 25.0/48.0;


double ShellDKGT::sg[4] ;
double ShellDKGT::tg[4] ;
double ShellDKGT::qg[4] ;
double ShellDKGT::wg[4] ;

 

//null constructor
ShellDKGT::ShellDKGT( ) :                            
Element( 0, ELE_TAG_ShellDKGT ),
connectedExternalNodes(3), load(0), Ki(0)
{ 
  for (int i = 0 ;  i < 4; i++ ) 
    materialPointers[i] = 0;

  sg[0] = one_over_three;
  sg[1] = one_over_five;
  sg[2] = three_over_five;
  sg[3] = one_over_five;
 

  tg[0] = one_over_three;
  tg[1] = three_over_five;
  tg[2] = one_over_five;
  tg[3] = one_over_five;

  qg[0] = one_over_three;
  qg[1] = one_over_five;
  qg[2] = one_over_five;
  qg[3] = three_over_five;
 

  wg[0] =wg1;
  wg[1] =wg2;
  wg[2] =wg2;
  wg[3] =wg2;

    b[0] = 0;
    b[1] = 0;
    b[2] = 0;
 
}


//*********************************************************************
//full constructor
ShellDKGT::ShellDKGT(  int tag,                      
                         int node1,
                         int node2,
   	                     int node3,
                       SectionForceDeformation &theMaterial, double b1, double b2, double b3 ) :
Element( tag, ELE_TAG_ShellDKGT ),
connectedExternalNodes(3), load(0), Ki(0)
{
  int i;
  connectedExternalNodes(0) = node1 ;           
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;

  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy( ) ;      

      if (materialPointers[i] == 0) {
	      opserr << "ShellDKGT::constructor - failed to get a material of type: ShellSection\n";
      } //end if
      
  } //end for i 

  sg[0] = one_over_three;
  sg[1] =  one_over_five;
  sg[2] = three_over_five;
  sg[3] =one_over_five;
 

  tg[0] = one_over_three;
  tg[1] = three_over_five;
  tg[2] = one_over_five;
  tg[3] = one_over_five;

  qg[0] = one_over_three;
  qg[1] = one_over_five;
  qg[2] = one_over_five;
  qg[3] = three_over_five;
 

  wg[0] = wg1;
  wg[1] = wg2;
  wg[2] = wg2;
  wg[3] = wg2;
  

    b[0] = b1;
    b[1] = b2;
    b[2] = b3;

 }
//******************************************************************

//destructor 
ShellDKGT::~ShellDKGT( )                         
{
  int i ;
  for ( i = 0 ;  i <4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ;   

  } //end for i

  for (i=0;i<3;i++){
  nodePointers[i] = 0 ;
  }

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
}
//**************************************************************************
//set domain
void  ShellDKGT::setDomain( Domain *theDomain )                             
{  
  int i;

  //node pointers
  for ( i = 0; i < 3; i++ ) {
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;               
     if (nodePointers[i] == 0) {
       opserr << "ShellDKGT::setDomain - no node " << connectedExternalNodes(i);
       opserr << " exists in the model\n";
     }
     const Vector &nodeDisp=nodePointers[i]->getTrialDisp();
     if (nodeDisp.Size() != 6) {
       opserr << "ShellDKGT::setDomain - node " << connectedExternalNodes(i);
       opserr << " NEEDS 6 dof - GARBAGE RESULTS or SEGMENTATION FAULT WILL FOLLOW\n";
     }       
  }

  //basis vectors and local coordinates
  computeBasis( ) ;

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  ShellDKGT::getNumExternalNodes( ) const
{
  return 3 ;
} 
 

//return connected external nodes
const ID&  ShellDKGT::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


Node **
ShellDKGT::getNodePtrs(void) 
{
  return nodePointers;
} 

//return number of dofs
int  ShellDKGT::getNumDOF( ) 
{
  return 18 ;
}


//commit state
int  ShellDKGT::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ShellDKGT::commitState () - failed in base class";
  }    
  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;

  return success ;
}
 


//revert to last commit 
int  ShellDKGT::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  ShellDKGT::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i <4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}


void  ShellDKGT::Print( OPS_Stream &s, int flag )
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ShellDKGQ\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
        s << "\t" << connectedExternalNodes(2) << "\t" << "\t0.00";
        s << endln;
        s << "PROP_3D\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << -1 << "\tSHELL\t1.0\0.0";
        s << endln;
    }
    
    else if (flag < -1) {
        
        int counter = (flag + 1) * -1;
        int eleTag = this->getTag();
        int i, j;
        for (i = 0; i < 4; i++) {
            const Vector &stress = materialPointers[i]->getStressResultant();
            
            s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
            for (j = 0; j < 6; j++)
                s << "\t" << stress(j);
            s << endln;
        }
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "DKGT Non-Locking Three Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << connectedExternalNodes(0) << endln;
        s << "Node 2 : " << connectedExternalNodes(1) << endln;
        s << "Node 3 : " << connectedExternalNodes(2) << endln;
        
        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);
        
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ShellDKGT\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        s << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << "], ";
        s << "\"section\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}

Response*
ShellDKGT::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");

  output.attr("eleType", "ShellDKGT");
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
  } 

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"Material") == 0) {
    if (argc < 2) {
      opserr << "ShellDKGT::setResponse() - need to specify more data\n";
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {                   
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",sg[pointNum-1]);
      output.attr("neta",tg[pointNum-1]);
      
      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();
    }

  } else if (strcmp(argv[0],"stresses") ==0) {

    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",sg[i]);
      output.attr("neta",tg[i]);
      
      output.tag("SectionForceDeformation");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      
      output.tag("ResponseType","p11");
      output.tag("ResponseType","p22");
      output.tag("ResponseType","p1212");
      output.tag("ResponseType","m11");
      output.tag("ResponseType","m22");
      output.tag("ResponseType","m12");
      output.tag("ResponseType","q1");
      output.tag("ResponseType","q2");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }
    
    theResponse =  new ElementResponse(this, 2, Vector(32));
  }
  
  else if (strcmp(argv[0],"strains") ==0) {

    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",sg[i]);
      output.attr("neta",tg[i]);
      
      output.tag("SectionForceDeformation");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());
      
      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","gamma12");
      output.tag("ResponseType","theta11");
      output.tag("ResponseType","theta22");
      output.tag("ResponseType","theta33");
      output.tag("ResponseType","gamma13");
      output.tag("ResponseType","gamma23");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }
    
    theResponse =  new ElementResponse(this, 3, Vector(32));
  }

  output.endTag();
  return theResponse;
}

int
ShellDKGT::getResponse(int responseID, Information &eleInfo)
{
  int cnt = 0;

  static Vector stresses(32);
  static Vector strains(32);

  switch (responseID) {
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
    break;

  case 2: // stresses
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStressResultant();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      stresses(cnt+3) = sigma(3);
      stresses(cnt+4) = sigma(4);
      stresses(cnt+5) = sigma(5);
      stresses(cnt+6) = sigma(6);
      stresses(cnt+7) = sigma(7);
      cnt += 8;
    }
    return eleInfo.setVector(stresses);
    break;
  case 3: //strain
    for (int i = 0; i < 4; i++) {

      // Get section deformation
      const Vector &deformation = materialPointers[i]->getSectionDeformation();
      strains(cnt) = deformation(0);
      strains(cnt+1) = deformation(1);
      strains(cnt+2) = deformation(2);
      strains(cnt+3) = deformation(3);
      strains(cnt+4) = deformation(4);
      strains(cnt+5) = deformation(5);
      strains(cnt+6) = deformation(6);
      strains(cnt+7) = deformation(7);
      cnt += 8;
    }
    return eleInfo.setVector(strains);
    break;
  default:
    return -1;
  }
  cnt=0;

  //return 0;
}


//return stiffness matrix 
const Matrix&  ShellDKGT::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ; 
 // opserr << stiff;

  return stiff ;
}    


const Matrix&  ShellDKGT::getInitialStiff( ) 
{
	if(Ki!=0)
		return *Ki;

	static const int ndf = 6; //two membrane plus three bending plus one drill

	static const int nstress = 8; //three membrane, three moment, two shear

	static const int ngauss = 4;

	static const int numnodes = 3;  ///node number

	int i,j,k,p,q;
	int jj,kk;

	int p1,q1;

	int p2,q2;

	//int success;

	double volume = 0.0;

	static double xsj; //determinant jacabian matrix

	static double dvol[ngauss]; //volume element

	//static Vector strain(nstress);//strain

	static double shp[3][numnodes]; //shape function 2d at a gauss point

//	static double shpM[3][numnodes];//shape function-membrane at a gausss point

	static double shpDrill[4][numnodes];//shape function-drilling dof(Nu,1&Nu,2&Nv,1&Nv,2) at a gauss point

	static double shpBend[6][9];//shape function -bending part(Hx,Hy,Hx-1,2&Hy-1,2) at a gauss point

	//static Vector residJ(ndf,ndf); //nodeJ residual, global coordinates

	static Matrix stiffJK(ndf,ndf);//nodeJK stiffness, global coordinates

	//static Vector residJlocal(ndf,ndf); // nodeJ residual, local coordinates

	static Matrix stiffJKlocal(ndf,ndf);//nodeJK stiffness, local coordinates

	static Matrix stiffJK1(ndf,ndf);

	static Matrix stiffJK2(ndf,ndf);

	static Matrix stiffJK3(ndf,ndf);

	//static Vector stress(nstress); //stress resultants

	static Matrix dd(nstress,nstress); // material tangent

	static double sx[2][2]; // inverse of Jacobian

	Matrix Tmat(6,6);
	Matrix TmatTran(6,6);
	Matrix Pmat(6,6);
	Matrix PmatTran(6,6);

	//--------------------B-matrices--------------------------------------

	static Matrix BJ(nstress, ndf); // B matrix node J

	static Matrix BJtran(ndf, nstress);

	static Matrix BK(nstress, ndf); // B matrix node K

	static Matrix BJtranD(ndf, nstress); //BJtran * dd


	static Matrix Bmembrane(3,3); //membrane B matrix

	static Matrix Bbend(3,3); //bending B matrix

	static Matrix Bshear(2,3); //shear B matrix (zero)

	static double saveB[nstress][ndf][numnodes];

	//--------------------------------------------------------

	stiff.Zero( );

	//define Pmat- transpose the dofs               
	Pmat.Zero();
	double one=1.00;
	Pmat(0,0) = one;
	Pmat(1,1) = one;
	Pmat(2,5) = one;
	Pmat(3,2) = one;
	Pmat(4,3) = one;
	Pmat(5,4) = one;
   //transpose PmatTran=transpose(Pmat)
	for (p1=0;p1<6;p1++){
		for (q1=0;q1<6;q1++){
			PmatTran(p1,q1) = Pmat(q1,p1);
		}
	}//end for p1

	//define Tmat xl=Tmat * x - from global to local coordinates
	Tmat.Zero();
	Tmat(0,0) = g1[0];
	Tmat(0,1) = g2[0];
	Tmat(0,2) = g3[0];
	Tmat(1,0) = g1[1];
	Tmat(1,1) = g2[1];
	Tmat(1,2) = g3[1];
	Tmat(2,0) = g1[2];
	Tmat(2,1) = g2[2];
	Tmat(2,2) = g3[2];

	Tmat(3,3) = g1[0];
	Tmat(3,4) = g2[0];
	Tmat(3,5) = g3[0];
	Tmat(4,3) = g1[1];
	Tmat(4,4) = g2[1];
	Tmat(4,5) = g3[1];
	Tmat(5,3) = g1[2];
	Tmat(5,4) = g2[2];
	Tmat(5,5) = g3[2];


   //transpose TmatTran=transpose(Tmat)
	for (p2=0;p2<6;p2++){
		for (q2=0;q2<6;q2++){
			TmatTran(p2,q2) = Tmat(q2,p2);
		}
	}//end for p2


    

	//------------gauss loop--------------------------
	for (i=0; i<ngauss;i++){

		//get shape functions
		shape2d(sg[i],tg[i],qg[i],xl,shp,xsj,sx);

		shapeDrill(sg[i],tg[i],qg[i],xl,sx,shpDrill);

		shapeBend(sg[i],tg[i],qg[i],xl,sx,shpBend);

		//volume element to be saved
		dvol[i] =0.5* wg[i] * xsj;
		volume += dvol[i];

		Bshear.Zero();

		// j-node loop to compute strain
		for (j = 0; j < numnodes; j++){

			//compute B matrix

			Bmembrane = computeBmembrane(j,shp,shpDrill);

			Bbend = computeBbend(j,shpBend);

			BJ = assembleB(Bmembrane, Bbend, Bshear);

			//save the B-matrix
			for (p=0;p<nstress; p++){
				for(q=0;q<ndf;q++){
					saveB[p][q][j] = BJ(p,q);
				}
			}//end for p

		}//end j-node loop

		dd = materialPointers[i]->getInitialTangent( );
		dd *= dvol[i];

		//tangent stiff matrix calculations node loops

		jj = 0;
		for (j=0; j<numnodes;j++){

			//extract BJ
			for(p=0;p<nstress;p++){
				for (q=0; q<ndf;q++)
					BJ(p,q)=saveB[p][q][j];
			}// end for p
			//multiply bending terms by -1.0 for correct statement of equilibrium
			for (p=3;p<6;p++){
				for (q=3; q<6; q++)
					BJ(p,q) *= (-1.0);
			}//end for p

			//transpose BJtran=transpose(BJ);
			for (p=0;p<ndf;p++){
				for(q=0;q<nstress;q++)
					BJtran(p,q)=BJ(q,p);
			}//end for p

			//BJtranD = BJtran *dd;
			BJtranD.addMatrixProduct(0.0,BJtran,dd,1.0);

			//k loop
			kk=0;
			for(k=0; k<numnodes; k++){

				//extract BK
				for (p=0; p<nstress;p++){
					for(q=0;q<ndf;q++)
						BK(p,q) = saveB[p][q][k];
				}// end for p

				//stiffJKlocal = BJtranD * BK
				stiffJKlocal.addMatrixProduct(0.0,BJtranD,BK,1.0);

				//transpose dof and coordinates
				stiffJK1.addMatrixProduct(0.0,PmatTran,stiffJKlocal,1.0);
				stiffJK2.addMatrixProduct(0.0,stiffJK1,Pmat,1.0);
				stiffJK3.addMatrixProduct(0.0,TmatTran,stiffJK2,1.0);
				stiffJK.addMatrixProduct(0.0,stiffJK3,Tmat,1.0);

				//generate stiff
				for(p=0;p<ndf;p++){
					for(q=0;q<ndf;q++){
						stiff(jj+p,kk+q) += stiffJK(p,q);
					}//end for q
				}//end for p

				kk += ndf;
			}//end for k loop

			jj += ndf;

		}// end for j loop

	}//end for i gauss loop

	Ki = new Matrix(stiff);

	return stiff;

}





//**********************************************************************
// shape function for drill dof


//return mass matrix
const Matrix&  ShellDKGT::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


void  ShellDKGT::zeroLoad( )
{
  if (load != 0)
    load->Zero();
 
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;
    appliedB[2] = 0.0;

  return ;
}


int 
ShellDKGT::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
        applyLoad = 1;
        appliedB[0] += loadFactor * data(0) ;
        appliedB[1] += loadFactor * data(1) ;
        appliedB[2] += loadFactor * data(2) ;
        return 0;
    }

  opserr << "ShellDKGT::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}



int 
ShellDKGT::addInertiaLoadToUnbalance(const Vector &accel)
{
  int tangFlag = 1 ;

  int i;

  int allRhoZero = 0;
  for (i=0; i<4; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      allRhoZero = 1;
  }

  if (allRhoZero == 0) 
    return 0;

  int count = 0;
  for (i=0; i<3; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<6; j++)
      resid(count++) = Raccel(i);
  }

  formInertiaTerms( tangFlag ) ;
  if (load == 0) 
    load = new Vector(18);
  load->addMatrixVector(1.0, mass, resid, -1.0);

  return 0;
}



//get residual
const Vector&  ShellDKGT::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  ShellDKGT::getResistingForceIncInertia( )
{
  static Vector res(18);
  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  formInertiaTerms( tang_flag ) ;

  res = resid;
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads 
  if (load != 0)
    res -= *load;

  return res;
}

//*********************************************************************
//form inertia terms

void   
ShellDKGT::formInertiaTerms( int tangFlag ) 

{

  //translational mass only
  //rotational inertia terms are neglected


  static const int ndf = 6 ; 

  static const int numberNodes = 3 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double sx[2][2]; //inverse jacobian matrix
 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;


  int i, j, k, p;
  int jj, kk ;

  double temp, rhoH, massJK ;


  //zero mass 
  mass.Zero( ) ;


  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i],qg[i],xl,shp,xsj,sx) ;

    //volume element to also be saved
    dvol = 0.5*wg[i]*xsj;  


    //node loop to compute accelerations
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += ( shp[massIndex][j] * nodePointers[j]->getTrialAccel() ) ;
      momentum.addVector(1.0,  
			 nodePointers[j]->getTrialAccel(),
			 shp[massIndex][j] ) ;

      
    //density
    rhoH = materialPointers[i]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rhoH ;


    //residual and tangent calculations node loops
    for ( j=0, jj=0; j<numberNodes; j++, jj+=ndf ) {

      temp = shp[massIndex][j] * dvol ;

      for ( p = 0; p < 3; p++ )
        resid( jj+p ) += ( temp * momentum(p) ) ;

      
      if ( tangFlag == 1 && rhoH != 0.0) {

	 //multiply by density
	 temp *= rhoH ;

	 //node-node translational mass
         for ( k=0, kk=0; k<numberNodes; k++, kk+=ndf ) {

	   massJK = temp * shp[massIndex][k] ;

	   for ( p = 0; p < 3; p++ ) 
	      mass( jj+p, kk+p ) +=  massJK ;
            
          } // end for k loop

      } // end if tang_flag 

    } // end for j loop

  } //end for i gauss loop 

}

//*********************************************************************

//form residual and tangent
void  
ShellDKGT::formResidAndTangent( int tang_flag ) 
{
	//
	//six(6) nodal dof's ordered:
	//-----------
	//|    u1    |   <---- plate membrane
	//|    u2    |
	//|----------|
	//|   w = u3 |  <----plate bending
	//|   theta1 |
	//|   theta2 |
	//|----------|
	//|   theta3 |  <- drill (tran from membrane)
	//|----------|
 
	static const int ndf = 6; //two membrane + 3 moment +drill

	static const int nstress = 8; //3 membrane , 3 moment, 2 shear

	static const int ngauss = 4;

	static const int numnodes = 3;

	int i,j,k,p,q;
	int jj,kk;

	int p1,q1;

	int p2,q2;

	int success;

	double volume = 0.0;

	static double xsj; //determinant jacobian matrix

	static double dvol[ngauss]; //volume element

	static Vector strain(nstress); //strain

	static double shp[3][numnodes]; //shape function 2d at a gauss point

	static double shpDrill[4][numnodes]; //shape function drilling dof at a gauss point

	static double shpBend[6][9]; //shape function - bending part at a gauss point

	static Vector residJ(ndf); //nodeJ residual, global coordinates

	static Matrix stiffJK(ndf,ndf);//nodeJK stiffness, global coordinates

	static Vector residJlocal(ndf); //nodeJ residual, local coordinates

	static Matrix stiffJKlocal(ndf,ndf); //nodeJK stiffness, local coordinates

	static Matrix stiffJK1(ndf,ndf);

	static Matrix stiffJK2(ndf,ndf);

	static Matrix stiffJK3(ndf,ndf);

	static Vector residJ1(ndf);

	static Vector stress(nstress); //stress resultants

	static Matrix dd(nstress,nstress);//material tangent

	//static Matrix J0(2,2); //Jacobian at center

	//static Matrix J0inv(2,2); //inverse of Jacobian at center
	static double sx[2][2];

	Matrix Tmat(6,6);  //local-global coordinates transform matrix

	Matrix TmatTran(6,6);

	Matrix Pmat(6,6); //transform dofs order

	Matrix PmatTran(6,6);

	//-------------------B-matrices---------------------------------

	static Matrix BJ(nstress, ndf); // B matrix node J

	static Matrix BJtran(ndf, nstress);
	static Matrix BJt(nstress, ndf);

	static Matrix BK(nstress, ndf); // B matrix node K

	static Matrix BJtranD(ndf, nstress); //BJtran * dd

	static Matrix BJP(nstress, ndf); //BJ * Pmat, transform the dof order

	static Matrix BJPT(nstress,ndf); //BJP * Tmmat, from global coordinates to local coordinates

	static Matrix Bmembrane(3,3); //membrane B matrix

	static Matrix Bbend(3,3); //bending B matrix

	static Matrix Bshear(2,3); //shear B matrix (zero)

	static double saveB[nstress][ndf][numnodes];
	//---------------------------------------------------------------

	//zero stiffness and residual
	stiff.Zero( );
	resid.Zero( );

	//start Yuli Huang & Xinzheng Lu
	//updateBasis( );
	computeBasis( );

	
	//end Yuli Huang & Xinzheng Lu

	//define Pmat- transpose the dofs
	Pmat.Zero();
	double one=1.00;
	Pmat(0,0) = one;
	Pmat(1,1) = one;
	Pmat(2,5) = one;
	Pmat(3,2) = one;
	Pmat(4,3) = one;
	Pmat(5,4) = one;
   //transpose PmatTran=transpose(Pmat)
	for (p1=0;p1<6;p1++){
		for (q1=0;q1<6;q1++){
			PmatTran(p1,q1) = Pmat(q1,p1);
		}
	}//end for p1

	//define Tmat xl=Tmat * x from global to local coordinates
	Tmat.Zero();
	Tmat(0,0) = g1[0];
	Tmat(0,1) = g1[1];
	Tmat(0,2) = g1[2];
	Tmat(1,0) = g2[0];
	Tmat(1,1) = g2[1];
	Tmat(1,2) = g2[2];
	Tmat(2,0) = g3[0];
	Tmat(2,1) = g3[1];
	Tmat(2,2) = g3[2];

	Tmat(3,3) = g1[0];
	Tmat(3,4) = g1[1];
	Tmat(3,5) = g1[2];
	Tmat(4,3) = g2[0];
	Tmat(4,4) = g2[1];
	Tmat(4,5) = g2[2];
	Tmat(5,3) = g3[0];
	Tmat(5,4) = g3[1];
	Tmat(5,5) = g3[2];

   //transpose TmatTran=transpose(Tmat)
	for (p2=0;p2<6;p2++){
		for (q2=0;q2<6;q2++){
			TmatTran(p2,q2) = Tmat(q2,p2);
		}
	}//end for p2
   
	//------------gauss loop--------------------------
	for (i=0; i<ngauss;i++){

		//get shape functions
		shape2d(sg[i],tg[i],qg[i],xl,shp,xsj,sx);

		shapeDrill(sg[i],tg[i],qg[i],xl,sx,shpDrill);

		shapeBend(sg[i],tg[i],qg[i],xl,sx,shpBend);

		//volume element to be saved
		dvol[i] =0.5 * wg[i] * xsj;
//		opserr<<dvol[i]<<endln;
		volume += dvol[i];

		Bshear.Zero( );

		//zero the strains
		strain.Zero( );
		
	
		// j-node loop to compute strain
		for (j = 0; j < numnodes; j++){

			//compute B matrix
			Bmembrane = computeBmembrane(j,shp,shpDrill);
			Bbend = computeBbend(j,shpBend);
			BJ = assembleB(Bmembrane, Bbend, Bshear);	


	         
			//save the B-matrix
			for (p=0;p<nstress; p++){
				for(q=0;q<ndf;q++){
					saveB[p][q][j] = BJ(p,q);
				}
			}//end for p


			//nodal "displacements"
			const Vector &ul = nodePointers[j]->getTrialDisp( );

			//compute the strain 
			//Note: transform the dof's order
			//BJP = BJ * P;
			BJP.addMatrixProduct(0.0, BJ, Pmat,1.0);
			BJPT.addMatrixProduct(0.0,BJP,Tmat,1.0);
			//strain += (BJ*ul);
			strain.addMatrixVector(1.0, BJPT, ul,1.0);
			
		}//end j-node loop
		
		//send the strain to the material
		success = materialPointers[i]->setTrialSectionDeformation(strain);

		//compute the stress
		stress = materialPointers[i]->getStressResultant( );


		//multiply by volume element
		stress *= dvol[i];

		if(tang_flag == 1){
			dd = materialPointers[i]->getSectionTangent( );
         
			dd *= dvol[i];
 
		}//end if tang_flag
	   
		//residual and tangent calculations node loops

		jj = 0;
		for (j=0; j<numnodes;j++){

			//extract BJ
			for(p=0;p<nstress;p++){
				for (q=0; q<ndf;q++)
					BJ(p,q)=saveB[p][q][j];
			}// end for p
			//multiply bending terms by -1.0 for correct statement of equilibrium
			for (p=3;p<6;p++){
				for (q=3; q<6; q++)
					BJ(p,q) *= (-1.0);
			}//end for p
		
			//transpose BJtran=transpose(BJ);
			for (p=0;p<ndf;p++){
				for(q=0;q<nstress;q++)
					BJtran(p,q)=BJ(q,p);
			}//end for p

			//compute residual force
			residJlocal.addMatrixVector(0.0,BJtran,stress,1.0);	
			residJ1.addMatrixVector(0.0,PmatTran,residJlocal,1.0);
			residJ.addMatrixVector(0.0,TmatTran,residJ1,1.0);

			for(p=0; p<ndf; p++)
            {
				resid(jj+p) += residJ(p);
            }


			//BJtranD = BJtran *dd;
			//BJtranD.addMatrixProduct(0.0,BJtran,dd,1.0);
			if(tang_flag == 1){

            BJtranD.addMatrixProduct(0.0,BJtran,dd,1.0);
				
			//k loop
			kk=0;
			for(k=0; k<numnodes; k++){

				//extract BK
				for (p=0; p<nstress;p++){
					for(q=0;q<ndf;q++)
						BK(p,q) = saveB[p][q][k];
				}// end for p


				//stiffJKlocal = BJtranD * BK
				stiffJKlocal.addMatrixProduct(0.0,BJtranD,BK,1.0);
				    

				//transpose dof and coordinates
				stiffJK1.addMatrixProduct(0.0,PmatTran,stiffJKlocal,1.0);
				stiffJK2.addMatrixProduct(0.0,stiffJK1,Pmat,1.0);
				stiffJK3.addMatrixProduct(0.0,TmatTran,stiffJK2,1.0);
				stiffJK.addMatrixProduct(0.0,stiffJK3,Tmat,1.0);

				
              //opserr<<Bbend<<endln;
				//generate stiff
				for(p=0;p<ndf;p++){
					for(q=0;q<ndf;q++){
						stiff(jj+p,kk+q) += stiffJK(p,q);
					}//end for q
				}//end for p

				kk += ndf;
			}//end for k loop

			}//end if tang_flag

			jj += ndf;			

		}// end for j loop

	}//end for i gauss loop

	
      if(applyLoad == 1)
      {
          const int numberGauss = 4 ;
          const int nShape = 3 ;
          const int numberNodes = 4 ;
          const int massIndex = nShape - 1 ;
          double temp, rhoH;
          //If defined, apply self-weight
          static Vector momentum(ndf) ;
          double ddvol = 0;
          for ( i = 0; i < numberGauss; i++ ) {

              //get shape functions    
              // shape2d( sg[i], tg[i], xl, shp, xsj ) ;
              shape2d( sg[i], tg[i], qg[i], xl, shp, xsj, sx) ;

              //volume element to also be saved
              ddvol = wg[i] * xsj ;  


              //node loop to compute accelerations
              momentum.Zero( ) ;
              momentum(0) = appliedB[0];
              momentum(1) = appliedB[1];
              momentum(2) = appliedB[2];

                
              //density
              rhoH = materialPointers[i]->getRho() ;

              //multiply acceleration by density to form momentum
              momentum *= rhoH ;


              //residual and tangent calculations node loops
              for ( j=0, jj=0; j<numberNodes; j++, jj+=ndf ) {

                temp = shp[massIndex][j] * ddvol ;

                for ( p = 0; p < 3; p++ )
                  resid( jj+p ) += ( temp * momentum(p) ) ;
              }
          }
      }

	return;

}


//************************************************************************
//compute local coordinates and basis

void   
ShellDKGT::computeBasis( ) 
{  
  //could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 } 
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  //and use those as basis vectors but this is easier 
  //and the shell is flat anyway.

  static Vector temp(3) ;

  static Vector v1(3) ;
  static Vector v2(3) ;
  static Vector v3(3) ;

  //get two vectors (v1, v2) in plane of shell by 
  // nodal coordinate differences

  const Vector &coor0 = nodePointers[0]->getCrds( ) ;

  const Vector &coor1 = nodePointers[1]->getCrds( ) ;

  const Vector &coor2 = nodePointers[2]->getCrds( ) ;
 

  v1.Zero( ) ;
  v1  = coor1;
  v1 -= coor0;
 
  
  v2.Zero( ) ;
  v2  = coor2 ;
  v2 -= coor0 ;
 
  //normalize v1 
  //double length = LovelyNorm( v1 ) ;
  double length = v1.Norm( ) ;
  v1 /= length ;

  //Gram-Schmidt process for v2 

  //double alpha = LovelyInnerProduct( v2, v1 ) ;
  double alpha = v2^v1 ;

  //v2 -= alpha*v1 ;
  temp = v1 ;
  temp *= alpha ;
  v2 -= temp ;

  //normalize v2 
  //length = LovelyNorm( v2 ) ;
  length = v2.Norm( ) ;
  v2 /= length ;

  //cross product for v3  
  v3 = LovelyCrossProduct( v1, v2 ) ;
  
  //local nodal coordinates in plane of shell

  int i ;
  for ( i = 0; i < 3; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( );
       xl[0][i] = coorI^v1 ;  
       xl[1][i] = coorI^v2 ;

  }  //end for i 

  //basis vectors stored as array of doubles
  for ( i = 0; i < 3; i++ ) {
      g1[i] = v1(i) ;
      g2[i] = v2(i) ;
      g3[i] = v3(i) ;
  }  //end for i 


}

//start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//************************************************************************
//compute local coordinates and basis

void   
ShellDKGT::updateBasis( ) 
{

  //could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 } 
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  //and use those as basis vectors but this is easier 
  //and the shell is flat anyway.

  static Vector temp(3) ;

  static Vector v1(3) ;
  static Vector v2(3) ;
  static Vector v3(3) ;

  //get two vectors (v1, v2) 
  // min plane of shell by 
  // nodal coordinate differences

  const Vector &coor0 = nodePointers[0]->getCrds( ) + nodePointers[0]->getTrialDisp();

  const Vector &coor1 = nodePointers[1]->getCrds( ) + nodePointers[1]->getTrialDisp();

  const Vector &coor2 = nodePointers[2]->getCrds( ) + nodePointers[2]->getTrialDisp();
  

 v1.Zero( ) ;
  //v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor1;
  v1 -= coor0 ;
 
  
  v2.Zero( ) ;
  //v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor2 ;
  v2 -= coor0 ;


  //normalize v1 
  //double length = LovelyNorm( v1 ) ;
  double length = v1.Norm( ) ;

  v1 /= length ;

  //Gram-Schmidt process for v2 

  //double alpha = LovelyInnerProduct( v2, v1 ) ;
  double alpha = v2^v1 ;

  //v2 -= alpha*v1 ;
  temp = v1 ;
  temp *= alpha ;
  v2 -= temp ;

  //normalize v2 
  //length = LovelyNorm( v2 ) ;
  length = v2.Norm( ) ;
  v2 /= length ;

  //cross product for v3  
  v3 = LovelyCrossProduct( v1, v2 ) ;
  
  //local nodal coordinates in plane of shell

  int i ;
  for ( i = 0; i < 3; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) + nodePointers[i]->getDisp( );//modify by Lisha Wang
       xl[0][i] = coorI^v1 ;  
       xl[1][i] = coorI^v2 ;

  }  //end for i 

  //basis vectors stored as array of doubles
  for ( i = 0; i < 3; i++ ) {
      g1[i] = v1(i) ;
      g2[i] = v2(i) ;
      g3[i] = v3(i) ;
  }  //end for i 


}
//end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)




//**********************************************************************

int  ShellDKGT::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = this->getTag();
  idData(9) = connectedExternalNodes(0);
  idData(10) = connectedExternalNodes(1);
  idData(11) = connectedExternalNodes(2);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellDKGT::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(4);               //????
  //vectData(0) = Ktt;
  vectData(0) = alphaM;
  vectData(1) = betaK;
  vectData(2) = betaK0;
  vectData(3) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellDKGT::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellDKGT::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}
    

int  ShellDKGT::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)      //idData����ı����ʲô��˼
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellDKGT::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(8));
  connectedExternalNodes(0) = idData(9);
  connectedExternalNodes(1) = idData(10);
  connectedExternalNodes(2) = idData(11);


  static Vector vectData(4);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellDKGT::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  //Ktt = vectData(0);
  alphaM = vectData(0);
  betaK = vectData(1);
  betaK0 = vectData(2);
  betaKc = vectData(3);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "ShellDKGT::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ShellDKGT::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewSection(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "ShellDKGT::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ShellDKGT::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}
//**************************************************************************

int
ShellDKGT::displaySelf(Renderer &theViewer, int displayMode, float fact)
{

    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();	
    const Vector &end3Crd = nodePointers[2]->getCrds();	

    static Matrix coords(3,3);
    static Vector values(4);
    static Vector P(32) ;

    for (int j=0; j<4; j++)
		values(j) = 0.0;

    if (displayMode >= 0) {
		// Display mode is positive:
		// display mode = 0 -> plot no contour
		// display mode = 1-8 -> plot 1-8 stress resultant

		// Get nodal displacements
		const Vector &end1Disp = nodePointers[0]->getDisp();
		const Vector &end2Disp = nodePointers[1]->getDisp();
		const Vector &end3Disp = nodePointers[2]->getDisp();

		// Get stress resultants
        if (displayMode <= 8 && displayMode > 0) {
			for (int i=0; i<4; i++) {
				const Vector &stress = materialPointers[i]->getStressResultant();
				values(i) = stress(displayMode-1);
			}
		}

		// Get nodal absolute position = OriginalPosition + (Displacement*factor)
		for (int i = 0; i < 3; i++) {
			coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
			coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
			coords(2,i) = end3Crd(i) + end3Disp(i)*fact;
		}
	} else {
		// Display mode is negative.
		// Plot eigenvectors
		int mode = displayMode * -1;
		const Matrix &eigen1 = nodePointers[0]->getEigenvectors();
		const Matrix &eigen2 = nodePointers[1]->getEigenvectors();
		const Matrix &eigen3 = nodePointers[2]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 3; i++) {
				coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
				coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
				coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;

			}    
		} else {
			for (int i = 0; i < 3; i++) {
				coords(0,i) = end1Crd(i);
				coords(1,i) = end2Crd(i);
				coords(2,i) = end3Crd(i);
				}
		}
	}

    int error = 0;
	
	// Draw a poligon with coordinates coords and values (colors) corresponding to values vector
    error += theViewer.drawPolygon (coords, values);

    return error;

}

void  
ShellDKGT::shape2d( double ss, double tt,double qq, 
		           const double x[2][3], 
		           double shp[3][3],
				   double &xsj ,double sx[2][2]
		           )

{ 

   int i ;

  double a[3],b[3],c[3];

  double area;
  double xs[2][2];
  
  a[0]=x[0][1]*x[1][2]-x[0][2]*x[1][1];
  a[1]=x[0][2]*x[1][0]-x[0][0]*x[1][2];
  a[2]=x[0][0]*x[1][1]-x[0][1]*x[1][0];

  b[0]=x[1][1]-x[1][2];
  b[1]=x[1][2]-x[1][0];
  b[2]=x[1][0]-x[1][1];

  c[0]=-x[0][1]+x[0][2];
  c[1]=-x[0][2]+x[0][0];
  c[2]=-x[0][0]+x[0][1];

  area=0.5*(x[0][0]*x[1][1]+x[0][1]*x[1][2]+x[0][2]*x[1][0]-x[0][0]*x[1][2]-x[0][1]*x[1][0]-x[0][2]*x[1][1]);

    shp[2][0]=ss;
	shp[2][1]=tt; 
	shp[2][2]=qq;

   xs[0][0]=x[0][1]-x[0][0];
   xs[0][1]=x[1][1]-x[1][0];
   xs[1][0]=x[0][2]-x[0][0];
   xs[1][1]=x[1][2]-x[1][0];

  xsj = (xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0]) ;
  
   //inverse jacobian
  double jinv ;
  jinv= 1/xsj ;
  sx[0][0] =  xs[1][1] * jinv ;
  sx[1][1] =  xs[0][0] * jinv ;
  sx[0][1] = -xs[0][1] * jinv ; 
  sx[1][0] = -xs[1][0] * jinv ; 

  for (i=0;i<3;i++)
  {
	
	  shp[0][i]=b[i]/2.0/area;
	  shp[1][i]=c[i]/2.0/area;
  }

  

  return ;
}



void
ShellDKGT::shapeDrill(double ss, double tt,double qq, 
					  const double x[2][3],
					  double sx[2][2], double shpDrill[4][3])
{
	


  double a[3],b[3],c[3];
  double l[3];
  double area;
  
  a[0]=x[0][1]*x[1][2]-x[0][2]*x[1][1];
  a[1]=x[0][2]*x[1][0]-x[0][0]*x[1][2];
  a[2]=x[0][0]*x[1][1]-x[0][1]*x[1][0];

  b[0]=x[1][1]-x[1][2];
  b[1]=x[1][2]-x[1][0];
  b[2]=x[1][0]-x[1][1];

  c[0]=-x[0][1]+x[0][2];
  c[1]=-x[0][2]+x[0][0];
  c[2]=-x[0][0]+x[0][1];

  area=0.5*(x[0][0]*x[1][1]+x[0][1]*x[1][2]+x[0][2]*x[1][0]-x[0][0]*x[1][2]-x[0][1]*x[1][0]-x[0][2]*x[1][1]);

 l[0]=ss;
 l[1]=tt;
 l[2]=qq;

	  shpDrill[0][0]=b[0]*(b[2]*l[1]-b[1]*l[2])/4.0/area; 
	  shpDrill[1][0]=c[0]*(b[2]*l[1]-b[1]*l[2])/4.0/area;
	  shpDrill[2][0]=b[0]*(c[2]*l[1]-c[1]*l[2])/4.0/area;
	  shpDrill[3][0]=c[0]*(c[2]*l[1]-c[1]*l[2])/4.0/area;

	  shpDrill[0][1]=b[1]*(b[0]*l[2]-b[2]*l[0])/4.0/area; 
	  shpDrill[1][1]=c[1]*(b[0]*l[2]-b[2]*l[0])/4.0/area;
	  shpDrill[2][1]=b[1]*(c[0]*l[2]-c[2]*l[0])/4.0/area;
	  shpDrill[3][1]=c[1]*(c[0]*l[2]-c[2]*l[0])/4.0/area;

	  shpDrill[0][2]=b[2]*(b[1]*l[0]-b[0]*l[1])/4.0/area; 
	  shpDrill[1][2]=c[2]*(b[1]*l[0]-b[0]*l[1])/4.0/area;
	  shpDrill[2][2]=b[2]*(c[1]*l[0]-c[0]*l[1])/4.0/area;
	  shpDrill[3][2]=c[2]*(c[1]*l[0]-c[0]*l[1])/4.0/area;
	  
	 
	return;
}


const Matrix&  
ShellDKGT::assembleB( const Matrix &Bmembrane,
                               const Matrix &Bbend, 
                               const Matrix &Bshear ) 
{
	static Matrix B(8,6);
	
	int p,q;

	int pp;

// For Shell : 
//
//---B Matrices in standard {1,2,3} mechanics notation---------
//
//            -                     _          
//           | Bmembrane  |     0    |
//           | --------------------- |     
//    B =    |     0      |  Bbend   |   (8x6) 
//           | --------------------- |
//           |            |  Bshear  |
//            -           -         -
//
//-------------------------------------------------------------
	B.Zero( );
	//assemble B from sub-matrices

	//membrane parts
	for(p=0; p<3; p++) {

		for(q=0; q<3; q++)
			B(p,q) = Bmembrane(p,q);
	}//end for p

	//bending parts
	for (p=3; p<6; p++){
		for(q=3; q<6;q++)
			B(p,q)=Bbend(p-3,q-3);
	}//end for p

	//shear parts
	for(p=0;p<2;p++){
       
		for(q=3; q<6; q++)
			B(p+6,q)=Bshear(p,q-3);
	}//end for p

   return B;

}

//***********************************************************************
//compute Bmembrane matrix

const Matrix&
ShellDKGT::computeBmembrane(int node, const double shp[3][3],
							const double shpDrill[4][3])
{
	static Matrix Bmembrane(3,3);

	// ------Bmembrane Matrix in standard {1,2,3} mechanics notation ---------------
	// 
	//               | N,1    0      Nu,1  |
	//   Bmembrane = | 0     N,2     Nv,2  |
	//               | N,2   N,1  Nv,1+Nu,2|
   // -----------------------------------------------------------------------------------

	Bmembrane.Zero( );

	Bmembrane(0,0) = shp[0][node];
	Bmembrane(0,2) = shpDrill[0][node];
	Bmembrane(1,1) = shp[1][node]; 
	Bmembrane(1,2) = shpDrill[3][node];
	Bmembrane(2,0) = shp[1][node];
	Bmembrane(2,1) = shp[0][node];
	Bmembrane(2,2) = shpDrill[1][node]+shpDrill[2][node];

	return Bmembrane;

}


//***********************************************************************
//compute Bbend matrix

const Matrix&
 ShellDKGT::computeBbend(int node ,const double shpBend[6][9])
 {
	 static Matrix Bbend(3,3);

	 int i,j,k;

	 //----------Bbend Matrix in standard {1,2,3}mechanics notation------------------
	 //
	 //               |
	 //      Bbend =  |  Hx[3*node],1  Hx[3*node + 1],1  Hx[3*node + 2],1  |
	 //               |  Hy[3*node],2  Hy[3*node + 1],2  Hy[3*node + 2],2  |
	 //               |  Hx[3*node],2  Hx[3*node + 1],2  Hx[3*node + 2],2  |
	 //               |  + Hy[3*node],1 +Hy[3*node + 1],1 +Hy[3*node + 2],1|
	 //------------------------------------------------------------------------------
	 i = 3*node;
	 j = 3*node + 1;
	 k = 3*node + 2;

	 Bbend.Zero( );

	 Bbend(0,0) = shpBend[2][i];
	 Bbend(0,1) = shpBend[2][j];
	 Bbend(0,2) = shpBend[2][k];

	 Bbend(1,0) = shpBend[5][i];
	 Bbend(1,1) = shpBend[5][j];
	 Bbend(1,2) = shpBend[5][k];

 	 Bbend(2,0) = shpBend[3][i] + shpBend[4][i];
	 Bbend(2,1) = shpBend[3][j] + shpBend[4][j];
	 Bbend(2,2) = shpBend[3][k] + shpBend[4][k];

	 return Bbend;
 }


//*********************************************************************
//shape function for bending plate


void 
	ShellDKGT::shapeBend(double ss, double tt, double qq, const double x[2][3],
					 double sx[2][2], double shpBend[6][9])
{
	static double N[3][6];
	static double temp[4][9];

	double a4, a5, a6;
	double b4, b5, b6;
	double c4, c5, c6;
	double d4, d5, d6;
	double e4, e5, e6;


	double p4,p5,p6,t4,t5,t6,q4,q5,q6,r4,r5,r6;

	double x12, x23, x31;
	double y12, y23, y31;
	double L12, L23, L31;

	double a[3],b[3],c[3];
    double l[3];
    double area;

	const double one_over_four = 1.0/4.0;

	//define xij,yij,lij
	x12 = x[0][0] - x[0][1];
	x23 = x[0][1] - x[0][2];
	x31 = x[0][2] - x[0][0];


	y12 = x[1][0] - x[1][1];
	y23 = x[1][1] - x[1][2];
	y31 = x[1][2] - x[1][0];


	L12 = x12 * x12 + y12 * y12;
	L23 = x23 * x23 + y23 * y23;
	L31 = x31 * x31 + y31 * y31;

    a[0]=x[0][1]*x[1][2]-x[0][2]*x[1][1];
    a[1]=x[0][2]*x[1][0]-x[0][0]*x[1][2];
    a[2]=x[0][0]*x[1][1]-x[0][1]*x[1][0];

    b[0]=x[1][1]-x[1][2];
    b[1]=x[1][2]-x[1][0];
    b[2]=x[1][0]-x[1][1];

    c[0]=-x[0][1]+x[0][2];
    c[1]=-x[0][2]+x[0][0];
    c[2]=-x[0][0]+x[0][1];

    area=0.5*(x[0][0]*x[1][1]+x[0][1]*x[1][2]+x[0][2]*x[1][0]-x[0][0]*x[1][2]-x[0][1]*x[1][0]-x[0][2]*x[1][1]);

    l[0]=ss;
    l[1]=tt;
    l[2]=qq;

	//a5-8,b5-8,c5-8,d5-8,e5-8
	a4 = - x23/L23;
	b4 = one_over_four * 3.0 * x23 * y23/L23;
	c4 = one_over_four *(x23*x23 - 2.0*y23*y23)/L23;
	d4 = - y23/L23;
	e4 = one_over_four *(y23*y23 - 2.0*x23*x23)/L23;

	a5 = - x31/L31;
	b5 = one_over_four * 3.0 * x31 * y31/L31;
	c5 = one_over_four *(x31*x31 - 2.0*y31*y31)/L31;
	d5 = - y31/L31;
	e5 = one_over_four *(y31*y31 - 2.0*x31*x31)/L31;

	a6 = - x12/L12;
	b6 = one_over_four * 3.0 * x12 * y12/L12;
	c6 = one_over_four *(x12*x12 - 2.0*y12*y12)/L12;
	d6 = - y12/L12;
	e6 = one_over_four *(y12*y12 - 2.0*x12*x12)/L12;

	p4=6*a4; p5=6*a5; p6=6*a6;
	t4=6*d4; t5=6*d5; t6=6*d6;
	q4=4*b4; q5=4*b5; q6=4*b6;

	r4=3*y23*y23/L23;
	r5=3*y31*y31/L31;
	r6=3*y12*y12/L12;

	//define 3-d isoparametric shape function
	N[0][0] = (2.0*l[0]-1.0)*l[0];
	N[0][1] = (2.0*l[1]-1.0)*l[1];
	N[0][2] = (2.0*l[2]-1.0)*l[2];
	N[0][3] = 4.0*l[2]*l[1];
	N[0][4] = 4.0*l[0]*l[2];
	N[0][5] = 4.0*l[1]*l[0];

	//Hx, Hy
	const double three_over_two = 3.0/2.0;
	//Hx
	shpBend[0][0] = (a6*N[0][5]-a5*N[0][4]) * three_over_two;
	shpBend[0][1] = b5*N[0][4] + b6*N[0][5];
	shpBend[0][2] = N[0][0] - c5*N[0][4] - c6*N[0][5];

	shpBend[0][3] = (a4*N[0][3]-a6*N[0][5]) * three_over_two;
	shpBend[0][4] = b6*N[0][5] + b4*N[0][3];
	shpBend[0][5] = N[0][1] - c6*N[0][5] - c4*N[0][3];

	shpBend[0][6] = (a5*N[0][4]-a4*N[0][3]) * three_over_two;
	shpBend[0][7] = b4*N[0][3] + b5*N[0][4];
	shpBend[0][8] = N[0][2] - c4*N[0][3] - c5*N[0][4];


	//Hy
	shpBend[1][0] = (d6*N[0][5]-d5*N[0][4])* three_over_two;
	shpBend[1][1] = -N[0][0] + e5*N[0][4] + e6*N[0][5];
	shpBend[1][2] = -b5*N[0][4] - b6*N[0][5];

	shpBend[1][3] = (d4*N[0][3]-d6*N[0][5])* three_over_two;
	shpBend[1][4] = -N[0][1] + e6*N[0][5] + e4*N[0][3];
	shpBend[1][5] = -b6*N[0][5] - b4*N[0][3];

	shpBend[1][6] = (d5*N[0][4]-d4*N[0][3])* three_over_two;
	shpBend[1][7] = -N[0][2] + e4*N[0][3] + e5*N[0][4];
	shpBend[1][8] = -b4*N[0][3] - b5*N[0][4];

	//Hx,xi

	temp[0][0]=p6*(1.0-2.0*l[1])+(p5-p6)*l[2];
	temp[0][1]=q6*(1.0-2.0*l[1])-(q5+q6)*l[2];
	temp[0][2]=-4.0+6.0*(l[1]+l[2])+r6*(1.0-2.0*l[1])-l[2]*(r5+r6);
	temp[0][3]=-p6*(1.0-2.0*l[1])+l[2]*(p4+p6);
	temp[0][4]=q6*(1.0-2.0*l[1])-l[2]*(q6-q4);
	temp[0][5]=-2.0+6.0*l[1]+r6*(1.0-2.0*l[1])+l[2]*(r4-r6);
	temp[0][6]=-l[2]*(p5+p4);
	temp[0][7]=l[2]*(q4-q5);
	temp[0][8]=-l[2]*(r5-r4);

	//Hx,eta
	temp[1][0]=-p5*(1.0-2.0*l[2])-l[1]*(p6-p5);
	temp[1][1]=q5*(1.0-2.0*l[2])-l[1]*(q5+q6);
	temp[1][2]=-4.0+6.0*(l[1]+l[2])+r5*(1.0-2.0*l[2])-l[1]*(r5+r6);
	temp[1][3]=l[1]*(p4+p6);
	temp[1][4]=l[1]*(q4-q6);
	temp[1][5]=-l[1]*(r6-r4);
	temp[1][6]=p5*(1.0-2.0*l[2])-l[1]*(p4+p5);
	temp[1][7]=q5*(1.0-2.0*l[2])+l[1]*(q4-q5);
	temp[1][8]=-2.0+6.0*l[2]+r5*(1.0-2.0*l[2])+l[1]*(r4-r5);


	//Hy,xi
	temp[2][0]=t6*(1.0-2.0*l[1])+l[2]*(t5-t6);
	temp[2][1]=1.0+r6*(1.0-2.0*l[1])-l[2]*(r5+r6);
	temp[2][2]=-q6*(1.0-2.0*l[1])+l[2]*(q5+q6);
	temp[2][3]=-t6*(1.0-2.0*l[1])+l[2]*(t4+t6);
	temp[2][4]=-1.0+r6*(1.0-2.0*l[1])+l[2]*(r4-r6);
	temp[2][5]=-q6*(1.0-2.0*l[1])-l[2]*(q4-q6);
	temp[2][6]=-l[2]*(t4+t5);
	temp[2][7]=l[2]*(r4-r5);
	temp[2][8]=-l[2]*(q4-q5);

	//Hy,eta
	temp[3][0]=-t5*(1.0-2.0*l[2])-l[1]*(t6-t5);
	temp[3][1]=1+r5*(1.0-2.0*l[2])-l[1]*(r5+r6);
	temp[3][2]=-q5*(1.0-2.0*l[2])+l[1]*(q5+q6);
	temp[3][3]=l[1]*(t4+t6);
	temp[3][4]=l[1]*(r4-r6);
	temp[3][5]=-l[1]*(q4-q6);
	temp[3][6]=t5*(1.0-2.0*l[2])-l[1]*(t4+t5);
	temp[3][7]=-1.0+r5*(1.0-2.0*l[2])+l[1]*(r4-r5);
	temp[3][8]=-q5*(1.0-2.0*l[2])-l[1]*(q4-q5);

		for (int i=0;i<9;i++) {
		
		shpBend[2][i] = temp[0][i]*y31/2.0/area +temp[1][i]*y12/2.0/area; //Hx,x
		shpBend[3][i] = temp[0][i]*(-x31)/2.0/area +temp[1][i]*(-x12)/2.0/area; //Hx,y
		shpBend[4][i] = temp[2][i]*y31/2.0/area +temp[3][i]*y12/2.0/area; //Hy,x
		shpBend[5][i] = temp[2][i]*(-x31)/2.0/area +temp[3][i]*(-x12)/2.0/area; //Hy,y
	}//end for i

	return;

}
