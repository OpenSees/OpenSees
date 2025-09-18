 /* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                     1               **
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
// $Date: 2011/03/10 22:51:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShellMITC4Thermal.cpp,v $

// Written: Leopoldo Tesser, Diego A. Talledo, V�ronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984
// Modified for SIF modelling by Jian Jiang, Liming Jiang [http://openseesforfire.github.io] 


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellMITC4Thermal.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <NodalThermalAction.h>
#include <ThermalActionWrapper.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#define min(a,b) ( (a)<(b) ? (a):(b) )

static int numShellMITC4Thermal = 0;

void *
OPS_ShellMITC4Thermal(void)
{
  if (numShellMITC4Thermal == 0) {
    opserr << "Using ShellMITC4Thermal - Developed by: Leopoldo Tesser, Diego A. Talledo, Veronique Le Corvec\n";
    numShellMITC4Thermal++;
  }

  Element *theElement = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4Thermal $tag $iNode $jNoe $kNode $lNode $secTag";
    return 0;	
  }
  
  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC4Thermal \n";
    return 0;
  }

  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(iData[5]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC4Thermal " << iData[0] << "section " << iData[5] << " not found\n";
    return 0;
  }
  
  theElement = new ShellMITC4Thermal(iData[0], iData[1], iData[2], iData[3],
			      iData[4], *theSection);

  return theElement;

}


//static data
Matrix  ShellMITC4Thermal::stiff(24,24) ;
Vector  ShellMITC4Thermal::resid(24) ;
Matrix  ShellMITC4Thermal::mass(24,24) ;

//quadrature data
const double  ShellMITC4Thermal::root3 = sqrt(3.0) ;
const double  ShellMITC4Thermal::one_over_root3 = 1.0 / root3 ;

double ShellMITC4Thermal::sg[4] ;
double ShellMITC4Thermal::tg[4] ;
double ShellMITC4Thermal::wg[4] ;

 

//null constructor
ShellMITC4Thermal::ShellMITC4Thermal( ) :
Element( 0, ELE_TAG_ShellMITC4Thermal ),
connectedExternalNodes(4), load(0), Ki(0)
{ 
  for (int i = 0 ;  i < 4; i++ ) 
    materialPointers[i] = 0;

  sg[0] = -one_over_root3;
  sg[1] = one_over_root3;
  sg[2] = one_over_root3;
  sg[3] = -one_over_root3;  

  tg[0] = -one_over_root3;
  tg[1] = -one_over_root3;
  tg[2] = one_over_root3;
  tg[3] = one_over_root3;  

  wg[0] = 1.0;
  wg[1] = 1.0;
  wg[2] = 1.0;
  wg[3] = 1.0;
}


//*********************************************************************
//full constructor
ShellMITC4Thermal::ShellMITC4Thermal(  int tag, 
                         int node1,
                         int node2,
   	                     int node3,
                         int node4,
	                     SectionForceDeformation &theMaterial ) :
Element( tag, ELE_TAG_ShellMITC4Thermal ),
connectedExternalNodes(4), load(0), Ki(0)
{
  int i;

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy( ) ;

      if (materialPointers[i] == 0) {
	      opserr << "ShellMITC4Thermal::constructor - failed to get a material of type: ShellSection\n";
      } //end if
      
  } //end for i 

  sg[0] = -one_over_root3;
  sg[1] = one_over_root3;
  sg[2] = one_over_root3;
  sg[3] = -one_over_root3;  

  tg[0] = -one_over_root3;
  tg[1] = -one_over_root3;
  tg[2] = one_over_root3;
  tg[3] = one_over_root3;  

  wg[0] = 1.0;
  wg[1] = 1.0;
  wg[2] = 1.0;
  wg[3] = 1.0;
  
  //added for ThermalAction
    dataMix = new double [18];
   for ( int i=0; i<18; i++)
   dataMix[i] = 0.0;
   
  for ( int i=0; i<8; i++)
   residThermal[i] = 0.0;
  
  counterTemperature = 0;

 }
//******************************************************************

//destructor 
ShellMITC4Thermal::~ShellMITC4Thermal( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;

  } //end for i

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
}
//**************************************************************************


//set domain
void  ShellMITC4Thermal::setDomain( Domain *theDomain ) 
{  
  int i, j ;
  static Vector eig(3) ;
  static Matrix ddMembrane(3,3) ;

  //node pointers
  for ( i = 0; i < 4; i++ ) {
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
     if (nodePointers[i] == 0) {
       opserr << "ShellMITC4Thermal::setDomain - no node " << connectedExternalNodes(i);
       opserr << " exists in the model\n";
     }
     const Vector &nodeDisp=nodePointers[i]->getTrialDisp();
     if (nodeDisp.Size() != 6) {
       opserr << "ShellMITC4Thermal::setDomain - node " << connectedExternalNodes(i);
       opserr << " NEEDS 6 dof - GARBAGE RESULTS or SEGMENTATION FAULT WILL FOLLOW\n";
     }       
  }

  //compute drilling stiffness penalty parameter
  const Matrix &dd = materialPointers[0]->getInitialTangent( ) ;

  //assemble ddMembrane ;
  for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ )
         ddMembrane(i,j) = dd(i,j) ;
  } //end for i 

  //eigenvalues of ddMembrane
  eig = LovelyEig( ddMembrane ) ;
  
  //set ktt 
  //Ktt = dd(2,2) ;  //shear modulus 
  Ktt = min( eig(2), min( eig(0), eig(1) ) ) ;
  //Ktt = dd(2,2);

  //basis vectors and local coordinates
  computeBasis( ) ;

  this->DomainComponent::setDomain(theDomain);

}


//get the number of external nodes
int  ShellMITC4Thermal::getNumExternalNodes( ) const
{
  return 4 ;
} 
 

//return connected external nodes
const ID&  ShellMITC4Thermal::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


Node **
ShellMITC4Thermal::getNodePtrs(void) 
{
  return nodePointers;
} 

//return number of dofs
int  ShellMITC4Thermal::getNumDOF( ) 
{
  return 24 ;
}


//commit state
int  ShellMITC4Thermal::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ShellMITC4Thermal::commitState () - failed in base class";
  }    

  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;

  return success ;
}
 


//revert to last commit 
int  ShellMITC4Thermal::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  ShellMITC4Thermal::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  ShellMITC4Thermal::Print(OPS_Stream &s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ShellMITC4Thermal\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
        s << "\t" << connectedExternalNodes(2) << "\t" << connectedExternalNodes(3) << "\t0.00";
        s << endln;
        s << "PROP_3D\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << -1 << "\tSHELL\t1.0\0.0";
        s << endln;
    }
    
    if (flag < -1) {
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
        s << "MITC4 Non-Locking Four Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << connectedExternalNodes(0) << endln;
        s << "Node 2 : " << connectedExternalNodes(1) << endln;
        s << "Node 3 : " << connectedExternalNodes(2) << endln;
        s << "Node 4 : " << connectedExternalNodes(3) << endln;

        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);

        s << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ShellMITC4Thermal\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << ", " << connectedExternalNodes(3) << "], ";
        s << "\"section\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}

Response*
ShellMITC4Thermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "ShellMITC4Thermal");
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
      opserr << "ShellMITC4Thermal::setResponse() - need to specify more data\n";
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
ShellMITC4Thermal::getResponse(int responseID, Information &eleInfo)
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
}


//return stiffness matrix 
const Matrix&  ShellMITC4Thermal::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;
  
#ifdef _SDEBUG
  if(this->getTag()==1)
    opserr<< "ShellMITC4Thermal: "<<this->getTag()<< "Tangent stiffness: "<<endln<<stiff<<endln;

#endif
  
  return stiff ;
}    

//return secant matrix 
const Matrix&  ShellMITC4Thermal::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  static const int nstress = 8 ; //three membrane, three moment, two shear

  static const int ngauss = 4 ;

  static const int numnodes = 4 ;

  int i,  j,  k, p, q ;
  int jj, kk ;

  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  //  static double Shape[3][numnodes][ngauss] ; //all the shape functions

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Matrix dd(nstress,nstress) ;  //material tangent

  static Matrix J0(2,2) ;  //Jacobian at center
 
  static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;


    static Matrix Bbend(3,3) ;  // bending B matrix

    static Matrix Bshear(2,3) ; // shear B matrix

    static Matrix Bmembrane(3,2) ; // membrane B matrix


    static double BdrillJ[ndf] ; //drill B matrix

    static double BdrillK[ndf] ;  

    double *drillPointer ;

    static double saveB[nstress][ndf][numnodes] ;

  //-------------------------------------------------------

  stiff.Zero( ) ;

 
  double dx34 = xl[0][2]-xl[0][3];
  double dy34 = xl[1][2]-xl[1][3];

  double dx21 = xl[0][1]-xl[0][0];
  double dy21 = xl[1][1]-xl[1][0];

  double dx32 = xl[0][2]-xl[0][1];
  double dy32 = xl[1][2]-xl[1][1];

  double dx41 = xl[0][3]-xl[0][0];
  double dy41 = xl[1][3]-xl[1][0];

  Matrix G(4,12);
  G.Zero();
  double one_over_four = 0.25;
  G(0,0)=-0.5;
  G(0,1)=-dy41*one_over_four;
  G(0,2)=dx41*one_over_four;
  G(0,9)=0.5;
  G(0,10)=-dy41*one_over_four;
  G(0,11)=dx41*one_over_four;
  G(1,0)=-0.5;
  G(1,1)=-dy21*one_over_four;
  G(1,2)=dx21*one_over_four;
  G(1,3)=0.5;
  G(1,4)=-dy21*one_over_four;
  G(1,5)=dx21*one_over_four;
  G(2,3)=-0.5;
  G(2,4)=-dy32*one_over_four;
  G(2,5)=dx32*one_over_four;
  G(2,6)=0.5;
  G(2,7)=-dy32*one_over_four;
  G(2,8)=dx32*one_over_four;
  G(3,6)=0.5;
  G(3,7)=-dy34*one_over_four;
  G(3,8)=dx34*one_over_four;
  G(3,9)=-0.5;
  G(3,10)=-dy34*one_over_four;
  G(3,11)=dx34*one_over_four;

  Matrix Ms(2,4);
  Ms.Zero();
  Matrix Bsv(2,12);
  Bsv.Zero();

  double Ax = -xl[0][0]+xl[0][1]+xl[0][2]-xl[0][3];
  double Bx =  xl[0][0]-xl[0][1]+xl[0][2]-xl[0][3];
  double Cx = -xl[0][0]-xl[0][1]+xl[0][2]+xl[0][3];

  double Ay = -xl[1][0]+xl[1][1]+xl[1][2]-xl[1][3];
  double By =  xl[1][0]-xl[1][1]+xl[1][2]-xl[1][3];
  double Cy = -xl[1][0]-xl[1][1]+xl[1][2]+xl[1][3];

  double alph = atan(Ay/Ax);
  double beta = 3.141592653589793/2-atan(Cx/Cy);
  Matrix Rot(2,2);
  Rot.Zero();
  Rot(0,0)=sin(beta);
  Rot(0,1)=-sin(alph);
  Rot(1,0)=-cos(beta);
  Rot(1,1)=cos(alph);
  Matrix Bs(2,12);
  
  double r1 = 0;
  double r2 = 0;
  double r3 = 0;

  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {
    
    r1 = Cx + sg[i]*Bx;
    r3 = Cy + sg[i]*By;
    r1 = r1*r1 + r3*r3;
	r1 = sqrt (r1);
	r2 = Ax + tg[i]*Bx;
	r3 = Ay + tg[i]*By;
	r2 = r2*r2 + r3*r3;
	r2 = sqrt (r2);

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;
    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  
    volume += dvol[i] ;

    Ms(1,0)=1-sg[i];
	Ms(0,1)=1-tg[i];
	Ms(1,2)=1+sg[i];
	Ms(0,3)=1+tg[i];
	Bsv = Ms*G;

    for ( j = 0; j < 12; j++ ) {
		Bsv(0,j)=Bsv(0,j)*r1/(8*xsj);
		Bsv(1,j)=Bsv(1,j)*r2/(8*xsj);
    }
    Bs=Rot*Bsv;
    
    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

      for ( p = 0; p < 3; p++) {
		  Bshear(0,p) = Bs(0,j*3+p);
		  Bshear(1,p) = Bs(1,j*3+p);
      }//end for p

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //save the B-matrix
      for (p=0; p<nstress; p++) {
	    for (q=0; q<ndf; q++ )
	      saveB[p][q][j] = BJ(p,q) ;
      }//end for p

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	    //BdrillJ[p] = *drillPointer++ ;
	    BdrillJ[p] = *drillPointer ; //set p-th component
	    drillPointer++ ;             //pointer arithmetic
      }//end for p
    } // end for j
  

    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;

    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numnodes; j++ ) {

      //extract BJ
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  BJ(p,q) = saveB[p][q][j]   ;
      }//end for p

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      for ( p = 3; p < 6; p++ ) {
	for ( q = 3; q < 6; q++ ) 
	  BJ(p,q) *= (-1.0) ;
      } //end for p


      //transpose 
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	BdrillJ[p] = *drillPointer ;
	drillPointer++ ;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0 ) ;
      
      for (p=0; p<ndf; p++) 
	BdrillJ[p] *= ( Ktt*dvol[i] ) ;
      
      kk = 0 ;
      for ( k = 0; k < numnodes; k++ ) {

	//extract BK
	for (p=0; p<nstress; p++) {
	  for (q=0; q<ndf; q++ )
	    BK(p,q) = saveB[p][q][k]   ;
	}//end for p
	
	
	//drilling B matrix
	drillPointer = computeBdrill( k, shp ) ;
	for (p=0; p<ndf; p++ ) {
	  BdrillK[p] = *drillPointer ;
	  drillPointer++ ;
	}//end for p
	
	//stiffJK = BJtranD * BK  ;
	// +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
	stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0 ) ;
	
	for ( p = 0; p < ndf; p++ )  {
	  for ( q = 0; q < ndf; q++ ) {
	    stiff( jj+p, kk+q ) += stiffJK(p,q) 
	      + ( BdrillJ[p]*BdrillK[q] ) ;
	  }//end for q
	}//end for p
	
	kk += ndf ;
      } // end for k loop
      
      jj += ndf ;
    } // end for j loop
    
  } //end for i gauss loop 

  Ki = new Matrix(stiff);
  
  
#ifdef _SDEBUG
  opserr<< "ShellMITC4Thermal: "<<this->getTag()<< "InitialTangent stiffness: "<<endln
  <<stiff<<endln;
#endif
  
  return stiff ;
}
    

//return mass matrix
const Matrix&  ShellMITC4Thermal::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


void  ShellMITC4Thermal::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}


int 
ShellMITC4Thermal::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
    const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_ShellThermalAction) {
  // ----- Real time Temperature is obtained from shellThermalAction 
	  /*
    dataMix[0] = data(0)*loadFactor;
    dataMix[2] = data(2)*loadFactor;
    dataMix[4] = data(4)*loadFactor;
    dataMix[6] = data(6)*loadFactor;
    dataMix[8] = data(8)*loadFactor;
    dataMix[10] = data(10)*loadFactor;
    dataMix[12] = data(12)*loadFactor;
    dataMix[14] = data(14)*loadFactor;
    dataMix[16] = data(16)*loadFactor;

    dataMix[1] = data(1);
    dataMix[3] = data(3);
    dataMix[5] = data(5);
    dataMix[7] = data(7);
    dataMix[9] = data(9);
    dataMix[11] = data(11);
    dataMix[13] = data(13);
    dataMix[15] = data(15);
    dataMix[17] = data(17);
    */

	Vector dataMixV=data;
	counterTemperature = 1;
    
	//opserr<<"ShellThermal :T data: "<<dataMixV<<endln;

  // Loop over the integration points
  
      
	for (int j=0; j<4; j++)  {
		const Vector &s = materialPointers[j]->getTemperatureStress(dataMixV);
		residThermal[2*j] = s(0);
		residThermal[2*j+1] = s(1);
	}//end for
/*
  const Vector &coor0 = nodePointers[0]->getCrds( ) ;

  const Vector &coor1 = nodePointers[1]->getCrds( ) ;

  const Vector &coor2 = nodePointers[2]->getCrds( ) ;
  
  const Vector &coor3 = nodePointers[3]->getCrds( ) ;
  // J.Jiang add to see coordinate of node
  double coorNode1[3],coorNode2[3],coorNode3[3],coorNode4[3];
 
 
  double L12,L23;
  L12 = fabs(coor1[0]-coor0[0])/2;
  L23 = fabs(coor2[1]-coor1[1])/2;

  //L12 = fabs(coor1[0]-coor0[0]);
  //L23 = fabs(coor2[1]-coor1[1]);

  //J.Jiang add to use resid per length
  //L12 = 1.0;
  //L23 = 1.0;
   
  
 //J.Jiang third try
  residThermal[0] = -ThermalLoad[0]*L23; residThermal[1] = -ThermalLoad[0]*L12;residThermal[3] = ThermalLoad[1]*L23;residThermal[4] = -ThermalLoad[1]*L12;
  residThermal[6] = ThermalLoad[0]*L23; residThermal[7] = -ThermalLoad[0]*L12;residThermal[9] = ThermalLoad[1]*L23;residThermal[10] = ThermalLoad[1]*L12;
  residThermal[12] = ThermalLoad[0]*L23; residThermal[13] = ThermalLoad[0]*L12;residThermal[15] = -ThermalLoad[1]*L23;residThermal[16] = ThermalLoad[1]*L12;
  residThermal[18] = -ThermalLoad[0]*L23; residThermal[19] = ThermalLoad[0]*L12;residThermal[21] = -ThermalLoad[1]*L23;residThermal[22] = -ThermalLoad[1]*L12;
 */
  }
  else if (type == LOAD_TAG_NodalThermalAction) {
 
	  //NodalLoad* theNodalThermal0,theNodalThermal1;	 
	 NodalThermalAction* theNodalThermal0 = nodePointers[0]->getNodalThermalActionPtr();
	 NodalThermalAction* theNodalThermal1 = nodePointers[1]->getNodalThermalActionPtr();	
	 NodalThermalAction* theNodalThermal2 = nodePointers[2]->getNodalThermalActionPtr();
	 NodalThermalAction* theNodalThermal3 = nodePointers[3]->getNodalThermalActionPtr();	


	 int type;
	 const Vector &data0 = theNodalThermal0->getData(type);
	 const Vector &data1 = theNodalThermal1->getData(type);
	 const Vector &data2 = theNodalThermal2->getData(type);
	 const Vector &data3 = theNodalThermal3->getData(type);

	 Vector Loc(9);
	 Vector NodalT0(9);
	 Vector NodalT1(9);
	 Vector NodalT2(9);
	 Vector NodalT3(9);

	 #ifdef _SHELLDBG
 	 opserr<<"NodalT0: "<<data0<<endln<<"NodalT1: "<<data1;
	 #endif

	 for(int i =0; i<9;i++){
		 if(data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
			 opserr<<"Warning:The NodalThermalAction in ShellMITC4Thermal "<<this->getTag()
			      << "incompatible loc input for datapoint "<< i << endln;
			 }
		 else{
			 Loc(i)=data0(2*i+1);
			 NodalT0(i)=data0(2*i);
			 NodalT1(i)=data1(2*i);
			 NodalT2(i)=data2(2*i);
			 NodalT3(i)=data3(2*i);
			 }
		 }
    
	double ThermalLoad[2];
	ThermalLoad[0] = 0.0;
    ThermalLoad[1] = 0.0;

	counterTemperature = 1;

  
  // Loop over the integration points
    
      
	for (int j=0; j<4; j++)  {
	
	 //Obtain interpolated temperature
		Vector dataMixV(18);
		double xi = sg[j];
		double eta = tg[j];
		
		for(int k=0; k<9 ;k++) {
			dataMixV(2*k)=  shapefn2d(xi, eta, 1)*NodalT0(k)+shapefn2d(xi, eta, 2)*NodalT1(k)+
				            shapefn2d(xi, eta, 3)*NodalT2(k)+shapefn2d(xi, eta, 4)*NodalT3(k);
			dataMixV(2*k+1) = Loc(k);

	    }

		const Vector &s = materialPointers[j]->getTemperatureStress(dataMixV);
		residThermal[2*j] = s(0);
		residThermal[2*j+1] = s(1);
      if (s == 0) {

	//opserr << "ShellMITC4Thermal::addLoad - failed to get a thermal load";
      } //end if
	}//end for

 }
  //end of NodalThermalAction
  else if (type == LOAD_TAG_ThermalActionWrapper) {
	counterTemperature = 1;
    // Loop over the integration points
     Vector theNode1Crds = nodePointers[0]->getCrds();
	 Vector theNode2Crds = nodePointers[1]->getCrds();
	 Vector theNode3Crds = nodePointers[2]->getCrds();
	 Vector theNode4Crds = nodePointers[3]->getCrds();

	 int ndm = theNode1Crds.Size();
	 Vector theIntCrds = Vector(ndm);

	for (int j=0; j<4; j++)  {
	
	 //Obtain interpolated temperature
		
		double xi = sg[j];
		double eta = tg[j];
		
		theIntCrds.Zero();
		for(int i=0;i<3;i++){
			theIntCrds(i) = shapefn2d(xi, eta, 1)*theNode1Crds(i)+shapefn2d(xi, eta, 2)*theNode2Crds(i)+
				            shapefn2d(xi, eta, 3)*theNode3Crds(i)+shapefn2d(xi, eta, 4)*theNode4Crds(i);
		}

#ifdef _SDEBUG
	  if(this->getTag()==1){
	  opserr<<"ShellMITC4Thermal::addLoad at the section "<<j<<" --theIntCrds"<<theIntCrds<<endln;
		  }
#endif

	Vector dataMixV = ((ThermalActionWrapper*) theLoad)->getIntData(theIntCrds);

     const Vector &s = materialPointers[j]->getTemperatureStress(dataMixV);
     residThermal[2*j] = s(0);
	 residThermal[2*j+1] = s(1);
    
	}//end of for loop

  }
  //end of thermalActionWrapper
  else {
    opserr << "ShellMITC4Thermal::ShellMITC4Thermal -- load type unknown for element with tag: "
	   << this->getTag() << "ShellMITC4Thermal::addLoad()\n"; 
			    
    return -1;
  }

  return 0;

}



int 
ShellMITC4Thermal::addInertiaLoadToUnbalance(const Vector &accel)
{
  int tangFlag = 1 ;
  static Vector r(24);

  int i;

  int allRhoZero = 0;
  for (i=0; i<4; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      allRhoZero = 1;
  }

  if (allRhoZero == 0) 
    return 0;

  formInertiaTerms( tangFlag ) ;

  int count = 0;
  for (i=0; i<4; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<6; j++)
      r(count++) = Raccel(j);
  }

  if (load == 0) 
    load = new Vector(24);

  load->addMatrixVector(1.0, mass, r, -1.0);

  return 0;
}



//get residual
const Vector&  ShellMITC4Thermal::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

#ifdef _SDEBUG
  opserr<< "ShellMITC4Thermal: "<<this->getTag()<< " Resid: "<<endln
	 <<resid<<endln;
#endif

  if (counterTemperature == 1)
{
	//for ( int i=0; i<24; i++)
//resid[i] +=residThermal[i];
}



  return resid ;   
}


//get residual with inertia terms
const Vector&  ShellMITC4Thermal::getResistingForceIncInertia( )
{
  static Vector res(24);
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
ShellMITC4Thermal::formInertiaTerms( int tangFlag ) 
{

  //translational mass only
  //rotational inertia terms are neglected


  static const int ndf = 6 ; 

  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

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
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //volume element to also be saved
    dvol = wg[i] * xsj ;  


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
int  
ShellMITC4Thermal::formResidAndTangent( int tang_flag ) 
{
  //
  //  six(6) nodal dof's ordered :
  //
  //    -        - 
  //   |    u1    |   <---plate membrane
  //   |    u2    |
  //   |----------|
  //   |  w = u3  |   <---plate bending
  //   |  theta1  | 
  //   |  theta2  | 
  //   |----------|
  //   |  theta3  |   <---drill 
  //    -        -  
  //
  // membrane strains ordered :
  //
  //            strain(0) =   eps00     i.e.   (11)-strain
  //            strain(1) =   eps11     i.e.   (22)-strain
  //            strain(2) =   gamma01   i.e.   (12)-shear
  //
  // curvatures and shear strains ordered  :
  //
  //            strain(3) =     kappa00  i.e.   (11)-curvature
  //            strain(4) =     kappa11  i.e.   (22)-curvature
  //            strain(5) =   2*kappa01  i.e. 2*(12)-curvature 
  //
  //            strain(6) =     gamma02  i.e.   (13)-shear
  //            strain(7) =     gamma12  i.e.   (23)-shear
  //
  //  same ordering for moments/shears but no 2 
  //  
  //  Then, 
  //              epsilon00 = -z * kappa00      +    eps00_membrane
  //              epsilon11 = -z * kappa11      +    eps11_membrane
  //  gamma01 = 2*epsilon01 = -z * (2*kappa01)  +  gamma01_membrane 
  //
  //  Shear strains gamma02, gamma12 constant through cross section
  //

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  static const int nstress = 8 ; //three membrane, three moment, two shear

  static const int ngauss = 4 ;

  static const int numnodes = 4 ;

  int i,  j,  k, p, q ;
  int jj, kk ;

  int success ;
  
  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static Vector strain(nstress) ;  //strain

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  //  static double Shape[3][numnodes][ngauss] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress resultants

  static Matrix dd(nstress,nstress) ;  //material tangent

  static Matrix J0(2,2) ;  //Jacobian at center
 
  static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  double epsDrill = 0.0 ;  //drilling "strain"

  double tauDrill = 0.0 ; //drilling "stress"

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;


    static Matrix Bbend(3,3) ;  // bending B matrix

    static Matrix Bshear(2,3) ; // shear B matrix

    static Matrix Bmembrane(3,2) ; // membrane B matrix


    static double BdrillJ[ndf] ; //drill B matrix

    static double BdrillK[ndf] ;  

    double *drillPointer ;

    static double saveB[nstress][ndf][numnodes] ;

  //------------------------------------------------------- 

  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;
  
//start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
  updateBasis( );
//end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)

  double dx34 = xl[0][2]-xl[0][3];
  double dy34 = xl[1][2]-xl[1][3];

  double dx21 = xl[0][1]-xl[0][0];
  double dy21 = xl[1][1]-xl[1][0];

  double dx32 = xl[0][2]-xl[0][1];
  double dy32 = xl[1][2]-xl[1][1];

  double dx41 = xl[0][3]-xl[0][0];
  double dy41 = xl[1][3]-xl[1][0];

  Matrix G(4,12);
  G.Zero();
  double one_over_four = 0.25;
  G(0,0)=-0.5;
  G(0,1)=-dy41*one_over_four;
  G(0,2)=dx41*one_over_four;
  G(0,9)=0.5;
  G(0,10)=-dy41*one_over_four;
  G(0,11)=dx41*one_over_four;
  G(1,0)=-0.5;
  G(1,1)=-dy21*one_over_four;
  G(1,2)=dx21*one_over_four;
  G(1,3)=0.5;
  G(1,4)=-dy21*one_over_four;
  G(1,5)=dx21*one_over_four;
  G(2,3)=-0.5;
  G(2,4)=-dy32*one_over_four;
  G(2,5)=dx32*one_over_four;
  G(2,6)=0.5;
  G(2,7)=-dy32*one_over_four;
  G(2,8)=dx32*one_over_four;
  G(3,6)=0.5;
  G(3,7)=-dy34*one_over_four;
  G(3,8)=dx34*one_over_four;
  G(3,9)=-0.5;
  G(3,10)=-dy34*one_over_four;
  G(3,11)=dx34*one_over_four;

  Matrix Ms(2,4);
  Ms.Zero();
  Matrix Bsv(2,12);
  Bsv.Zero();

  double Ax = -xl[0][0]+xl[0][1]+xl[0][2]-xl[0][3];
  double Bx =  xl[0][0]-xl[0][1]+xl[0][2]-xl[0][3];
  double Cx = -xl[0][0]-xl[0][1]+xl[0][2]+xl[0][3];

  double Ay = -xl[1][0]+xl[1][1]+xl[1][2]-xl[1][3];
  double By =  xl[1][0]-xl[1][1]+xl[1][2]-xl[1][3];
  double Cy = -xl[1][0]-xl[1][1]+xl[1][2]+xl[1][3];

  double alph = atan(Ay/Ax);
  double beta = 3.141592653589793/2-atan(Cx/Cy);
  Matrix Rot(2,2);
  Rot.Zero();
  Rot(0,0)=sin(beta);
  Rot(0,1)=-sin(alph);
  Rot(1,0)=-cos(beta);
  Rot(1,1)=cos(alph);
  Matrix Bs(2,12);
  
  double r1 = 0;
  double r2 = 0;
  double r3 = 0;

  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {

    r1 = Cx + sg[i]*Bx;
    r3 = Cy + sg[i]*By;
    r1 = r1*r1 + r3*r3;
	r1 = sqrt (r1);
	r2 = Ax + tg[i]*Bx;
	r3 = Ay + tg[i]*By;
	r2 = r2*r2 + r3*r3;
	r2 = sqrt (r2);

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;
    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  
    volume += dvol[i] ;

    Ms(1,0)=1-sg[i];
	Ms(0,1)=1-tg[i];
	Ms(1,2)=1+sg[i];
	Ms(0,3)=1+tg[i];
	Bsv = Ms*G;

    for ( j = 0; j < 12; j++ ) {
		Bsv(0,j)=Bsv(0,j)*r1/(8*xsj);
		Bsv(1,j)=Bsv(1,j)*r2/(8*xsj);
    }
    Bs=Rot*Bsv;

    //zero the strains
    strain.Zero( ) ;
    epsDrill = 0.0 ;


    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;
 #ifdef _SDEBUG
	if(this->getTag()==1&&i==3)
			opserr<<"Shell "<<this->getTag()<< " Bbend  "<<Bbend<<endln;
#endif
      for ( p = 0; p < 3; p++) {
		  Bshear(0,p) = Bs(0,j*3+p);
		  Bshear(1,p) = Bs(1,j*3+p);
      }//end for p

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //save the B-matrix
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  saveB[p][q][j] = BJ(p,q) ;
      }//end for p


      //nodal "displacements" 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

#ifdef _SDEBUG
	if((this->getTag())==1&&j==2)
			opserr<<"Node "<<j<<" TrialDisp "<<ul<<endln;
#endif
	  if(ul.Norm()>1e6)
				return(-1);

      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0, BJ,ul,1.0 ) ;

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	      BdrillJ[p] = *drillPointer ; //set p-th component
	      drillPointer++ ;             //pointer arithmetic
      }//end for p

      //drilling "strain" 
      for ( p = 0; p < ndf; p++ )
	      epsDrill +=  BdrillJ[p]*ul(p) ;
    } // end for j
  

    if(counterTemperature !=1&&counterTemperature !=2)
       success = materialPointers[i]->setTrialSectionDeformation( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStressResultant( ) ;
  //Liming
	if(counterTemperature==1||counterTemperature ==2){
		double Memstress11 = stress(0);
		double Memstress22 = stress(1);
		stress(0) = Memstress11- residThermal[2*i];
	   stress(1) = Memstress22- residThermal[2*i];
	   double Bendstress1 = stress(3);
		double Bendstress2 = stress(4);
	   stress(3) = Bendstress1- residThermal[2*i+1];
	   stress(4) = Bendstress2- residThermal[2*i+1];
	}
	
#ifdef _sDEBUG
		if(this->getTag()==1&&i==2)
			opserr<<"Shellthermal "<<this->getTag()<< " stress "<<stress<<endln;
#endif	
    //drilling "stress" 
    tauDrill = Ktt * epsDrill ;

    //multiply by volume element
    stress   *= dvol[i] ;
    tauDrill *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getSectionTangent( ) ;

      dd *= dvol[i] ;
    } //end if tang_flag


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numnodes; j++ ) {

      //extract BJ
      for (p=0; p<nstress; p++) {
	    for (q=0; q<ndf; q++ )
	      BJ(p,q) = saveB[p][q][j]   ;
      }//end for p

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      for ( p = 3; p < 6; p++ ) {
	    for ( q = 3; q < 6; q++ ) 
	      BJ(p,q) *= (-1.0) ;
      } //end for p

      //transpose 
      for (p=0; p<ndf; p++) {
	    for (q=0; q<nstress; q++) 
	      BJtran(p,q) = BJ(q,p) ;
      }//end for p

      residJ.addMatrixVector(0.0, BJtran,stress,1.0 ) ;

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	    BdrillJ[p] = *drillPointer ;
	    drillPointer++ ;
      }//end for p

      //residual including drill
      for ( p = 0; p < ndf; p++ )
        resid( jj + p ) += ( residJ(p) + BdrillJ[p]*tauDrill ) ;

      if ( tang_flag == 1 ) {

	    BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0 ) ;

	    for (p=0; p<ndf; p++) 
	      BdrillJ[p] *= ( Ktt*dvol[i] ) ;

        kk = 0 ;
        for ( k = 0; k < numnodes; k++ ) {

	      //extract BK
	      for (p=0; p<nstress; p++) {
	        for (q=0; q<ndf; q++ )
	          BK(p,q) = saveB[p][q][k]   ;
	      }//end for p
	  
	      //drilling B matrix
	      drillPointer = computeBdrill( k, shp ) ;
	      for (p=0; p<ndf; p++ ) {
	        BdrillK[p] = *drillPointer ;
	        drillPointer++ ;
	      }//end for p
 
          //stiffJK = BJtranD * BK  ;
	      // +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
	      stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0 ) ;

          for ( p = 0; p < ndf; p++ )  {
	        for ( q = 0; q < ndf; q++ ) {
	          stiff( jj+p, kk+q ) += stiffJK(p,q) 
		                   + ( BdrillJ[p]*BdrillK[q] ) ;
	        }//end for q
          }//end for p

          kk += ndf ;
        } // end for k loop

      } // end if tang_flag 

    jj += ndf ;
    } // end for j loop

  } //end for i gauss loop 
  
  counterTemperature++;
  return 0 ;
}


//************************************************************************
//compute local coordinates and basis

void   
ShellMITC4Thermal::computeBasis( ) 
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
  
  const Vector &coor3 = nodePointers[3]->getCrds( ) ;

  v1.Zero( ) ;
  //v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor2 ;
  v1 += coor1 ;
  v1 -= coor3 ;
  v1 -= coor0 ;
  v1 *= 0.50 ;
  
  v2.Zero( ) ;
  //v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor3 ;
  v2 += coor2 ;
  v2 -= coor1 ;
  v2 -= coor0 ;
  v2 *= 0.50 ;
 
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
  for ( i = 0; i < 4; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;
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
ShellMITC4Thermal::updateBasis( ) 
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

  const Vector &coor0 = nodePointers[0]->getCrds( ) + nodePointers[0]->getTrialDisp();

  const Vector &coor1 = nodePointers[1]->getCrds( ) + nodePointers[1]->getTrialDisp();

  const Vector &coor2 = nodePointers[2]->getCrds( ) + nodePointers[2]->getTrialDisp();
  
  const Vector &coor3 = nodePointers[3]->getCrds( ) + nodePointers[3]->getTrialDisp();

  v1.Zero( ) ;
  //v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor2 ;
  v1 += coor1 ;
  v1 -= coor3 ;
  v1 -= coor0 ;
  v1 *= 0.50 ;
  
  v2.Zero( ) ;
  //v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor3 ;
  v2 += coor2 ;
  v2 -= coor1 ;
  v2 -= coor0 ;
  v2 *= 0.50 ;
 
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
  for ( i = 0; i < 4; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;
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


//*************************************************************************
//compute Bdrill

double*
ShellMITC4Thermal::computeBdrill( int node, const double shp[3][4] )
{

  //static Matrix Bdrill(1,6) ;
  static double Bdrill[6] ;
  static double B1 ;
  static double B2 ;
  static double B6 ;


//---Bdrill Matrix in standard {1,2,3} mechanics notation---------
//
//             -                                       -
//   Bdrill = | -0.5*N,2   +0.5*N,1    0    0    0   -N |   (1x6) 
//             -                                       -  
//
//----------------------------------------------------------------

  //  Bdrill.Zero( ) ;

  //Bdrill(0,0) = -0.5*shp[1][node] ;

  //Bdrill(0,1) = +0.5*shp[0][node] ;

  //Bdrill(0,5) =     -shp[2][node] ;


  B1 =  -0.5*shp[1][node] ; 

  B2 =  +0.5*shp[0][node] ;

  B6 =  -shp[2][node] ;

  Bdrill[0] = B1*g1[0] + B2*g2[0] ;
  Bdrill[1] = B1*g1[1] + B2*g2[1] ; 
  Bdrill[2] = B1*g1[2] + B2*g2[2] ;

  Bdrill[3] = B6*g3[0] ;
  Bdrill[4] = B6*g3[1] ; 
  Bdrill[5] = B6*g3[2] ;
 
  return Bdrill ;

}


//********************************************************************
//assemble a B matrix

const Matrix&  
ShellMITC4Thermal::assembleB( const Matrix &Bmembrane,
                               const Matrix &Bbend, 
                               const Matrix &Bshear ) 
{

  //Matrix Bbend(3,3) ;  // plate bending B matrix

  //Matrix Bshear(2,3) ; // plate shear B matrix

  //Matrix Bmembrane(3,2) ; // plate membrane B matrix


    static Matrix B(8,6) ;

    static Matrix BmembraneShell(3,3) ; 
    
    static Matrix BbendShell(3,3) ; 

    static Matrix BshearShell(2,6) ;
 
    static Matrix Gmem(2,3) ;

    static Matrix Gshear(3,6) ;

    int p, q ;
    int pp ;

//    
// For Shell : 
//
//---B Matrices in standard {1,2,3} mechanics notation---------
//
//            -                     _          
//           | Bmembrane  |     0    |
//           | --------------------- |     
//    B =    |     0      |  Bbend   |   (8x6) 
//           | --------------------- |
//           |         Bshear        |
//            -           -         -
//
//-------------------------------------------------------------
//
//


    //shell modified membrane terms
    
    Gmem(0,0) = g1[0] ;
    Gmem(0,1) = g1[1] ;
    Gmem(0,2) = g1[2] ;

    Gmem(1,0) = g2[0] ;
    Gmem(1,1) = g2[1] ;
    Gmem(1,2) = g2[2] ;

    //BmembraneShell = Bmembrane * Gmem ;
    BmembraneShell.addMatrixProduct(0.0, Bmembrane,Gmem,1.0 ) ;


    //shell modified bending terms 

    Matrix &Gbend = Gmem ;

    //BbendShell = Bbend * Gbend ;
    BbendShell.addMatrixProduct(0.0, Bbend,Gbend,1.0 ) ; 


    //shell modified shear terms 

    Gshear.Zero( ) ;

    Gshear(0,0) = g3[0] ;
    Gshear(0,1) = g3[1] ;
    Gshear(0,2) = g3[2] ;

    Gshear(1,3) = g1[0] ;
    Gshear(1,4) = g1[1] ;
    Gshear(1,5) = g1[2] ;

    Gshear(2,3) = g2[0] ;
    Gshear(2,4) = g2[1] ;
    Gshear(2,5) = g2[2] ;

    //BshearShell = Bshear * Gshear ;
    BshearShell.addMatrixProduct(0.0, Bshear,Gshear,1.0 ) ;
   

  B.Zero( ) ;

  //assemble B from sub-matrices 

  //membrane terms 
  for ( p = 0; p < 3; p++ ) {

      for ( q = 0; q < 3; q++ ) 
          B(p,q) = BmembraneShell(p,q) ;

  } //end for p


  //bending terms
  for ( p = 3; p < 6; p++ ) {
    pp = p - 3 ;
    for ( q = 3; q < 6; q++ ) 
        B(p,q) = BbendShell(pp,q-3) ; 
  } // end for p
    

  //shear terms 
  for ( p = 0; p < 2; p++ ) {
      pp = p + 6 ;
      
      for ( q = 0; q < 6; q++ ) {
          B(pp,q) = BshearShell(p,q) ; 
      } // end for q
 
  } //end for p
  
  return B ;

}

//***********************************************************************
//compute Bmembrane matrix

const Matrix&   
ShellMITC4Thermal::computeBmembrane( int node, const double shp[3][4] ) 
{

  static Matrix Bmembrane(3,2) ;

//---Bmembrane Matrix in standard {1,2,3} mechanics notation---------
//
//                -             -
//               | +N,1      0   | 
// Bmembrane =   |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
//  three(3) strains and two(2) displacements (for plate)
//-------------------------------------------------------------------

  Bmembrane.Zero( ) ;

  Bmembrane(0,0) = shp[0][node] ;
  Bmembrane(1,1) = shp[1][node] ;
  Bmembrane(2,0) = shp[1][node] ;
  Bmembrane(2,1) = shp[0][node] ;

  return Bmembrane ;

}

//***********************************************************************
//compute Bbend matrix

const Matrix&   
ShellMITC4Thermal::computeBbend( int node, const double shp[3][4] )
{

    static Matrix Bbend(3,2) ;

//---Bbend Matrix in standard {1,2,3} mechanics notation---------
//
//            -             -
//   Bbend = |    0    -N,1  | 
//           |  +N,2     0   |    (3x2)
//           |  +N,1   -N,2  |
//            -             -  
//
//  three(3) curvatures and two(2) rotations (for plate)
//----------------------------------------------------------------

    Bbend.Zero( ) ;

    Bbend(0,1) = -shp[0][node] ;
    Bbend(1,0) =  shp[1][node] ;
    Bbend(2,0) =  shp[0][node] ;
    Bbend(2,1) = -shp[1][node] ; 

    return Bbend ;
}



//************************************************************************
//shape function routine for four node quads

void  
ShellMITC4Thermal::shape2d( double ss, double tt, 
		           const double x[2][4], 
		           double shp[3][4], 
		           double &xsj            )

{ 

  int i, j, k ;

  double temp ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double xs[2][2] ;
  static double sx[2][2] ;

  for ( i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      xs[i][j] = 0.0 ;

      for ( k = 0; k < 4; k++ )
	  xs[i][j] +=  x[i][k] * shp[j][k] ;

    } //end for j
  }  // end for i 

  xsj = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;

  //inverse jacobian
  double jinv = 1.0 / xsj ;
  sx[0][0] =  xs[1][1] * jinv ;
  sx[1][1] =  xs[0][0] * jinv ;
  sx[0][1] = -xs[0][1] * jinv ; 
  sx[1][0] = -xs[1][0] * jinv ; 


  //form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx[0][0] + shp[1][i]*sx[1][0] ;
    shp[1][i] = shp[0][i]*sx[0][1] + shp[1][i]*sx[1][1] ;
    shp[0][i] = temp ;
  } // end for i

  return ;
}
	  

double
ShellMITC4Thermal::shapefn2d( double ss, double tt ,int i)

{ 
   double shpVal = 0;

   switch (i){
	   case 1:
		   shpVal = 0.25*(1-ss)*(1-tt);break;
	   case 2:
		   shpVal = 0.25*(1+ss)*(1-tt);break;
	   case 3:
		   shpVal = 0.25*(1+ss)*(1+tt);break;
	   case 4:
		   shpVal = 0.25*(1-ss)*(1+tt);break;
	   default:
		   opserr<<"ShellMITC4Thermal::shapefn2d received an invalid i: "<<i <<endln ;
	   }
   return shpVal;
}

	  
//**********************************************************************

Matrix  
ShellMITC4Thermal::transpose( int dim1, 
                                       int dim2, 
		                       const Matrix &M ) 
{
  int i ;
  int j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {
     for ( j = 0; j < dim2; j++ ) 
         Mtran(j,i) = M(i,j) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  ShellMITC4Thermal::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(13);
  
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
  idData(12) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4Thermal::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(5);
  vectData(0) = Ktt;
  vectData(1) = alphaM;
  vectData(2) = betaK;
  vectData(3) = betaK0;
  vectData(4) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4Thermal::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellMITC4Thermal::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}
    
int  ShellMITC4Thermal::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(13);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4Thermal::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(8));
  connectedExternalNodes(0) = idData(9);
  connectedExternalNodes(1) = idData(10);
  connectedExternalNodes(2) = idData(11);
  connectedExternalNodes(3) = idData(12);

  static Vector vectData(5);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4Thermal::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  Ktt = vectData(0);
  alphaM = vectData(1);
  betaK = vectData(2);
  betaK0 = vectData(3);
  betaKc = vectData(4);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "ShellMITC4Thermal::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ShellMITC4Thermal::recvSelf() - material " << i << "failed to recv itself\n";
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
	  opserr << "ShellMITC4Thermal::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ShellMITC4Thermal::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}
//**************************************************************************

int
ShellMITC4Thermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(4, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
    }

    // Display mode is positive:
    // display mode = 0 -> plot no contour
    // display mode = 1-8 -> plot 1-8 stress resultant
    static Vector values(4);
    if (displayMode < 8 && displayMode > 0) {
        for (int i = 0; i < 4; i++) {
            const Vector& stress = materialPointers[i]->getStressResultant();
            values(i) = stress(displayMode - 1);
        }
    }
    else {
        for (int i = 0; i < 4; i++)
            values(i) = 0.0;
    }

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}
