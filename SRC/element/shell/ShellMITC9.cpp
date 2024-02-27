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
                                                                        
// Written: Leopoldo Tesser, Diego Talledo
// 9-node lagrandian shell element with membrane and drill

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
#include <ShellMITC9.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#define min(a,b) ( (a)<(b) ? (a):(b) )

static int numShellMITC9 = 0;

void *
OPS_ShellMITC9(void)
{
  if (numShellMITC9 == 0) {
    opserr << "Using ShellMITC9 - Developed by: Leopoldo Tesser and Diego A. Talledo\n";
    numShellMITC9++;
  }

  Element *theElement = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 11) {
    opserr << "Want: element ShellMITC9 $tag $node1 $node2 .... $node9 $secTag";
    return 0;	
  }
  
  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC9\n";
    return 0;
  }

  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(iData[10]);

  if (theSection == 0) {
    opserr << "ERROR:  element ShellMITC9 " << iData[0] << "section " << iData[10] << " not found\n";
    return 0;
  }
  
  theElement = new ShellMITC9(iData[0], iData[1], iData[2], iData[3],
			   iData[4], iData[5], iData[6], iData[7],
			   iData[8], iData[9], *theSection);

  return theElement;
}


//static data
Matrix  ShellMITC9::stiff(54,54) ;
Vector  ShellMITC9::resid(54) ; 
Matrix  ShellMITC9::mass(54,54) ;

//quadrature data
const double  ShellMITC9::root3 = sqrt(3.0) ;
const double  ShellMITC9::root3_over_root5 = root3 / sqrt(5.0) ;

double ShellMITC9::sg[9] ;
double ShellMITC9::tg[9] ;
double ShellMITC9::wg[9] ;

//null constructor
ShellMITC9::ShellMITC9( ) :
Element( 0, ELE_TAG_ShellMITC9 ),
connectedExternalNodes(9), load(0), Ki(0)
{ 
  for (int i = 0 ;  i < 9; i++ )
    materialPointers[i] = 0;

  sg[0] = -root3_over_root5 ;
  sg[1] = 0 ;
  sg[2] = root3_over_root5 ;
  sg[3] = root3_over_root5 ;
  sg[4] = root3_over_root5 ;
  sg[5] = 0 ;
  sg[6] = -root3_over_root5 ;
  sg[7] = -root3_over_root5 ;
  sg[8] = 0 ;

  tg[0] = -root3_over_root5 ;
  tg[1] = -root3_over_root5 ;
  tg[2] = -root3_over_root5 ;
  tg[3] = 0 ; 
  tg[4] = root3_over_root5 ;
  tg[5] = root3_over_root5 ;
  tg[6] = root3_over_root5 ; 
  tg[7] = 0 ;
  tg[8] = 0 ;

  wg[0] = 25.0 / 81.0 ;
  wg[1] = 40.0 / 81.0 ;
  wg[2] = 25.0 / 81.0 ;
  wg[3] = 40.0 / 81.0 ;
  wg[4] = 25.0 / 81.0 ;
  wg[5] = 40.0 / 81.0 ;
  wg[6] = 25.0 / 81.0 ;
  wg[7] = 40.0 / 81.0 ;
  wg[8] = 64.0 / 81.0 ;
}


//*********************************************************************
//full constructor
ShellMITC9::ShellMITC9(  int tag,
                   int node1,
                   int node2,
		   int node3,
                   int node4,
                   int node5,
                   int node6,
		   int node7,
                   int node8,
		   int node9,
		   SectionForceDeformation &theMaterial ) :
Element( tag, ELE_TAG_ShellMITC9 ),
connectedExternalNodes(9), load(0), Ki(0)
{
  int i;

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;
  connectedExternalNodes(4) = node5 ;
  connectedExternalNodes(5) = node6 ;
  connectedExternalNodes(6) = node7 ;
  connectedExternalNodes(7) = node8 ;
  connectedExternalNodes(8) = node9 ;

  for ( i = 0 ;  i < 9; i++ ) {
    materialPointers[i] = theMaterial.getCopy( ) ;
    if (materialPointers[i] == 0) {
      opserr << "ShellMITC9::constructor - failed to get a material of type: ShellSection\n";
    } //end if
  } //end for i


  sg[0] = -root3_over_root5 ;
  sg[1] = 0 ;
  sg[2] = root3_over_root5 ;
  sg[3] = root3_over_root5 ;
  sg[4] = root3_over_root5 ;
  sg[5] = 0 ;
  sg[6] = -root3_over_root5 ;
  sg[7] = -root3_over_root5 ;
  sg[8] = 0 ;

  tg[0] = -root3_over_root5 ;
  tg[1] = -root3_over_root5 ;
  tg[2] = -root3_over_root5 ;
  tg[3] = 0 ; 
  tg[4] = root3_over_root5 ;
  tg[5] = root3_over_root5 ;
  tg[6] = root3_over_root5 ; 
  tg[7] = 0 ;
  tg[8] = 0 ;

  wg[0] = 25.0 / 81.0 ;
  wg[1] = 40.0 / 81.0 ;
  wg[2] = 25.0 / 81.0 ;
  wg[3] = 40.0 / 81.0 ;
  wg[4] = 25.0 / 81.0 ;
  wg[5] = 40.0 / 81.0 ;
  wg[6] = 25.0 / 81.0 ;
  wg[7] = 40.0 / 81.0 ;
  wg[8] = 64.0 / 81.0 ;
}
//******************************************************************

//destructor 
ShellMITC9::~ShellMITC9( )
{
  int i ;
  
  for ( i = 0 ;  i < 9; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ;
  } //end for i

  for ( i = 0 ;  i < 9; i++ ) {
    nodePointers[i] = 0 ;
  } //end for i

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
}
//**************************************************************************


//set domain
void  ShellMITC9::setDomain( Domain *theDomain ) 
{
  int i,j ;
  static Vector eig(3) ;
  static Matrix ddMembrane(3,3) ;

  //node pointers
  for ( i = 0; i < 9; i++ ) {
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
     
     if (nodePointers[i] == 0) {
       opserr << "ShellMITC9::setDomain - no node " << connectedExternalNodes(i);
       opserr << " exists in the model\n";
     }
  }

  //compute drilling stiffness penalty parameter
  const Matrix &dd = materialPointers[0]->getInitialTangent( ) ;

  //assemble ddMembrane ;
  for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ ){
         ddMembrane(i,j) = dd(i,j) ;
	  } //end for j
  } //end for i 

  //eigenvalues of ddMembrane
  eig = LovelyEig( ddMembrane ) ;
  
  //set ktt 
  //Ktt = dd(2,2) ;  //shear modulus 
  Ktt = min( eig(2), min( eig(0), eig(1) ) ) ;
  //basis vectors and local coordinates
  computeBasis( ) ;

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  ShellMITC9::getNumExternalNodes( ) const
{
  return 9 ;
} 
 

//return connected external nodes
const ID&  ShellMITC9::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


Node **
ShellMITC9::getNodePtrs(void) 
{
  return nodePointers;
} 

//return number of dofs
int  ShellMITC9::getNumDOF( ) 
{
  return 54 ;
}


//commit state
int  ShellMITC9::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ShellMITC9::commitState () - failed in base class";
  }    

  for (int i = 0; i < 9; i++ )
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  ShellMITC9::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 9; i++ )
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  ShellMITC9::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 9; i++ )
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  ShellMITC9::Print( OPS_Stream &s, int flag )
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ShellMITC9\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
        s << "\t" << connectedExternalNodes(2) << "\t" << connectedExternalNodes(3);
        s << "\t" << connectedExternalNodes(4) << "\t" << connectedExternalNodes(5);
        s << "\t" << connectedExternalNodes(6) << "\t" << connectedExternalNodes(7);
        s << "\t" << connectedExternalNodes(8) << "\t0.00";
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
        for (i = 0; i < 9; i++) {
            const Vector &stress = materialPointers[i]->getStressResultant();
            
            s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
            for (j = 0; j < 6; j++)
                s << "\t" << stress(j);
            s << endln;
        }
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "NL Nine Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << connectedExternalNodes(0) << endln;
        s << "Node 2 : " << connectedExternalNodes(1) << endln;
        s << "Node 3 : " << connectedExternalNodes(2) << endln;
        s << "Node 4 : " << connectedExternalNodes(3) << endln;
        s << "Node 5 : " << connectedExternalNodes(4) << endln;
        s << "Node 6 : " << connectedExternalNodes(5) << endln;
        s << "Node 7 : " << connectedExternalNodes(6) << endln;
        s << "Node 8 : " << connectedExternalNodes(7) << endln;
        s << "Node 9 : " << connectedExternalNodes(8) << endln;
        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);
        
        s << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ShellMITC9\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << ", " << connectedExternalNodes(3) << ", ";
        s << connectedExternalNodes(4) << ", " << connectedExternalNodes(5) << ", ";
        s << connectedExternalNodes(6) << ", " << connectedExternalNodes(7) << ", ";
        s << connectedExternalNodes(8) << "], ";
        s << "\"section\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}

Response*
ShellMITC9::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "ShellMITC9");
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
      opserr << "ShellMITC9::setResponse() - need to specify more data\n";
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 9) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",sg[pointNum-1]);
      output.attr("neta",tg[pointNum-1]);
      
      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();
    }

  } else if (strcmp(argv[0],"stresses") ==0) {

	  for (int i=0; i<9; i++) {
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
    
    theResponse =  new ElementResponse(this, 2, Vector(72));

  } else if (strcmp(argv[0],"strains") ==0) {

	  for (int i=0; i<9; i++) {
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
    
    theResponse =  new ElementResponse(this, 3, Vector(72));
  }

  output.endTag();
  return theResponse;
}

int
ShellMITC9::getResponse(int responseID, Information &eleInfo)
{
  int i;
  int cnt = 0;
  static Vector stresses(84);
  static Vector strains(84);

  switch (responseID) {
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
    break;

  case 2: // stresses
    for (i = 0; i < 9; i++) {

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

	case 3: // strains
    for (i = 0; i < 9; i++) {

      // Get material stress response
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
}


//return stiffness matrix 
const Matrix&  ShellMITC9::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

//return secant matrix 
const Matrix&  ShellMITC9::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  static const int nstress = 8 ; //three membrane, three moment, two shear

  static const int ngauss = 9 ; 

  static const int numnodes = 9 ;

  int i,  j,  k, p, q ;
  int jj, kk ;
  int node ;

  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Matrix dd(nstress,nstress) ;  //material tangent

  //static Matrix J0(2,2) ;  //Jacobian at center
 
  //static Matrix J0inv(2,2) ; //inverse of Jacobian at center

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


  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  //computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {

    //get shape functions
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  

    volume += dvol[i] ;

    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

	  Bshear = computeBshear( j, shp ) ;

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
  return stiff ;
}

//return mass matrix
const Matrix&  ShellMITC9::getMass( ) 
{

  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 

void  ShellMITC9::zeroLoad( )
{

  if (load != 0)
    load->Zero();

  return ;
}

int 
ShellMITC9::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ShellMITC9::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int 
ShellMITC9::addInertiaLoadToUnbalance(const Vector &accel)
{
  static Vector r(54);
  int tangFlag = 1 ;

  int i;

  int allRhoZero = 0;
  for (i=0; i<9; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      allRhoZero = 1;
  }

  if (allRhoZero == 0) 
    return 0;

  formInertiaTerms( tangFlag ) ;

  int count = 0;
  for (i=0; i<9; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<6; j++)
      r(count++) = Raccel(j);
  }

  if (load == 0) 
    load = new Vector(54);
  load->addMatrixVector(1.0, mass, r, -1.0);

  return 0;
}

//get residual
const Vector&  ShellMITC9::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}

//get residual with inertia terms
const Vector&  ShellMITC9::getResistingForceIncInertia( )
{
  static Vector res(54);
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
ShellMITC9::formInertiaTerms( int tangFlag ) 
{
  //translational mass only
  //rotational inertia terms are neglected


  static const int ndf = 6 ; 

  static const int numberNodes = 9 ;

  static const int numberGauss = 9 ;

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
    //jj = 0 ;
    for ( j=0, jj=0; j<numberNodes; j++, jj+=ndf ) {

      temp = shp[massIndex][j] * dvol ;

      for ( p = 0; p < 3; p++ )
        resid( jj+p ) += ( temp * momentum(p) ) ;
      
      if ( tangFlag == 1 && rhoH != 0.0) {

	    //multiply by density
	    temp *= rhoH ;

	    //node-node translational mass
        //kk = 0 ;
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
ShellMITC9::formResidAndTangent( int tang_flag ) 
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

  static const int ngauss = 9 ;

  static const int numnodes = 9 ;

  int i,  j,  k, p, q ;
  int jj, kk ;
  int node ;

  int success ;
  
  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static Vector strain(nstress) ;  //strain

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress resultants

  static Matrix dd(nstress,nstress) ;  //material tangent

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

  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  
  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {

    //get shape functions
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;
	volume += dvol[i] ;

    //zero the strains
    strain.Zero( ) ;
    epsDrill = 0.0 ;

    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

	  Bshear = computeBshear( j, shp ) ;
	  
      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

	  //save the B-matrix
	  for (p=0; p<nstress; p++) {
		for (q=0; q<ndf; q++ ) {
		  saveB[p][q][j] = BJ(p,q) ;
		}//end for q
	  }//end for p

      //nodal "displacements" 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0, BJ,ul,1.0 ) ;

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	    //BdrillJ[p] = *drillPointer++ ;
	    BdrillJ[p] = *drillPointer ; //set p-th component
	    drillPointer++ ;             //pointer arithmetic
      }//end for p

      //drilling "strain" 
      for ( p = 0; p < ndf; p++ )
	    epsDrill +=  BdrillJ[p]*ul(p) ;
	} // end for j
  

    //send the strain to the material 
    success = materialPointers[i]->setTrialSectionDeformation( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStressResultant( ) ;

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
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	    for (q=0; q<nstress; q++) 
	      BJtran(p,q) = BJ(q,p) ;
      }//end for p

      //residJ = BJtran * stress ;
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

        //BJtranD = BJtran * dd ;
	    BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0 ) ;
        
	    for (p=0; p<ndf; p++) {
	      BdrillJ[p] *= ( Ktt*dvol[i] ) ;
        }//end for p

        kk = 0 ;
        for ( k = 0; k < numnodes; k++ ) {
          //extract BK
	      for (p=0; p<nstress; p++) {
	        for (q=0; q<ndf; q++ ){
	          BK(p,q) = saveB[p][q][k];
              
			}//end for q
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

  return ;
}


//************************************************************************
//compute local coordinates and basis
void   
ShellMITC9::computeBasis( ) 
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
  for ( i = 0; i < 9; i++ ) {
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

//*************************************************************************
//compute Bdrill
double*
ShellMITC9::computeBdrill( int node, const double shp[3][9] )
{
  static double Bdrill[6] ;
  static double B1 ;
  static double B2 ;
  static double B6 ;

//---Bdrill Matrix in standard {1,2,3} mechanics notation---------
//             -                                       -
//   Bdrill = | -0.5*N,2   +0.5*N,1    0    0    0   -N |   (1x6) 
//             -                                       -  
//----------------------------------------------------------------

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

//*************************************************************************
//assemble a B matrix
const Matrix&  
ShellMITC9::assembleB( const Matrix &Bmembrane,
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

// For Shell : 
//
//---B Matrices in standard {1,2,3} mechanics notation---------
//            -                     _          
//           | Bmembrane  |     0    |
//           | --------------------- |     
//    B =    |     0      |  Bbend   |   (8x6) 
//           | --------------------- |
//           |         Bshear        |
//            -           -         -
//-------------------------------------------------------------

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
ShellMITC9::computeBmembrane( int node, const double shp[3][9] ) 
{
  static Matrix Bmembrane(3,2) ;

//---Bmembrane Matrix in standard {1,2,3} mechanics notation---------
//                -             -
//               | +N,1      0   | 
// Bmembrane =   |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
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
ShellMITC9::computeBbend( int node, const double shp[3][9] )
{
  static Matrix Bbend(3,2) ;

//---Bbend Matrix in standard {1,2,3} mechanics notation---------
//            -             -
//   Bbend = |    0    -N,1  | 
//           |  +N,2     0   |    (3x2)
//           |  +N,1   -N,2  |
//            -             -  
//  three(3) curvatures and two(2) rotations (for plate)
//----------------------------------------------------------------

  Bbend.Zero( ) ;

  Bbend(0,1) = -shp[0][node] ;
  Bbend(1,0) =  shp[1][node] ;
  Bbend(2,0) =  shp[0][node] ;
  Bbend(2,1) = -shp[1][node] ; 

  return Bbend ;
}

//***********************************************************************
//compute standard Bshear matrix
const Matrix&  
ShellMITC9::computeBshear( int node, const double shp[3][9] )
{
  static Matrix Bshear(2,3) ;

//---Bshear Matrix in standard {1,2,3} mechanics notation------
//             -                -
//   Bshear = | +N,1      0    +N |  (2x3)
//            | +N,2     -N     0 |
//             -                -  
//  two(2) shear strains, one(1) displacement and two(2) rotations for plate 
//-----------------------------------------------------------------------

  Bshear.Zero( ) ;

  Bshear(0,0) =  shp[0][node] ;
  Bshear(0,2) =  shp[2][node] ;
  Bshear(1,0) =  shp[1][node] ;
  Bshear(1,1) = -shp[2][node] ;

  return Bshear ;
}  


//************************************************************************
//shape function routine for four node quads
void  
ShellMITC9::shape2d( double ss, double tt,const double x[2][9], 
		          double shp[3][9],double &xsj)

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

  shp[2][8] = (1-ss*ss)*(1-tt*tt) ;

  shp[2][4] = (0.5+t[1]*tt)*(1-ss*ss)-shp[2][8]/2 ;
  shp[2][5] = (0.5+s[1]*ss)*(1-tt*tt)-shp[2][8]/2 ;
  shp[2][6] = (0.5+t[3]*tt)*(1-ss*ss)-shp[2][8]/2 ;
  shp[2][7] = (0.5+s[3]*ss)*(1-tt*tt)-shp[2][8]/2 ;

  shp[2][0] = shp[2][0]-(shp[2][4]+shp[2][7])/2-shp[2][8]/4 ;
  shp[2][1] = shp[2][1]-(shp[2][5]+shp[2][4])/2-shp[2][8]/4 ;
  shp[2][2] = shp[2][2]-(shp[2][6]+shp[2][5])/2-shp[2][8]/4 ;
  shp[2][3] = shp[2][3]-(shp[2][7]+shp[2][6])/2-shp[2][8]/4 ;

  shp[0][8] = -2*ss*(1-tt*tt) ;
  shp[1][8] = -2*tt*(1-ss*ss) ;

  shp[0][4] = -2*ss*(0.5+t[1]*tt)-shp[0][8]/2 ;
  shp[0][5] = (1-tt*tt)*s[1]     -shp[0][8]/2 ;
  shp[0][6] = -2*ss*(0.5+t[3]*tt)-shp[0][8]/2 ;
  shp[0][7] = (1-tt*tt)*s[3]     -shp[0][8]/2 ;

  shp[1][4] = (1-ss*ss)*t[1]     -shp[1][8]/2 ;
  shp[1][5] = -2*tt*(0.5+s[1]*ss)-shp[1][8]/2 ;
  shp[1][6] = (1-ss*ss)*t[3]     -shp[1][8]/2 ;
  shp[1][7] = -2*tt*(0.5+s[3]*ss)-shp[1][8]/2 ;

  shp[0][0] = shp[0][0]-(shp[0][4]+shp[0][7])/2-shp[0][8]/4 ;
  shp[0][1] = shp[0][1]-(shp[0][5]+shp[0][4])/2-shp[0][8]/4 ;
  shp[0][2] = shp[0][2]-(shp[0][6]+shp[0][5])/2-shp[0][8]/4 ;
  shp[0][3] = shp[0][3]-(shp[0][7]+shp[0][6])/2-shp[0][8]/4 ;

  shp[1][0] = shp[1][0]-(shp[1][4]+shp[1][7])/2-shp[1][8]/4 ;
  shp[1][1] = shp[1][1]-(shp[1][5]+shp[1][4])/2-shp[1][8]/4 ;
  shp[1][2] = shp[1][2]-(shp[1][6]+shp[1][5])/2-shp[1][8]/4 ;
  shp[1][3] = shp[1][3]-(shp[1][7]+shp[1][6])/2-shp[1][8]/4 ;
  
  // Construct jacobian and its inverse
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      xs[i][j] = 0.0 ;
      for ( k = 0; k < 9; k++ ) {
	    xs[i][j] +=  x[i][k] * shp[j][k] ;
	  } //end for k
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
  for ( i = 0; i < 9; i++ ) {
    temp      = shp[0][i]*sx[0][0] + shp[1][i]*sx[1][0] ;
    shp[1][i] = shp[0][i]*sx[0][1] + shp[1][i]*sx[1][1] ;
    shp[0][i] = temp ;
  } // end for i
  
  return ;
}
	   
//**********************************************************************
Matrix  
ShellMITC9::transpose( int dim1,int dim2,const Matrix &M ) 
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
int  ShellMITC9::sendSelf (int commitTag,Channel &theChannel)
{
  int res = 0;
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  // Now quad sends the ids of its materials
  int matDbTag;
  static ID idData(28);
  int i;

  for (i = 0; i < 9; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+9) = matDbTag;
  }
  
  idData(18) = this->getTag();
  idData(19) = connectedExternalNodes(0);
  idData(20) = connectedExternalNodes(1);
  idData(21) = connectedExternalNodes(2);
  idData(22) = connectedExternalNodes(3);
  idData(23) = connectedExternalNodes(4);
  idData(24) = connectedExternalNodes(5);
  idData(25) = connectedExternalNodes(6);
  idData(26) = connectedExternalNodes(7);
  idData(27) = connectedExternalNodes(8);
  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC9::sendSelf() - " << this->getTag() << " failed to send ID\n";
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
    opserr << "WARNING ShellMITC9::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 9; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellMITC9::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  return res;
}
    
int  ShellMITC9::recvSelf (int commitTag,Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  static ID idData(28);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC9::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(18));
  connectedExternalNodes(0) = idData(19);
  connectedExternalNodes(1) = idData(20);
  connectedExternalNodes(2) = idData(21);
  connectedExternalNodes(3) = idData(22);
  connectedExternalNodes(4) = idData(23);
  connectedExternalNodes(5) = idData(24);
  connectedExternalNodes(6) = idData(25);
  connectedExternalNodes(7) = idData(26);
  connectedExternalNodes(8) = idData(27);
  static Vector vectData(5);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC9::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  Ktt = vectData(0);
  alphaM = vectData(1);
  betaK = vectData(2);
  betaK0 = vectData(3);
  betaKc = vectData(4);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+9);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
	     opserr << "ShellMITC9::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	     return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	     opserr << "ShellMITC9::recvSelf() - material " << i << "failed to recv itself\n";
	     return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+9);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	    delete materialPointers[i];
	    materialPointers[i] = theBroker.getNewSection(matClassTag);
	    if (materialPointers[i] == 0) {
	      opserr << "ShellMITC9::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	      exit(-1);
		}
	  }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	    opserr << "ShellMITC9::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }
  return res;
}
//**************************************************************************

int
ShellMITC9::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get the end point display coords (don't include node 9)
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
    nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
    nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
    nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
    nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(8, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v5(i);
        coords(2, i) = v2(i);
        coords(3, i) = v6(i);
        coords(4, i) = v3(i);
        coords(5, i) = v7(i);
        coords(6, i) = v4(i);
        coords(7, i) = v8(i);
    }

    // set the quantity to be displayed at the nodes;
    static Vector values(8);
    if (displayMode < 8 && displayMode > 0) {
        const Vector& stress1 = materialPointers[0]->getStressResultant();
        const Vector& stress2 = materialPointers[1]->getStressResultant();
        const Vector& stress3 = materialPointers[2]->getStressResultant();
        const Vector& stress4 = materialPointers[3]->getStressResultant();
        const Vector& stress5 = materialPointers[4]->getStressResultant();
        const Vector& stress6 = materialPointers[5]->getStressResultant();
        const Vector& stress7 = materialPointers[6]->getStressResultant();
        const Vector& stress8 = materialPointers[7]->getStressResultant();
        values(0) = stress1(displayMode - 1);
        values(1) = stress5(displayMode - 1);
        values(2) = stress2(displayMode - 1);
        values(3) = stress6(displayMode - 1);
        values(4) = stress3(displayMode - 1);
        values(5) = stress7(displayMode - 1);
        values(6) = stress4(displayMode - 1);
        values(7) = stress8(displayMode - 1);
    }
    else {
        for (int i = 0; i < 8; i++) {
            values(i) = 0.0;
        }
    }

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}
