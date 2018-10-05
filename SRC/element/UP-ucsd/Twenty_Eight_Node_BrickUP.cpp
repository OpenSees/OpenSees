/* *********************************************************************

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

//

// by Jinchi Lu and Zhaohui Yang (May 2004)

//

// 20-8 Noded TwentyEightNodeBrickUP element

//



#include <stdio.h>

#include <stdlib.h>

#include <math.h>



#include <ID.h>

#include <Vector.h>

#include <Matrix.h>

#include <Element.h>

#include <Node.h>

#include <Domain.h>

#include <ErrorHandler.h>

#include <Twenty_Eight_Node_BrickUP.h>

#include <shp3d.h>

#include <shp3dv.h>

#include <Renderer.h>

#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <Information.h>
#include <Parameter.h>

#include <Channel.h>

#include <FEM_ObjectBroker.h>

#include <elementAPI.h>

void* OPS_TwentyEightNodeBrickUP()
{
    if (OPS_GetNDM() != 3 ) {
	opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with 20_8_BrickUP element\n";
	return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < 27) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element 20_8_BrickUP eleTag? Node1? ... Node20? thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
	return 0;
    }

    // brickUPId, Node[20], matID
    int tags[22];
    int num = 22;
    if (OPS_GetIntInput(&num,tags) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(tags[21]);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "material tag: " << tags[21];
	opserr << "\nBrick element: " << tags[0] << endln;
    }

    // bk, r, perm1, perm2, perm3
    double data[5];
    num = 5;
    if (OPS_GetDoubleInput(&num,data) < 0) {
	opserr<<"WARNING: invalid double input\n";
	return 0;
    }

    double opt[3] = {0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 3) {
	num = 3;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,opt) < 0) {
	    opserr<<"WARNING: invalid double input\n";
	    return 0;
	}
    }

    return new TwentyEightNodeBrickUP(tags[0],tags[1],tags[2],tags[3],tags[4],
				      tags[5],tags[6],tags[7],tags[8],tags[9],
				      tags[10],tags[11],tags[12],tags[13],tags[14],
				      tags[15],tags[16],tags[17],tags[18],tags[19],tags[20],
				      *mat,data[0],data[1],data[2],data[3],data[4],
				      opt[0],opt[1],opt[2]);
}


//static data

double  TwentyEightNodeBrickUP::xl[3][20] ;



Matrix  TwentyEightNodeBrickUP::stiff(68,68) ;

Vector  TwentyEightNodeBrickUP::resid(68) ;

Matrix  TwentyEightNodeBrickUP::mass(68,68) ;

Matrix  TwentyEightNodeBrickUP::damp(68,68) ;



const int TwentyEightNodeBrickUP::nintu=27;

const int TwentyEightNodeBrickUP::nintp=8;

const int TwentyEightNodeBrickUP::nenu=20;

const int TwentyEightNodeBrickUP::nenp=8;

double TwentyEightNodeBrickUP::shgu[4][20][27];

double TwentyEightNodeBrickUP::shgp[4][8][8];

double TwentyEightNodeBrickUP::shgq[4][20][8];

double TwentyEightNodeBrickUP::shlu[4][20][27];

double TwentyEightNodeBrickUP::shlp[4][8][8];

double TwentyEightNodeBrickUP::shlq[4][20][8];

double TwentyEightNodeBrickUP::wu[27];

double TwentyEightNodeBrickUP::wp[8];

double TwentyEightNodeBrickUP::dvolu[27];

double TwentyEightNodeBrickUP::dvolp[8];

double TwentyEightNodeBrickUP::dvolq[8];



//null constructor

TwentyEightNodeBrickUP::TwentyEightNodeBrickUP( ) :

Element( 0, ELE_TAG_Twenty_Eight_Node_BrickUP ), materialPointers(0),

connectedExternalNodes(20), applyLoad(0), load(0), Ki(0), kc(0), rho(0)

{

  for (int i=0; i<20; i++ ) {

    nodePointers[i] = 0;

  }

  b[0] = b[1] = b[2] = 0.;

  perm[0] = perm[1] = perm[2] = 0.;



  // calculate local shape functions and derivatives

  compuLocalShapeFunction();

}





//*********************************************************************

//full constructor

TwentyEightNodeBrickUP::TwentyEightNodeBrickUP(int tag,

					       int node1,

					       int node2,

					       int node3,

					       int node4,

					       int node5,

					       int node6,

					       int node7,

					       int node8,

					       int node9,

					       int node10,

					       int node11,

					       int node12,

					       int node13,

					       int node14,

					       int node15,

					       int node16,

					       int node17,

					       int node18,

					       int node19,

					       int node20,

					       NDMaterial &theMaterial, double bulk, double rhof,

					       double p1, double p2, double p3,

					       double b1, double b2, double b3) :Element( tag, ELE_TAG_Twenty_Eight_Node_BrickUP ),

connectedExternalNodes(20), applyLoad(0), load(0), Ki(0), kc(bulk), rho(rhof)

{

	connectedExternalNodes(0) = node1 ;

	connectedExternalNodes(1) = node2 ;

	connectedExternalNodes(2) = node3 ;

	connectedExternalNodes(3) = node4 ;

	connectedExternalNodes(4) = node5 ;

	connectedExternalNodes(5) = node6 ;

	connectedExternalNodes(6) = node7 ;

	connectedExternalNodes(7) = node8 ;

	connectedExternalNodes(8) = node9 ;

	connectedExternalNodes(9) = node10 ;

	connectedExternalNodes(10) = node11 ;

	connectedExternalNodes(11) = node12 ;

	connectedExternalNodes(12) = node13 ;

	connectedExternalNodes(13) = node14 ;

	connectedExternalNodes(14) = node15 ;

	connectedExternalNodes(15) = node16 ;

	connectedExternalNodes(16) = node17 ;

	connectedExternalNodes(17) = node18 ;

	connectedExternalNodes(18) = node19 ;

	connectedExternalNodes(19) = node20 ;



	int i ;

    // Allocate arrays of pointers to NDMaterials

    materialPointers = new NDMaterial *[nintu];



    if (materialPointers == 0) {

      opserr << "TwentyEightNodeBrickUP::TwentyEightNodeBrickUP - failed allocate material model pointer\n";

      exit(-1);

    }

    for ( i=0; i<nintu; i++ ) {



      materialPointers[i] = theMaterial.getCopy("ThreeDimensional") ;



      if (materialPointers[i] == 0) {

	opserr <<"TwentyEightNodeBrickUP::constructor - failed to get a material of type: ThreeDimensional\n";

	exit(-1);

      } //end if



    } //end for i



    // Body forces

    b[0] = b1;

    b[1] = b2;

    b[2] = b3;

    // Permeabilities

    perm[0] = p1;

    perm[1] = p2;

    perm[2] = p3;

    //printf("b %15.6e %15.6e %15.6e perm %15.6e %15.6e %15.6e\n", b1, b2,b3,p1,p2,p3);

    // calculate local shape functions and derivatives

    compuLocalShapeFunction();



}

//******************************************************************





//destructor

TwentyEightNodeBrickUP::~TwentyEightNodeBrickUP( )

{

  int i ;

  for ( i = 0; i < nintu; i++) {

    if (materialPointers[i])

      delete materialPointers[i];

  }



  // Delete the array of pointers to NDMaterial pointer arrays

  if (materialPointers)

    delete [] materialPointers;



  for ( i=0 ; i<nenu; i++ ) {

    nodePointers[i] = 0 ;

  } //end for i



  if (load != 0)

    delete load;



  if (Ki != 0)

    delete Ki;

}





//set domain

void  TwentyEightNodeBrickUP::setDomain( Domain *theDomain )

{

  int i,dof ;



  // Check Domain is not null - invoked when object removed from a domain

  if (theDomain == 0) {

    for ( i=0; i<nenu; i++ )

      nodePointers[i] = 0;

    return;

  }



  //node pointers

  for ( i=0; i<nenu; i++ ) {

    nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

    if (nodePointers[i] == 0) {

      opserr << "FATAL ERROR TwentyEightNodeBrickUP ("<<this->getTag()<<"): node not found in domain"<<endln;

      return;

    }



    dof = nodePointers[i]->getNumberDOF();

    if( (i<nenp && dof != 4) || (i>=nenp && dof != 3) ) {

      opserr << "FATAL ERROR TwentyEightNodeBrickUP ("<<this->getTag()<<"): has wrong number of DOFs at its nodes"<<endln;

      return;

    }

  }

  this->DomainComponent::setDomain(theDomain);

}





//get the number of external nodes

int  TwentyEightNodeBrickUP::getNumExternalNodes( ) const

{

  return nenu ;

}





//return connected external nodes

const ID&  TwentyEightNodeBrickUP::getExternalNodes( )

{

	return connectedExternalNodes ;

}



//return connected external node

Node **

TwentyEightNodeBrickUP::getNodePtrs(void)

{

	return nodePointers ;

}





//return number of dofs

int  TwentyEightNodeBrickUP::getNumDOF( )

{

	return 68 ;

}





//commit state

int  TwentyEightNodeBrickUP::commitState( )

{

	int success = 0 ;



	// call element commitState to do any base class stuff

	if ((success = this->Element::commitState()) != 0) {

		opserr << "TwentyEightNodeBrickUP::commitState () - failed in base class";

	}



	for (int i=0; i<nintu; i++ )

		success += materialPointers[i]->commitState( ) ;



	return success ;

}







//revert to last commit

int  TwentyEightNodeBrickUP::revertToLastCommit( )

{

	int i ;

	int success = 0 ;



	for ( i=0; i<nintu; i++ )

		success += materialPointers[i]->revertToLastCommit( ) ;



	return success ;

}





//revert to start

int  TwentyEightNodeBrickUP::revertToStart( )

{

	int i ;

	int success = 0 ;



	for ( i=0; i<nintu; i++ )

		success += materialPointers[i]->revertToStart( ) ;



	return success ;

}



//print out element data

void  TwentyEightNodeBrickUP::Print( OPS_Stream &s, int flag )

{



	if (flag == 2) {



		s << "#20_8_BrickUP\n";



		int i;

		const int numNodes = 20;

		const int nstress = 6 ;



		for (i=0; i<numNodes; i++) {

			const Vector &nodeCrd = nodePointers[i]->getCrds();

			const Vector &nodeDisp = nodePointers[i]->getDisp();

			s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)

				<< " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;

		}



		// spit out the section location & invoke print on the scetion

		const int numMaterials = nintu;



		static Vector avgStress(7);

		static Vector avgStrain(nstress);

		avgStress.Zero();

		avgStrain.Zero();

		for (i=0; i<numMaterials; i++) {

			avgStress += materialPointers[i]->getStress();

			avgStrain += materialPointers[i]->getStrain();

		}

		avgStress /= numMaterials;

		avgStrain /= numMaterials;



		s << "#AVERAGE_STRESS ";

		for (i=0; i<7; i++)

			s << avgStress(i) << " " ;

		s << endln;



		s << "#AVERAGE_STRAIN ";

		for (i=0; i<nstress; i++)

			s << avgStrain(i) << " " ;

		s << endln;



		/*

		for (i=0; i<numMaterials; i++) {

		s << "#MATERIAL\n";

		//      materialPointers[i]->Print(s, flag);

		s << materialPointers[i]->getStress();

		}

		*/



	} else {



		s << endln ;

		s << "20-8 Noded TwentyEightNodeBrickUP \n" ;

		s << "Element Number: " << this->getTag() << endln ;

		s << "Node 1 : " << connectedExternalNodes(0) << endln ;

		s << "Node 2 : " << connectedExternalNodes(1) << endln ;

		s << "Node 3 : " << connectedExternalNodes(2) << endln ;

		s << "Node 4 : " << connectedExternalNodes(3) << endln ;

		s << "Node 5 : " << connectedExternalNodes(4) << endln ;

		s << "Node 6 : " << connectedExternalNodes(5) << endln ;

		s << "Node 7 : " << connectedExternalNodes(6) << endln ;

		s << "Node 8 : " << connectedExternalNodes(7) << endln ;

		s << "Node 9 : " << connectedExternalNodes(8) << endln ;

		s << "Node 10 : " << connectedExternalNodes(9) << endln ;

		s << "Node 11 : " << connectedExternalNodes(10) << endln ;

		s << "Node 12 : " << connectedExternalNodes(11) << endln ;

		s << "Node 13 : " << connectedExternalNodes(12) << endln ;

		s << "Node 14 : " << connectedExternalNodes(13) << endln ;

		s << "Node 15 : " << connectedExternalNodes(14) << endln ;

		s << "Node 16 : " << connectedExternalNodes(15) << endln ;

		s << "Node 17 : " << connectedExternalNodes(16) << endln ;

		s << "Node 18 : " << connectedExternalNodes(17) << endln ;

		s << "Node 19 : " << connectedExternalNodes(18) << endln ;

		s << "Node 20 : " << connectedExternalNodes(19) << endln ;



		s << "Material Information : \n " ;

		materialPointers[0]->Print( s, flag ) ;



		s << endln ;

	}

}



int

TwentyEightNodeBrickUP::update()

{

	int i, j, k, k1;

	static double u[3][20];

	static double xsj;

	static Matrix B(6, 3);

	double volume = 0.;



	for (i = 0; i < nenu; i++) {

	     const Vector &disp = nodePointers[i]->getTrialDisp();

	     u[0][i] = disp(0);

	     u[1][i] = disp(1);

	     u[2][i] = disp(2);

    }



	static Vector eps(6);



	int ret = 0;



	//compute basis vectors and local nodal coordinates

	computeBasis( ) ;



	for( i = 0; i < nintu; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 0);

		//volume element to also be saved

		dvolu[i] = wu[i] * xsj ;

		volume += dvolu[i];

	} // end for i

    //printf("volume = %f\n", volume);



	// Loop over the integration points

	for (i = 0; i < nintu; i++) {



		// Interpolate strains

		//eps = B*u;

		//eps.addMatrixVector(0.0, B, u, 1.0);

		eps.Zero();

		for ( j = 0; j < nenu; j++) {





			B(0,0) = shgu[0][j][i];

			B(0,1) = 0.;

			B(0,2) = 0.;

			B(1,0) = 0.;

			B(1,1) = shgu[1][j][i];

			B(1,2) = 0.;

			B(2,0) = 0.;

			B(2,1) = 0.;

			B(2,2) = shgu[2][j][i];

			B(3,0) = shgu[1][j][i];

			B(3,1) = shgu[0][j][i];

			B(3,2) = 0.;

			B(4,0) = 0.;

			B(4,1) = shgu[2][j][i];

			B(4,2) = shgu[1][j][i];

			B(5,0) = shgu[2][j][i];

			B(5,1) = 0.;

			B(5,2) = shgu[0][j][i];



			//BJ = computeB( j, shp ) ;



			//nodal displacements

			const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

			Vector ul3(3);

			ul3(0) = ul(0);

			ul3(1) = ul(1);

			ul3(2) = ul(2);

			//compute the strain

			//strain += (BJ*ul) ;

			eps.addMatrixVector(1.0,B,ul3,1.0 ) ;



			/* for( k = 0; k < 6; k++) {

			for( k1 = 0; k1 < 3; k1++) {

			eps[k] += BJ(k, k1)*u[k1][j];

			}

		}	*/





		}



		// Set the material strain

		ret += materialPointers[i]->setTrialStrain(eps);

	}



	return ret;

}



//return tangent stiffness matrix

const Matrix&  TwentyEightNodeBrickUP::getTangentStiff( )

{

	return getStiff( 1 );

}



// return initial stiffness matrix

const Matrix&  TwentyEightNodeBrickUP::getInitialStiff( )

{

	return getStiff( 0 );

}



// compute stiffness matrix

const Matrix&  TwentyEightNodeBrickUP::getStiff( int flag )

{

	if (flag != 0 && flag != 1) {

		opserr << "FATAL TwentyEightNodeBrickUP::getStiff() - illegal use\n";

		exit(-1);

	}



	if (flag == 0 && Ki != 0)

		return *Ki;



	int i, j ;



	static double xsj ;  // determinant jacaobian matrix

	double volume = 0.;

	//-------------------------------------------------------

	int j3, j3m1, j3m2, ik, ib, jk, jb;

	static Matrix B(6,nenu*3);

	static Matrix BTDB(nenu*3,nenu*3);

	static Matrix D(6, 6);

	B.Zero();

	BTDB.Zero();

	stiff.Zero();



	//compute basis vectors and local nodal coordinates

	computeBasis( ) ;



	for( i = 0; i < nintu; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 0);

		//volume element to also be saved

		dvolu[i] = wu[i] * xsj ;

		volume += dvolu[i];

	} // end for i

    //printf("volume = %f\n", volume);



//	for( i = 0; i < nintu; i++ ) {

//		for(int j = 0; j < nenu; j++ ) {

//			printf("%5d %5d %15.6e %15.6e %15.6e %15.6e\n", i,j,

//				shgu[0][j][i], shgu[1][j][i], shgu[2][j][i], shgu[3][j][i]);

//		}

//	}

//	exit(-1);



	// Loop over the integration points

	for (i = 0; i < nintu; i++) {



		// Get the material tangent

		if( flag == 0 )

			D = materialPointers[i]->getInitialTangent();

		else

			D = materialPointers[i]->getTangent();

		//const Matrix &D = materialPointers[i]->getTangent();





		for (j=0; j<nenu; j++) {



			j3   = 3*j+2;

			j3m1 = j3 - 1;

			j3m2 = j3 - 2;



			B(0,j3m2) = shgu[0][j][i];

			B(0,j3m1) = 0.;

			B(0,j3  ) = 0.;



			B(1,j3m2) = 0.;

			B(1,j3m1) = shgu[1][j][i];

			B(1,j3  ) = 0.;



			B(2,j3m2) = 0.;

			B(2,j3m1) = 0.;

			B(2,j3  ) = shgu[2][j][i];



			B(3,j3m2) = shgu[1][j][i];

			B(3,j3m1) = shgu[0][j][i];

			B(3,j3  ) = 0.;



			B(4,j3m2) = 0.;

			B(4,j3m1) = shgu[2][j][i];

			B(4,j3  ) = shgu[1][j][i];



			B(5,j3m2) = shgu[2][j][i];

			B(5,j3m1) = 0.;

			B(5,j3  ) = shgu[0][j][i];



		}



		// Perform numerical integration

		//K = K + (B^ D * B) * intWt(i) * detJ;

		BTDB.addMatrixTripleProduct(1.0, B, D, dvolu[i]);

	}



	for (i = 0; i < nenu; i++) {

		if (i<nenp)

			ik = i*4;

		else

			ik = nenp*4 + (i-nenp)*3;

		ib = i*3;



		for (j = 0; j < nenu; j++) {

			if (j<nenp)

				jk = j*4;

			else

				jk = nenp*4 + (j-nenp)*3;

			jb = j*3;

			for( int i1 = 0; i1 < 3; i1++)

				for(int j1 = 0; j1 < 3; j1++) {

					stiff(ik+i1, jk+j1) = BTDB(ib+i1,jb+j1);

				}

		}

	}

	if( flag == 1) {

		return stiff;

	}

	Ki = new Matrix(stiff);

	if (Ki == 0) {

		opserr << "FATAL TwentyEightNodeBrickUP::getStiff() -";

		opserr << "ran out of memory\n";

		exit(-1);

	}



	return *Ki;

}





//return mass matrix

const Matrix&  TwentyEightNodeBrickUP::getMass( )

{

	int tangFlag = 1 ;



	formInertiaTerms( tangFlag ) ;



	return mass ;

}





//return mass matrix

const Matrix&  TwentyEightNodeBrickUP::getDamp( )

{

	int tangFlag = 1 ;



	formDampingTerms( tangFlag ) ;



	return damp ;

}



void TwentyEightNodeBrickUP::formDampingTerms( int tangFlag )

{

	static double xsj ;  // determinant jacaobian matrix

	int i, j, k, m, ik, jk;

	double volume = 0.;

	//zero damp

	damp.Zero( ) ;



	if (betaK != 0.0)

		damp.addMatrix(1.0, this->getTangentStiff(), betaK);

	if (betaK0 != 0.0)

		damp.addMatrix(1.0, this->getInitialStiff(), betaK0);

	if (betaKc != 0.0)

		damp.addMatrix(1.0, *Kc, betaKc);



	if (alphaM != 0.0) {

		this->getMass();

		for( i = 0; i < nenu; i++ ) {

			if( i < nenp)

				ik = i*4;

			else

				ik = nenp*4 + (i - nenp) * 3;

			for( j = 0; j < nenu; j++) {

				if( j < nenp)

					jk = j * 4;

				else

					jk = nenp * 4 + (j-nenp) * 3;

				for( k = 0; k < 3; k++)

					damp(ik + k, jk + k) += mass(ik + k, jk + k) * alphaM;

			}

		}

	}



	//compute basis vectors and local nodal coordinates

	computeBasis( ) ;

	//gauss loop to compute and save shape functions

	for( i = 0; i < nintp; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 1);

		//volume element to also be saved

		dvolp[i] = wp[i] * xsj ;

		volume += dvolp[i];

	} // end for i

	//printf("volume = %f\n", volume);



	volume = 0.;

	for( i = 0; i < nintp; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 2);

		//volume element to also be saved

		dvolq[i] = wp[i] * xsj ;

		volume += dvolq[i];

	} // end for i

	//printf("volume = %f\n", volume);



	// Compute coupling matrix

	for( i = 0; i < nenu; i++ ) {

		if( i < nenp)

			ik = i * 4;

		else

			ik = nenp * 4 + (i-nenp)*3;

		for( j = 0; j < nenp; j++) {

			jk = j * 4 + 3;

			for( m = 0; m < nintp; m++) {

				for( k = 0; k < 3; k++) {

					damp(ik+k,jk) += -dvolq[m]*shgq[k][i][m]*shgp[3][j][m];

				}

			}

			for( k = 0; k < 3; k++ ) {

				damp(jk, ik+k) = damp(ik+k, jk);

			}

		}

	}

	// Compute permeability matrix

	for( i = 0; i < nenp; i++ ) {

		ik = i*4 + 3;

		for( j = 0; j < nenp; j++ ) {

			jk = j * 4 + 3;

			for( m = 0; m < nintp; m++ ) {

				damp(ik,jk) += - dvolp[m]*(perm[0]*shgp[0][i][m]*shgp[0][j][m] +

					perm[1]*shgp[1][i][m]*shgp[1][j][m]+

					perm[2]*shgp[2][i][m]*shgp[2][j][m]);

			}

		}

	}

}





void  TwentyEightNodeBrickUP::zeroLoad( )

{
	if (load != 0) {
		load->Zero();
    }
  	applyLoad = 0;

  	appliedB[0] = 0.0;
 	appliedB[1] = 0.0;
 	appliedB[2] = 0.0;

	return ;
}





int

TwentyEightNodeBrickUP::addLoad(ElementalLoad *theLoad, double loadFactor)

{
	// Added option for applying body forces in load pattern: C.McGann, U.Washington
  	int type;
  	const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_BrickSelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor * b[0];
        appliedB[1] += loadFactor * b[1];
        appliedB[2] += loadFactor * b[2];
      return 0;
    } else if (type == LOAD_TAG_SelfWeight) {
        // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
        applyLoad = 1;
	    appliedB[0] += loadFactor*data(0)*b[0];
	    appliedB[1] += loadFactor*data(1)*b[1];
	    appliedB[2] += loadFactor*data(2)*b[2];
	    return 0;
  	} else {
    	opserr << "TwentyEightNodeBrickUP::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    	return -1;
  	}

  	return -1;
}



int

TwentyEightNodeBrickUP::addInertiaLoadToUnbalance(const Vector &accel)

{

	static Vector ra(68);

	int i, j, ik;

	ra.Zero();



	for( i = 0; i < nenu; i++) {

		const Vector &Raccel = nodePointers[i]->getRV(accel);

		if ((i<nenp && 4 != Raccel.Size()) || (i>=nenp && 3 != Raccel.Size())) {

			opserr << "TwentyEightNodeBrickUP::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";

			return -1;

		}



		if (i<nenp)

			ik = i*4;

		else

			ik = nenp*4 + (i-nenp)*3;



		ra[ik] = Raccel(0);

		ra[ik+1] = Raccel(1);

		ra[ik+2] = Raccel(2);

	}



	// Compute mass matrix

	int tangFlag = 1 ;

	formInertiaTerms( tangFlag ) ;



	// create the load vector if one does not exist

	if (load == 0)

	  load = new Vector(68);



	// add -M * RV(accel) to the load vector

	load->addMatrixVector(1.0, mass, ra, -1.0);

	//for( i = 0; i < 68; i++) {

	//	for( j = 0; j < 68; j++)

	//				load(i) += -mass(i,j)*ra[j];

	//}



	return 0;

}





//get residual

const Vector&  TwentyEightNodeBrickUP::getResistingForce( )

{

	int i, j, jk, k, k1;

	double xsj;

	static Matrix B(6, 3);

	double volume = 0.;



//	printf("calling getResistingForce()\n");

	resid.Zero();



	//compute basis vectors and local nodal coordinates

	computeBasis( ) ;

	//gauss loop to compute and save shape functions

	for( i = 0; i < nintu; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 0);

		//volume element to also be saved

		dvolu[i] = wu[i] * xsj ;

		volume += dvolu[i];

	} // end for i

	//printf("volume = %f\n", volume);

	volume = 0.;

	for( i = 0; i < nintp; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 1);

		//volume element to also be saved

		dvolp[i] = wp[i] * xsj ;

		volume += dvolp[i];

	} // end for i

	//printf("volume = %f\n", volume);



	// Loop over the integration points

	for (i = 0; i < nintu; i++) {



		// Get material stress response

		const Vector &sigma = materialPointers[i]->getStress();



		// Perform numerical integration on internal force

		//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;

		//P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);

		for (j = 0; j < nenu; j++) {

			if (j<nenp)

				jk = j*4;

			else

				jk = nenp*4 + (j-nenp)*3;



			B(0,0) = shgu[0][j][i];

			B(0,1) = 0.;

			B(0,2) = 0.;

			B(1,0) = 0.;

			B(1,1) = shgu[1][j][i];

			B(1,2) = 0.;

			B(2,0) = 0.;

			B(2,1) = 0.;

			B(2,2) = shgu[2][j][i];

			B(3,0) = shgu[1][j][i];

			B(3,1) = shgu[0][j][i];

			B(3,2) = 0.;

			B(4,0) = 0.;

			B(4,1) = shgu[2][j][i];

			B(4,2) = shgu[1][j][i];

			B(5,0) = shgu[2][j][i];

			B(5,1) = 0.;

			B(5,2) = shgu[0][j][i];





			for(k = 0; k < 3; k++) {

				for(k1 = 0; k1 < 6; k1++)

					resid(jk+k) += dvolu[i]*(B(k1,k)*sigma(k1));

			}

			// Subtract equiv. body forces from the nodes

			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;

			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);

			double r = mixtureRho(i);

			if (applyLoad == 0) {
				resid(jk) -= dvolu[i]*(shgu[3][j][i]*r*b[0]);
				resid(jk+1) -= dvolu[i]*(shgu[3][j][i]*r*b[1]);
				resid(jk+2) -= dvolu[i]*(shgu[3][j][i]*r*b[2]);
			} else {
				resid(jk) -= dvolu[i]*(shgu[3][j][i]*r*appliedB[0]);
				resid(jk+1) -= dvolu[i]*(shgu[3][j][i]*r*appliedB[1]);
				resid(jk+2) -= dvolu[i]*(shgu[3][j][i]*r*appliedB[2]);
			}

		}

	}



	// Subtract fluid body force

	for (j = 0; j < nenp; j++) {

		jk = j*4+3;

		for (i = 0; i < nintp; i++) {
			
			if (applyLoad == 0) {
				resid(jk) += dvolp[i]*rho*(perm[0]*b[0]*shgp[0][j][i] +
							 perm[1]*b[1]*shgp[1][j][i] + perm[2]*b[2]*shgp[2][j][i]);
			} else {
				resid(jk) += dvolp[i]*rho*(perm[0]*appliedB[0]*shgp[0][j][i] +
							 perm[1]*appliedB[1]*shgp[1][j][i] + perm[2]*appliedB[2]*shgp[2][j][i]);
			}

		}

	}



	// Subtract other external nodal loads ... P_res = P_int - P_ext

//	opserr<<"resid before:"<<resid<<endln;



	if (load != 0)

		resid -= *load;



//	opserr<<"resid "<<resid<<endln;



	return resid ;

}





//get residual with inertia terms

const Vector&  TwentyEightNodeBrickUP::getResistingForceIncInertia( )

{

	static Vector res(68);



	int i, j, ik;

	static double a[68];



	for (i=0; i<nenu; i++) {

		const Vector &accel = nodePointers[i]->getTrialAccel();

		if ((i<nenp && 4 != accel.Size()) || (i>=nenp && 3 != accel.Size())) {

			opserr << "TwentyEightNodeBrickUP::getResistingForceIncInertia matrix and vector sizes are incompatible\n";

			exit(-1);

		}



		if (i<nenp)

			ik = i*4;

		else

			ik = nenp*4 + (i-nenp)*3;

		a[ik] = accel(0);

		a[ik+1] = accel(1);

		a[ik+2] = accel(2);

		if (i<nenp) a[ik+3] = accel(3);

	}

	// Compute the current resisting force

	this->getResistingForce();

//	opserr<<"K "<<resid<<endln;



	// Compute the mass matrix

	this->getMass();



	for (i = 0; i < 68; i++) {

		for (j = 0; j < 68; j++){

			resid(i) += mass(i,j)*a[j];

		}

	}

//	printf("\n");

	//opserr<<"K+M "<<P<<endln;





	for (i=0; i<nenu; i++) {

		const Vector &vel = nodePointers[i]->getTrialVel();

		if ((i<nenp && 4 != vel.Size()) || (i>=nenp && 3 != vel.Size())) {

			opserr << "TwentyEightNodeBrickUP::getResistingForceIncInertia matrix and vector sizes are incompatible\n";

			exit(-1);

		}



		if (i<nenp)

			ik = i*4;

		else

			ik = nenp*4 + (i-nenp)*3;

		a[ik] = vel(0);

		a[ik+1] = vel(1);

		a[ik+2] = vel(2);

		if (i<nenp) a[ik+3] = vel(3);

	}



	this->getDamp();



	for (i = 0; i < 68; i++) {

		for (j = 0; j < 68; j++) {

			resid(i) += damp(i,j)*a[j];

		}

	}



	res = resid;

//	opserr<<"res "<<res<<endln;



	return res;

}





//*********************************************************************

//form inertia terms



void   TwentyEightNodeBrickUP::formInertiaTerms( int tangFlag )

{

	static double xsj ;  // determinant jacaobian matrix

	int i, j, k, ik, m, jk;

	double Nrho;



	//zero mass

	mass.Zero( ) ;



	//compute basis vectors and local nodal coordinates

	computeBasis( ) ;

	//gauss loop to compute and save shape functions

	for( i = 0; i < nintu; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 0);

		//volume element to also be saved

		dvolu[i] = wu[i] * xsj ;

	} // end for i

	for( i = 0; i < nintp; i++ ) {

		// compute Jacobian and global shape functions

		Jacobian3d(i, xsj, 1);

		//volume element to also be saved

		dvolp[i] = wp[i] * xsj ;

	} // end for i



	// Compute consistent mass matrix

	for (i = 0; i < nenu; i++) {

		if (i<nenp)

			ik = i*4;

		else

			ik = nenp*4 + (i-nenp)*3;



		for (j = 0; j < nenu; j++) {

			if (j<nenp)

				jk = j*4;

			else

				jk = nenp*4 + (j-nenp)*3;



			for (m = 0; m < nintu; m++) {

				Nrho = dvolu[m]*mixtureRho(m)*shgu[3][i][m]*shgu[3][j][m];

				for( k = 0; k < 3; k++) {

					mass(ik+k,jk+k) += Nrho;

				}

			}

		}

	}



	// Compute compressibility matrix

	double oneOverKc = 1./kc;

	for (i = 0; i < nenp; i++) {

		ik = i*4+3;



		for (j = 0; j < nenp; j++) {

			jk = j*4+3;



			for (m = 0; m < nintp; m++) {

				mass(ik,jk) += -dvolp[m]*oneOverKc*shgp[3][i][m]*shgp[3][j][m];

			}

		}

	}

}





double TwentyEightNodeBrickUP::mixtureRho(int i)

{

	double rhoi;



	rhoi= materialPointers[i]->getRho();

	//e = 0.7;  //theMaterial[i]->getVoidRatio();

	//n = e / (1.0 + e);

	//return n * rho + (1.0-n) * rhoi;

	return rhoi;

}



//************************************************************************

//compute local coordinates and basis



void   TwentyEightNodeBrickUP::computeBasis( )

{



	//nodal coordinates



	int i ;

	for ( i = 0; i < nenu; i++ ) {



		const Vector &coorI = nodePointers[i]->getCrds( ) ;



		xl[0][i] = coorI(0) ;

		xl[1][i] = coorI(1) ;

		xl[2][i] = coorI(2) ;



	}  //end for i



}





//**********************************************************************



int  TwentyEightNodeBrickUP::sendSelf (int commitTag, Channel &theChannel)

{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // TwentyEightNodeBrickUP packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(13);
  data(0) = this->getTag();
  data(1) = rho;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = b[2];

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;

  data(9) = kc;
  data(10) = perm[0];
  data(11) = perm[1];
  data(12) = perm[2];

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING TwentyEightNodeBrickUP::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }

  // Now TwentyEightNodeBrickUP sends the ids of its materials
  int matDbTag;

  static ID idData(74);

  int i;
  for (i = 0; i < nintu; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
		  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+nintu) = matDbTag;
 }

 for( i = 0; i < 20; i++)
   idData(54+i) = connectedExternalNodes(i);

 res += theChannel.sendID(dataTag, commitTag, idData);
 if (res < 0) {
   opserr << "WARNING TwentyEightNodeBrickUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
   return res;
 }

 // Finally, TwentyEightNodeBrickUP asks its material objects to send themselves
 for (i = 0; i < nintu; i++) {
   res += materialPointers[i]->sendSelf(commitTag, theChannel);
   if (res < 0) {
     opserr << "WARNING TwentyEightNodeBrickUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
     return res;
   }
 }

 return res;

}



int  TwentyEightNodeBrickUP::recvSelf (int commitTag,

									   Channel &theChannel,

									   FEM_ObjectBroker &theBroker)

{
  int res = 0;

  int dataTag = this->getDbTag();

  // TwentyEightNodeBrickUP creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(13);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING TwentyEightNodeBrickUP::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  rho = data(1);
  b[0] = data(2);
  b[1] = data(3);
  b[2] = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  kc = data(9);
  perm[0] = data(10);
  perm[1] = data(11);
  perm[2] = data(12);

  static ID idData(74);
  // TwentyEightNodeBrickUP now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING TwentyEightNodeBrickUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  int i;

  for( i = 0; i < 20; i++)
    connectedExternalNodes(i) = idData(54+i);

  if (materialPointers == 0) {
    // Allocate new materials
    materialPointers = new NDMaterial* [nintu];
    if (materialPointers == 0) {
      opserr << "TwentyEightNodeBrickUP::recvSelf() - Could not allocate NDMaterial array\n";
      return -1;
    }
    for (i = 0; i < nintu; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nintu);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	    opserr << "TwentyEightNodeBrickUP::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	    return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	    opserr << "TwentyEightNodeBrickUP::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }
  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (i = 0; i < nintu; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nintu);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	    delete materialPointers[i];
	    materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	    if (materialPointers[i] == 0) {
	      opserr << "TwentyEightNodeBrickUP::recvSelf() - Broker could not create NDMaterial of class type " <<
	      matClassTag << endln;
	      exit(-1);
	    }
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	    opserr << "TwentyEightNodeBrickUP::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  return res;

}

//**************************************************************************



int

TwentyEightNodeBrickUP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)

{

   return 0;

}


int
TwentyEightNodeBrickUP::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // permeability in horizontal direction
  if (strcmp(argv[0],"hPerm") == 0) {
    return param.addObject(3, this);

  // permeability in vertical direction
  } else if (strcmp(argv[0],"vPerm") == 0) {
    return param.addObject(4, this);

  } else {

    int matRes = res;
    for (int i=0; i<nintu; i++) {
        matRes =  materialPointers[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
  }

    return res;
}
    
int
TwentyEightNodeBrickUP::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;
  switch( parameterID ) {
	case 3:
		perm[0] = info.theDouble;
		this->getDamp();	// update mass matrix
		return 0;
	case 4:
		perm[1] = info.theDouble;
		perm[2] = info.theDouble;
		this->getDamp();	// update mass matrix
		return 0;
	default:
		return -1;
  }

  return -1;
}


Response*
TwentyEightNodeBrickUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","Twenty_Eight_Node_BrickUP");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=20; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData, nodePointers[i-1]->getTag());
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {


    for (int i=1; i<=20; i++) {
      sprintf(outputData,"P1_",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_",i);
      output.tag("ResponseType",outputData);
      if (i <= nenp) {
	sprintf(outputData,"Pp_",i);
	output.tag("ResponseType",outputData);
      }
    }

    theResponse = new ElementResponse(this, 1, resid);

  } else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)

    theResponse = new ElementResponse(this, 2, stiff);



  else if (strcmp(argv[0],"mass") == 0)

    theResponse = new ElementResponse(this, 3, mass);



  else if (strcmp(argv[0],"damp") == 0)

    theResponse = new ElementResponse(this, 4, damp);



  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

    int pointNum = atoi(argv[1]);

    if (pointNum > 0 && pointNum <= nintu) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag(); // GaussPoint
    }
  } else if (strcmp(argv[0],"stresses") ==0) {

    for (int i=0; i<nintu; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma13");
      output.tag("ResponseType","sigma23");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }

    theResponse = new ElementResponse(this, 5, Vector(nintu*6));

  }


  output.endTag(); // ElementOutput
  return theResponse;
}



int

TwentyEightNodeBrickUP::getResponse(int responseID, Information &eleInfo)

{

	static Vector stresses(nintu*6);



	if (responseID == 1)

		return eleInfo.setVector(this->getResistingForce());



	else if (responseID == 2)

		return eleInfo.setMatrix(this->getTangentStiff());



    else if (responseID == 3)

        return eleInfo.setMatrix(this->getMass());



    else if (responseID == 4)

        return eleInfo.setMatrix(this->getDamp());



	else if (responseID == 5) {



		// Loop over the integration points

		int cnt = 0;

		for (int i = 0; i < nintu; i++) {

			// Get material stress response

			const Vector &sigma = materialPointers[i]->getStress();

			stresses(cnt++) = sigma(0);

			stresses(cnt++) = sigma(1);

			stresses(cnt++) = sigma(2);

			stresses(cnt++) = sigma(3);

			stresses(cnt++) = sigma(4);

			stresses(cnt++) = sigma(5);

		}

		return eleInfo.setVector(stresses);



	}

	else



		return -1;

}



// calculate local shape functions

void

TwentyEightNodeBrickUP::compuLocalShapeFunction() {



	int i, k, j;

	static double shl[4][20][27], w[27];

	// solid phase

	brcshl(shl, w, nintu, nenu);

	for(k = 0; k < nintu; k++) {

		wu[k] = w[k];

		for( j = 0; j < nenu; j++)

			for( i = 0; i < 4; i++)

				shlu[i][j][k] = shl[i][j][k];

	}

	// fluid phase

	brcshl(shl, w, nintp, nenu);

	for(k = 0; k < nintp; k++) {

		wp[k] = w[k];

		for( j = 0; j < nenu; j++)

			for( i = 0; i < 4; i++)

				shlq[i][j][k] = shl[i][j][k];

	}

	// coupling term

	brcshl(shl, w, nintp, nenp);

	for(k = 0; k < nintp; k++) {

		wp[k] = w[k];

		for( j = 0; j < nenp; j++)

			for( i = 0; i < 4; i++)

				shlp[i][j][k] = shl[i][j][k];

	}



}



void

TwentyEightNodeBrickUP::Jacobian3d(int gaussPoint, double& xsj, int mode)

{

	int i, j, k, nint, nen;

	double rxsj, c1, c2, c3;

	static double xs[3][3];

	static double ad[3][3];

	static double shp[4][20];



	if( mode == 0 ) { // solid

		nint = nintu;

		nen = nenu;

	}

	else if( mode == 1) { // fluid

		nint = nintp;

		nen = nenp;

	}

	else if( mode == 2 ) { // coupling

		nint = nintp;

		nen = nenu;

	}

	else {

		opserr <<"TwentyEightNodeBrickUP::Jacobian3d - illegal mode: " << mode << "\n";

		exit(-1);

	} //end if



	for( j = 0; j < nen; j++) {

		for( i = 0; i < 4; i++) {

			if( mode == 0 )

				shp[i][j] = shlu[i][j][gaussPoint];

			else if( mode == 1 )

				shp[i][j] = shlp[i][j][gaussPoint];

			else if( mode == 2 )

				shp[i][j] = shlq[i][j][gaussPoint];

			else {

				opserr <<"TwentyEightNodeBrickUP::Jacobian3d - illegal mode: " << mode << "\n";

				exit(-1);

			} //end if

		}

	}









	//Compute jacobian transformation



	for ( j=0; j<3; j++ ) {

		for( k = 0; k < 3; k++ ) {

			xs[j][k] = 0;

			for( i = 0; i < nen; i++ ) {

				xs[j][k] += xl[j][i] * shp[k][i];

			}

		}

	}





	//Compute adjoint to jacobian



	ad[0][0] = xs[1][1]*xs[2][2] - xs[1][2]*xs[2][1] ;

	ad[0][1] = xs[2][1]*xs[0][2] - xs[2][2]*xs[0][1] ;

	ad[0][2] = xs[0][1]*xs[1][2] - xs[0][2]*xs[1][1] ;



	ad[1][0] = xs[1][2]*xs[2][0] - xs[1][0]*xs[2][2] ;

	ad[1][1] = xs[2][2]*xs[0][0] - xs[2][0]*xs[0][2] ;

	ad[1][2] = xs[0][2]*xs[1][0] - xs[0][0]*xs[1][2] ;



	ad[2][0] = xs[1][0]*xs[2][1] - xs[1][1]*xs[2][0] ;

	ad[2][1] = xs[2][0]*xs[0][1] - xs[2][1]*xs[0][0] ;

	ad[2][2] = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;



	//Compute determinant of jacobian



	xsj  = xs[0][0]*ad[0][0] + xs[0][1]*ad[1][0] + xs[0][2]*ad[2][0] ;

	if (xsj<=0) {

		opserr <<"TwentyEightNodeBrickUP::Jacobian3d - Non-positive Jacobian: " << xsj << "\n";

		for( i = 0; i < nen; i++ ) {

			printf("%5d %15.6e %15.6e %15.6e %15.6e\n", i,

				shp[0][i], shp[1][i], shp[2][i], shp[3][i]);

		}



		exit(-1);

	}



	rxsj = 1.0/xsj ;



	//Compute jacobian inverse



	for ( j=0; j<3; j++ ) {



        for ( i=0; i<3; i++ )

			xs[i][j] = ad[i][j]*rxsj ;



	} //end for j





	//Compute derivatives with repect to global coords.



	for ( k=0; k<nen; k++) {



		c1 = shp[0][k]*xs[0][0] + shp[1][k]*xs[1][0] + shp[2][k]*xs[2][0] ;

        c2 = shp[0][k]*xs[0][1] + shp[1][k]*xs[1][1] + shp[2][k]*xs[2][1] ;

        c3 = shp[0][k]*xs[0][2] + shp[1][k]*xs[1][2] + shp[2][k]*xs[2][2] ;



        shp[0][k] = c1 ;

        shp[1][k] = c2 ;

        shp[2][k] = c3 ;



	} //end for k



	for( j = 0; j < nen; j++) {

		for( i = 0; i < 4; i++) {

			if( mode == 0 )

				shgu[i][j][gaussPoint] = shp[i][j];

			else if( mode == 1 )

				shgp[i][j][gaussPoint] = shp[i][j];

			else if( mode == 2 )

				shgq[i][j][gaussPoint] = shp[i][j];

			else {

				opserr <<"TwentyEightNodeBrickUP::Jacobian3d - illegal mode: " << mode << "\n";

				exit(-1);

			} //end if

		}

	}





}





