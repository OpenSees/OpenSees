//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for Block2D.
//
// Written: Ed Love
// Created: 07/99
//
#include <Vector.h>
#include <Matrix.h>
#include <Vector3D.h>
#include <Block2D.h>


Block2D::Block2D(int numx, int numy, 
                 const ID& nodeID, 
                 const Matrix& coorArray,
                 int numNode) 
: 
  element(numNode),
  numNodesElement(numNode),
  nx(numx), ny(numy),
  errorFlag(0)
{

  if (numNodesElement == 9) {
    if (((numx % 2) != 0) || ((numy % 2) != 0)) {
      opserr << "ERROR: Block2D::Block2D - numX & numY for nine noded elements must be even\n";
      errorFlag = 1;
    }
  } 

  if (numNodesElement != 9 && numNodesElement != 4) {
      opserr << "ERROR: Block2D::Block2D - numNode must be either 4 or 9\n";
      errorFlag = 1;
  }

  if (this->setUpXl( nodeID, coorArray ) != 0)
    errorFlag = 1;

}


//destructor
Block2D::~Block2D( )
{ 

}


//set up xl array
int
Block2D::setUpXl( const ID &nodeID, const Matrix &coorArray ) 
{


  for (int i=0; i<4; i++ ){
    if ( nodeID(i) == -1 ) {
      opserr << "Warning : in Block2D, block node " 
                   << i << " is not defined.\n";
      return -1;
    }
  }

  // local storage xl = transpose(coorArray)
  for (int i=0; i<3; i++ ) {
    for (int j=0; j<9; j++ )
      xl[i][j] = coorArray(j,i);
  }


  if ( nodeID(4) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][4] = 0.5*( xl[i][0] + xl[i][1] );
  }

  if ( nodeID(5) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][5] = 0.5*( xl[i][1] + xl[i][2] );
  }

  if ( nodeID(6) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][6] = 0.5*( xl[i][2] + xl[i][3] );
  }

  if ( nodeID(7) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][7] = 0.5*( xl[i][3] + xl[i][0] );
  }

  if ( nodeID(8) == -1 ) {
    for (int i=0; i<3; i++ ) 
      xl[i][8]  = 0.25*( xl[i][0] + xl[i][1] + xl[i][2] + xl[i][3] ) ;
  }

  
  return 0;
}


//generate node
Vector3D
Block2D::getNodalCoords( int i, int j )
{

  double hx = 2.0 / nx;
  double hy = 2.0 / ny;
  double x = -1.0 + (i*hx);
  double y = -1.0 + (j*hy);

  Vector3D coor;
  coor[0] = x;
  coor[1] = y;
  coor[2] = 0.0;

  this->transformNodalCoordinates(coor);
  
  return coor;
}


//generate element
const ID&
Block2D::getElementNodes( int i, int j )  
{

  if (errorFlag == 1)
    return element;

  else if (numNodesElement == 4) {
    int nenx = nx + 1;
    // int neny = ny + 1;


    int node1 = i + j*nenx;
    int node2 = node1 + 1;

    int node3 = node2 + nenx;
    int node4 = node1 + nenx;
    
    element(0) = node1;
    element(1) = node2;
    element(2) = node3;
    element(3) = node4;

  } else {

    int nenx = nx + 1;
    // int neny = ny + 1;

    int node1, node2, node3, node4, node5, node6, node7, node8, node9;

    node1 = i*2 + j*2*nenx;
    node5 = node1 + 1;
    node2 = node1 + 2;

    node4 = node1 + 2*nenx;
    node7 = node4 + 1;
    node3 = node7 + 1;

    node8 = node1 + nenx;
    node9 = node8 + 1;
    node6 = node9 + 1;
    
    element(0) = node1;
    element(1) = node2;
    element(2) = node3;
    element(3) = node4;
    element(4) = node5;
    element(5) = node6;
    element(6) = node7;
    element(7) = node8;
    element(8) = node9;

  }

  return element;
}



//transform to real coordiantes
void
Block2D::transformNodalCoordinates(Vector3D& coor )
{

  static double shape[9]; 
  static double natCoor[2];

  natCoor[0] = coor[0];
  natCoor[1] = coor[1];

  coor.zero( );

  this->shape2d( natCoor[0], natCoor[1], shape );

  for ( int j=0; j<9; j++ ) {
      
    for ( int dim=0; dim<3; dim++ )
      coor[dim] += shape[j]*xl[dim][j];

  }

  return;

}


// shape functions
void  Block2D::shape2d( double x, double y, 
                        double shape[9]     ) 
{
  static double Nx[3];
  static double Ny[3];

  Nx[0] = 0.5 * x * ( x - 1.0 );
  Nx[1] = 1.0 - (x*x);
  Nx[2] = 0.5 * x * ( x + 1.0 );

  Ny[0] = 0.5 * y * ( y - 1.0 );
  Ny[1] = 1.0 - (y*y);
  Ny[2] = 0.5 * y * ( y + 1.0 );

  shape[0] = Nx[0]*Ny[0];
  shape[1] = Nx[2]*Ny[0];
  shape[2] = Nx[2]*Ny[2];
  shape[3] = Nx[0]*Ny[2];

  shape[4] = Nx[1]*Ny[0];
  shape[5] = Nx[2]*Ny[1];
  shape[6] = Nx[1]*Ny[2];
  shape[7] = Nx[0]*Ny[1];

  shape[8] = Nx[1]*Ny[1];

  return;
}
