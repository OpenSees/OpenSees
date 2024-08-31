//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of Block3D.
//
// Written: Ed Love
// Created: 07/01
//
#include <Block3D.h>


//constructor
Block3D::Block3D(int numx, int numy, int numz,
                 const ID& nodeID, 
                 const Matrix& coorArray ) 
: nx(numx), ny(numy), nz(numz),
  coor(3), 
  element(8) 
{
  this->setUpXl( nodeID, coorArray );
}


// destructor
Block3D::~Block3D( )
{ 
}


//set up xl array
void  Block3D::setUpXl( const ID &nodeID, const Matrix &coorArray ) 
{

  for (int i=0; i<8; i++ ){
    if ( nodeID(i) == -1 ) {
      opserr << "Warning : in Block3D, block node " 
           << i 
           << " is not defined.  No Generation will take place."
           << endln;
      break; 
    }
  }


  // xl = tranpose coorArray(27,3) 
  for (int i=0; i<3; i++ ) {
    for (int j=0; j<27; j++ )
      xl[i][j] = coorArray(j,i);
  }


  if ( nodeID(8) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][8] = 0.5*( xl[i][0] + xl[i][4] );
  }

  if ( nodeID(9) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][9] = 0.5*( xl[i][1] + xl[i][5] );
  }

  if ( nodeID(10) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][10] = 0.5*( xl[i][2] + xl[i][6] );
  }

  if ( nodeID(11) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][11] = 0.5*( xl[i][3] + xl[i][7] );
  }


  if ( nodeID(12) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][12] = 0.5*( xl[i][0] + xl[i][1] );
  }

  if ( nodeID(13) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][13] = 0.5*( xl[i][1] + xl[i][2] );
  }

  if ( nodeID(14) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][14] = 0.5*( xl[i][2] + xl[i][3] );
  }

  if ( nodeID(15) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][15] = 0.5*( xl[i][0] + xl[i][3] );
  }


  if ( nodeID(16) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][16] = 0.25*( xl[i][0] + xl[i][1] + xl[i][2] + xl[i][3] );
  }


  if ( nodeID(17) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][17] = 0.5*( xl[i][4] + xl[i][5] );
  }

  if ( nodeID(18) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][18] = 0.5*( xl[i][5] + xl[i][6] );
  }

  if ( nodeID(19) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][19] = 0.5*( xl[i][6] + xl[i][7] );
  }

  if ( nodeID(20) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][20] = 0.5*( xl[i][4] + xl[i][7] );
  }


  if ( nodeID(21) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][21] = 0.25*( xl[i][4] + xl[i][5] + xl[i][6] + xl[i][7] );
  }


  if ( nodeID(22) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][22] = 0.25*( xl[i][0] + xl[i][1] + xl[i][5] + xl[i][4] );
  }

  if ( nodeID(23) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][23] = 0.25*( xl[i][1] + xl[i][2] + xl[i][6] + xl[i][5] );
  }


  if ( nodeID(24) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][24] = 0.25*( xl[i][3] + xl[i][2] + xl[i][6] + xl[i][7] );
  }

  if ( nodeID(25) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][25] = 0.25*( xl[i][0] + xl[i][3] + xl[i][7] + xl[i][4] );
  }



  if ( nodeID(26) == -1 ) {
    for (int i=0; i<3; i++ )
      xl[i][26] = 0.125*( xl[i][0] + xl[i][1] + xl[i][2] + xl[i][3] +
                          xl[i][4] + xl[i][5] + xl[i][6] + xl[i][7]   );
  }

  return;
}


//generate node
const Vector&
Block3D::getNodalCoords( int i, int j, int k )
{

  /* loop as follows (in pseudocode)
     for ( k = 0, nz ) {
       for ( j = 0, ny ) {
         for ( i = 0, nx ) 
           call getNodalCoords(i,j,k);
       } 
     }
  */

  double hx = 2.0 / nx;

  double hy = 2.0 / ny;

  double hhz = 2.0 / nz;

  double x = -1.0 + (i*hx);

  double y = -1.0 + (j*hy);

  double z = -1.0 + (k*hhz);

  coor(0) = x;
  coor(1) = y;
  coor(2) = z;

  this->transformNodalCoordinates( );
  
  return coor;
}


//generate element
const ID&
Block3D::getElementNodes( int i, int j, int k )  
{

  int nenx = nx + 1;

  int neny = ny + 1;

  int nInXYplane = nenx * neny;


  int node1, node2, node3, node4;
  int node5, node6, node7, node8;


  node1 = i + (j*nenx) + (k*nInXYplane);
  node2 = node1 + 1;
  node3 = node2 + nenx;
  node4 = node1 + nenx;

  node5 = node1 + nInXYplane;
  node6 = node2 + nInXYplane;
  node7 = node3 + nInXYplane;
  node8 = node4 + nInXYplane;

  element(0) = node1;
  element(1) = node2;
  element(2) = node3;
  element(3) = node4;

  element(4) = node5;
  element(5) = node6;
  element(6) = node7;
  element(7) = node8;

  return element;
}



//transform to real coordinates
void  Block3D::transformNodalCoordinates( )
{

  double shape[27]; 
  static double natCoor[3];

  natCoor[0] = coor(0);
  natCoor[1] = coor(1);
  natCoor[2] = coor(2);

  coor.Zero( );

  this->shape3d( natCoor[0], natCoor[1], natCoor[2], shape );

  for (int j=0; j<27; j++ ) {
      
    for (int dim=0; dim<3; dim++ )
      coor(dim) += shape[j]*xl[dim][j];

  }

  return;

}



//shape functions
void  Block3D::shape3d( double r, double s, double t,
                        double shape[27]     ) 
/*
 * Adapted from:
      subroutine shp04(shp,glu,glo,gu,eu,to,xjac,detj,r,s,t,xl,ul)
c-----------------------------------------------------------------------
c.....compute shape functions and their derivatives for linear,quadratic
c.....lagrangian and serendipity isoparametric  3-d elements
c.....global coordinate system x,y,z
c.....local coordinate system xsi,eta,zeta
c-----------------------------------------------------------------------
*/
{

  static constexpr int ri[] = {-1, 1, 1,-1, -1, 1, 1,-1,  -1, 1, 1,-1,   0, 1, 0,-1, 0, 0, 1, 0,-1, 0,   0, 1, 0,-1, 0};

  static constexpr int si[] = {-1,-1, 1, 1, -1,-1, 1, 1,  -1,-1, 1, 1,  -1, 0, 1, 0, 0, -1, 0, 1, 0, 0,  -1, 0, 1, 0, 0};

  static constexpr int ti[] = {-1,-1,-1,-1,  1, 1, 1, 1,   0, 0, 0, 0,  -1,-1,-1,-1,-1, 1, 1, 1, 1, 1,   0, 0, 0, 0, 0};


  static const double d1 = 1.0;
  static const double d2 = 0.5;   // = 1.0/2.0;
  static const double d4 = 0.25;  // = 1.0/4.0;
  static const double d8 = 0.125; // = 1.0/8.0;


  double rr = r*r;
  double ss = s*s;
  double tt = t*t;

  int kk;

  //shape functions for 27-node element
  for (int k=1; k<=27; k++ ) {

    kk = k-1; //C-style numbering 

    double r0 = r*ri[kk];
    double s0 = s*si[kk];
    double t0 = t*ti[kk];

    //corner nodes top/bottom
    if ( k>=1 && k<=8 ) 
       shape[kk] = d8*(rr+r0)     *(ss+s0)          *(tt+t0);

    //corner nodes midside
    if ( k>=9 && k<=12 )  
       shape[kk] = d4*(rr+r0)     *(ss+s0)     *(d1-tt);

    //midside nodes top/bottom  r-dir
    if ( k==13 ||  k==15 || k==18 || k==20 )
       shape[kk] = d4*(d1-rr)*(ss+s0)     *(tt+t0);

    //midside nodes top/bottom  s-dir
    if ( k==14 ||  k==16 ||  k==19 ||  k==21 )
      shape[kk] = d4*(rr+r0)     *(d1-ss)*(tt+t0);

    //midside nodes mid plane  r-dir
    if ( k==23 ||  k==25 )
      shape[kk] = d2*(d1-rr)*(ss+s0)     *(d1-tt);

    //midside nodes mid plane  s-dir
    if ( k==24 ||  k==26 ) 
      shape[kk] = d2*(rr+r0)     *(d1-ss)*(d1-tt);

    //central nodes top/bottom
    if ( k==17 ||  k==22 )
      shape[kk] = d2*(d1-rr)*(d1-ss)*(tt+t0);

    if ( k==27 )
      shape[kk] = (d1-rr)*(d1-ss)*(d1-tt);

  }

  return;
}

