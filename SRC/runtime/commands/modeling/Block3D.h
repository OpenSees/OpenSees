//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: Ed Love
// Created: 07/01
//
// Description: This file contains the implementation of Block3D.
//
#ifndef Block3D_h
#define Block3D_h
//

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 

class Block3D {

 public:

  //constructor
  Block3D(int numx, int numy, int numz,
	  const ID& nodeID, 
	  const Matrix& coorArray);

  //destructor
  virtual ~Block3D();

  //generate node 
  const Vector &getNodalCoords(int i, int j, int k);

  //generate element
  const ID &getElementNodes(int i, int j, int k);

 protected:

 private:

  int nx; //number of elements x-direction

  int ny; //number of elements y-direction

  int nz; //number of elements z-direction

  double xl[3][27]; //block coordinates 

  Vector coor; //coordinates of a node

  ID element; //ID-array of an element

  //set up xl array
  void setUpXl(const ID &nodeID, const Matrix &coorArray);
  
  //transform to real coordiantes
  void transformNodalCoordinates();

  //shape functions
  void shape3d(double x1, 
	       double x2, 
	       double x3,
	       double shape[27]);

};

#endif
