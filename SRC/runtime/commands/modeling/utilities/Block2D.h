//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of Block2D.
//
// Written: Ed Love
// Created: 07/01
//
#ifndef Block2D_h
#define Block2D_h

#include <math.h>
#include <Vector3D.h>
#include <ID.h> 
class Matrix;
class Vector;

class Block2D {

public:

  Block2D(int numx, int numy, 
	  const ID& nodeID, 
	  const Matrix& coorArray,
	  int numNodeElement);


  // generate nodes
  Vector3D getNodalCoords(int i, int j);

  // generate elements
  const ID &getElementNodes(int i, int j);


 private:
  // set up xl array
  int setUpXl(const ID &nodeID, const Matrix &coorArray);
  
  // transform to real coordiantes
  void transformNodalCoordinates(Vector3D&);

  // shape functions
  void shape2d(double x1, 
                double x2, 
                double shape[9]);

  int nx; // number of elements x-direction
  int ny; // number of elements y-direction

  double xl[3][9]; //block coordinates 

  ID element; //ID-array of an element


  int numNodesElement; // 4 or 9
  int errorFlag;       // flag indicating if odd nx and ny ok for 9-noded elements
};

#endif
