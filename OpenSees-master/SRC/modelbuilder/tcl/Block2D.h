#ifndef Block2D_h
#define Block2D_h

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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/Block2D.h,v $
                                                                        
// Written: Ed Love
// Created: 07/01
//
// Description: This file contains the implementation of Block2D.

//
// What: "@(#) Block2D.h, revA"

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 

class Block2D {

 public:

  //constructor
  Block2D(int numx, int numy, 
	  const ID& nodeID, 
	  const Matrix& coorArray,
	  int numNodeElement);

  //destructor
  ~Block2D();

  //generate node 
  const Vector &getNodalCoords(int i, int j);

  //generate element
  const ID &getElementNodes(int i, int j);

 protected:

 private:

  int nx; //number of elements x-direction

  int ny; //number of elements y-direction

  double xl[3][9]; //block coordinates 

  Vector coor; //coordinates of a node

  ID element; //ID-array of an element

  //set up xl array
  void setUpXl(const ID &nodeID, const Matrix &coorArray);
  
  //transform to real coordinates
  void transformNodalCoordinates();

  //shape functions
  void shape2d(double x1, 
	       double x2, 
	       double shape[9]);

  int numNodesElement; // 4 or 9
  int errorFlag;       // flag indicating if odd nx and ny ok for 9-noded elements
};

#endif
