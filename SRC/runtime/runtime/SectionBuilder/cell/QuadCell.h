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
//
// File: QuadCell.h
//
// Written by Remo M. de Souza
// December 1998


#ifndef QuadCell_h 
#define QuadCell_h 

#include <Cell.h>
#include <Vector.h>

class Matrix;
class Vector;


class QuadCell: public Cell
{
  public:

    QuadCell();
    QuadCell(const Matrix &vertexCoords);
        
    ~QuadCell();
    
    // edition functions

    void setVertCoords (const Matrix &vertexCoords);

    // reinforcing bar inquiring functions
    
    double getArea() const;
    double getdValue() const;
    const  Matrix &getVertCoords() const;
    const  Vector &getCentroidPosition();

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const QuadCell &quadCell);    
    
  protected:
    
  private:
    Matrix vertCoord;
    Vector Centroid;
//    double area;
};


#endif

