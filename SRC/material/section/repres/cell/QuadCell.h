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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/cell/QuadCell.h,v $
                                                                        
                                                                        
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
    
    double getArea                     (void) const;
    const  Matrix &getVertCoords       (void) const;
    const  Vector &getCentroidPosition (void);

    void Print(ostream &s, int flag =0) const;   
    friend ostream &operator<<(ostream &s, const QuadCell &quadCell);    
    
  protected:
    
  private:
    Matrix vertCoord;
    Vector Centroid;
//    double area;
};


#endif

