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
                                                                        
#ifndef CircSectionCell_h 
#define CircSectionCell_h 

#include <Cell.h>
#include <Vector.h>

class Matrix;
class Vector;


class CircSectionCell: public Cell
{
  public:

    CircSectionCell();
    CircSectionCell(double r2, double r1, double alpha, double theta, double centerX, double centerY);
        
    ~CircSectionCell();
    
    // edition functions

    void setVertCoords (const Matrix &vertexCoords);

    // reinforcing bar inquiring functions
    
    double getArea                     (void) const;
    double getdValue                   (void) const;    
    const  Matrix &getVertCoords       (void) const;
    const  Vector &getCentroidPosition (void);

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const CircSectionCell &quadCell);    
    
  protected:
    
  private:
    double r1, r2; // r1 inner and r2 outer radii
    double alpha;  // inner angle of section
    double theta;  // angle of centerline about z axis

    double A;
    Vector Centroid;
    double offsetX, offsetY;
//    double area;
};


#endif

