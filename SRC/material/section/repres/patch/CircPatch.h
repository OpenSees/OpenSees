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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/CircPatch.h,v $
                                                                        
                                                                        
// File: CircPatch.h
// Written by Remo M. de Souza
// December 1998

#ifndef CircPatch_h 
#define CircPatch_h 

#include <Patch.h>
#include <Vector.h>

class Cell;
class Matrix;

class CircPatch: public Patch
{
  public:

    CircPatch();
    CircPatch(int materialID, int numSubdivCircunf, int numSubdivRadial,
              const Vector &centerPosition, double internRadius, 
              double externRadius, double initialAngle, double finalAngle);

    ~CircPatch();
    
    // edition functions

    void setMaterialID     (int materialID);
    void setDiscretization (int numSubdivCircunf, int numSubdivRadial);
    void setCenterPosition (const Vector & centerPosition);
    void setRadii          (double internRadius, double externRadius);
    void setAngles         (double initialAngle, double finalAngle);

    // reinforcing bar inquiring functions
    
    int     getMaterialID         (void) const; 
    int     getNumCells           (void) const;
    Cell  **getCells              (void) const;
    Patch  *getCopy               (void) const;

    void   getDiscretization      (int &numSubdivCircunf, int &numSubdivRadial) const;
    void   getRadii               (double &internRadius, double &externRadius)  const;
    void   getAngles              (double &initialAngle, double &finalAngle)    const;
    const  Vector &getCenterPosition (void) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, CircPatch &CircPatch);    
    
  protected:
    
  private:
    int    matID;
    int    nDivCirc, nDivRad;
    Vector centerPosit;
    double intRad, extRad;
    double initAng, finalAng;
};


#endif

 
