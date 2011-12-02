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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/QuadPatch.h,v $
                                                                        
                                                                        
// File: QuadPatch.h
// Written by Remo M. de Souza
// December 1998

#ifndef QuadPatch_h 
#define QuadPatch_h 

#include <Patch.h>

class Cell;
class Matrix;

class QuadPatch: public Patch
{
  public:

    QuadPatch();
    QuadPatch(int materialID, int numSubdivIJ, int numSubdivJK,
              const Matrix &vertexCoords);
        
    ~QuadPatch();
    
    // edition functions

    void setMaterialID     (int materialID);
    void setDiscretization (int numSubdivIJ, int numSubdivJK);
    void setVertCoords     (const Matrix &vertexCoords);

    // reinforcing bar inquiring functions
    
    int     getMaterialID         (void) const; 
    int     getNumCells           (void) const;
    Cell  **getCells              (void) const;
    Patch  *getCopy               (void) const;

    void   getDiscretization     (int &numSubdivIJ, int &numSubdivJK) const;
    const  Matrix &getVertCoords (void) const;

    void Print(ostream &s, int flag =0) const;   
    friend ostream &operator<<(ostream &s, QuadPatch &quadPatch);    
    
  protected:
    
  private:
    int    matID;
    int    nDivIJ, nDivJK;
    Matrix vertCoord;
};


#endif

 
