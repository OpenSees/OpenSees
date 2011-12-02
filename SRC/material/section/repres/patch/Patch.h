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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/Patch.h,v $
                                                                        
                                                                        
// File: Patch.h
// Written by Remo M. de Souza
// December 1998


#ifndef Patch_h 
#define Patch_h 


#include <OPS_Globals.h>

class Cell;

class Patch
{
  public:

    Patch();
    virtual ~Patch();
    
    // edition functions

    virtual void setMaterialID (int materialID) = 0;
    
    // inquiring functions

    virtual int     getMaterialID (void) const = 0; 
    virtual int     getNumCells   (void) const = 0;
    virtual Cell  **getCells      (void) const = 0;
    virtual Patch  *getCopy       (void) const = 0;

    virtual void Print(OPS_Stream &s, int flag =0) const =0;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const Patch &patch);    

  protected:
    
  private:
};


#endif

