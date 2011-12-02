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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/reinfLayer/ReinfLayer.h,v $
                                                                        
                                                                        
// File: ReinfLayer.h
// Written by Remo M. de Souza
// December 1998

#ifndef ReinfLayer_h 
#define ReinfLayer_h 

#include <iostream.h>

class ReinfBar;

class ReinfLayer
{
  public:

    ReinfLayer();
    virtual ~ReinfLayer();
    
    // edition functions

    virtual void setNumReinfBars     (int numReinfBars)        = 0;
    virtual void setMaterialID       (int materialID)          = 0;
    virtual void setReinfBarDiameter (double reinfBarDiemater) = 0;
    virtual void setReinfBarArea     (double reinfBarArea)     = 0;

    // reinforcing layer inquiring functions
    
    virtual int         getNumReinfBars     (void) const = 0;
    virtual int         getMaterialID       (void) const = 0; 
    virtual int         getReinfBarDiameter (void) const = 0;
    virtual double      getReinfBarArea     (void) const = 0;
    virtual ReinfLayer *getCopy             (void) const = 0;
    virtual ReinfBar   *getReinfBars        (void) const = 0;     
   
    virtual void Print(ostream &s, int flag =0) const = 0;   
    friend ostream &operator<<(ostream &s, const ReinfLayer &ReinfLayer);    
    
  protected:
    
  private:
};


#endif

