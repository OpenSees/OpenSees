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
// File: ReinfLayer.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef ReinfLayer_h 
#define ReinfLayer_h 

// #include <OPS_Globals.h>

class ReinfBar;
class OPS_Stream;

class ReinfLayer
{
  public:
 
    // edition functions
    virtual void setNumReinfBars     (int numReinfBars)        = 0;
    virtual void setMaterialID       (int materialID)          = 0;
    virtual void setReinfBarDiameter (double reinfBarDiemater) = 0;
    virtual void setReinfBarArea     (double reinfBarArea)     = 0;

    // reinforcing layer inquiring functions
    
    virtual int         getNumReinfBars() const = 0;
    virtual int         getMaterialID() const = 0; 
    virtual double      getReinfBarDiameter() const = 0;
    virtual double      getReinfBarArea() const = 0;
    virtual ReinfLayer *getCopy() const = 0;
    virtual ReinfBar   *getReinfBars() const = 0;     
   
    virtual void Print(OPS_Stream &s, int flag =0) const = 0;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const ReinfLayer &ReinfLayer);    
    
  protected:
    
  private:
};


#endif
