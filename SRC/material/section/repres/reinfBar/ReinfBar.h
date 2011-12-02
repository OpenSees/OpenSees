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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/reinfBar/ReinfBar.h,v $
                                                                        
                                                                        
// File: ReinfBar.h
//
// Written by Remo M. de Souza
// November 1998


#ifndef ReinfBar_h 
#define ReinfBar_h 

#include <iostream.h>

class Vector;

class ReinfBar
{
  public:

    ReinfBar();
    ReinfBar(double barArea, int materialID, const Vector &position);
        
    virtual ~ReinfBar();
    
    // edition functions

    void setDiameter (double barDiameter);
    void setArea     (double barArea);
    void setMaterial (int materialID);
    void setPosition (const Vector &position);

    // reinforcing bar inquiring functions
    
    double getDiameter (void) const;
    double getArea     (void) const;
    int    getMaterial (void) const; 

    const Vector & getPosition (void) const;

    virtual void Print(ostream &s, int flag =0) const;   
    friend ostream &operator<<(ostream &s, const ReinfBar &reinfBar);    
    
  protected:
    
  private:
    int    matID;
    double diameter;
    double area;
    Vector posit;
};


#endif

