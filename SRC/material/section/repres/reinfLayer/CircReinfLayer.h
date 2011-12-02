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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/reinfLayer/CircReinfLayer.h,v $
                                                                        
                                                                        
// File: CircReinfLayer.h
// Written by Remo M. de Souza
// December 1998


#ifndef CircReinfLayer_h 
#define CircReinfLayer_h 

#include <ReinfLayer.h>

class ReinfBar;

class CircReinfLayer : public ReinfLayer
{
  public:

    CircReinfLayer();
    CircReinfLayer(int materialID, int numReinfBars, double  reinfBarArea,
                   const Vector &centerPosition, double arcRadius, double
                   initialAngle, double finalAngle);

    ~CircReinfLayer();
    
    // edition functions

    void setNumReinfBars     (int numReinfBars);
    void setMaterialID       (int materialID);  
    void setReinfBarDiameter (double reinfBarDiameter);
    void setReinfBarArea     (double reinfBarArea);

    void setCenterPosition   (const Vector &centerPosition);
    void setArcRadius        (double arcRadius);
    void setInitAngle        (double initialAngle);
    void setFinalAngle       (double finalAngle);  

    // inquiring functions

    int           getNumReinfBars     (void) const;
    int           getMaterialID       (void) const;
    int           getReinfBarDiameter (void) const;
    double        getReinfBarArea     (void) const;
    ReinfBar     *getReinfBars        (void) const;

  
    ReinfLayer   *getCopy             (void) const;
    const Vector &getCenterPosition   (void) const;
    double        getArcRadius        (void) const;    
    double        getInitAngle        (void) const;
    double        getFinalAngle       (void) const;  
    

    void Print(ostream &s, int flag =0) const;   
    friend ostream &operator<<(ostream &s, const CircReinfLayer &CircReinfLayer);    
    
  protected:
    
  private:
    int    nReinfBars;
    int    matID;
    double barDiam;
    double area;
    Vector centerPosit;
    double arcRad;
    double initAng;
    double finalAng;
};


#endif




