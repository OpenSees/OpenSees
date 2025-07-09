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
// File: StraightReinfLayer.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef StraightReinfLayer_h 
#define StraightReinfLayer_h 

#include <ReinfLayer.h>
#include <Vector.h>

class ReinfBar;

class StraightReinfLayer : public ReinfLayer
{
  public:

    StraightReinfLayer(int materialID, int numReinfBars, double  reinfBarArea,
                       const Vector &initialPosition, 
                       const Vector &finalPosition);

    ~StraightReinfLayer();
    
    // edition functions

    void setNumReinfBars     (int numReinfBars);
    void setMaterialID       (int materialID);  
    void setReinfBarDiameter (double reinfBarDiameter);
    void setReinfBarArea     (double reinfBarArea);

    void setInitialPosition (const Vector &initialPosition);
    void setFinalPosition   (const Vector &finalPosition);

    // inquiring functions

    int           getNumReinfBars     (void) const;
    int           getMaterialID       (void) const;
    double        getReinfBarDiameter (void) const;
    double        getReinfBarArea     (void) const;
    ReinfBar     *getReinfBars        (void) const;

  
    ReinfLayer   *getCopy             (void) const;
    const Vector &getInitialPosition  (void) const;
    const Vector &getFinalPosition    (void) const;

    void Print(OPS_Stream &s, int flag =0) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const StraightReinfLayer &straightReinfLayer);    
    
  protected:
    
  private:
    int    nReinfBars;
    int    matID;
    double barDiam;
    double area;
    Vector initPosit;
    Vector finalPosit;
};


#endif

