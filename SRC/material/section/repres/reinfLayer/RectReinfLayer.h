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
                                                                        
#ifndef RectReinfLayer_h 
#define RectReinfLayer_h 

#include <ReinfLayer.h>
#include <Vector.h>

class ReinfBar;

class RectReinfLayer : public ReinfLayer
{
public:
  
  RectReinfLayer();
  RectReinfLayer(int materialID, int numy, int numz, double reinfBarArea,
		 double yc, double zc, double dy, double dz);
  
  ~RectReinfLayer();
  
  // edition functions
  
  void setNumReinfBars     (int numReinfBars);
  void setMaterialID       (int materialID);  
  void setReinfBarDiameter (double reinfBarDiameter);
  void setReinfBarArea     (double reinfBarArea);
  
  void setCenterPosition   (const Vector &centerPosition);
  
  // inquiring functions
  
  int           getNumReinfBars     (void) const;
  int           getMaterialID       (void) const;
  double        getReinfBarDiameter (void) const;
  double        getReinfBarArea     (void) const;
  ReinfBar     *getReinfBars        (void) const;
  
  
  ReinfLayer   *getCopy             (void) const;
  
  void Print(OPS_Stream &s, int flag =0) const;   
  friend OPS_Stream &operator<<(OPS_Stream &s, const RectReinfLayer &CircReinfLayer);
  
protected:
  
private:
  int    nBarsy;
  int    nBarsz;
  int    matID;
  double barDiam;
  double area;
  double yc, zc; // center coordinates
  double dy, dz; // "diameter" in y,z
};


#endif




