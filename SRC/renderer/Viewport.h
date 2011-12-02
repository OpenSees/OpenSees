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
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Viewport.h,v $
                                                                        
                                                                        
#ifndef Viewport_H
#define Viewport_H

#include "db.H"
class Device;

class Viewport {
 public:
  Viewport();
  ~Viewport();
  
  int update(void);
  FACE &transform(FACE &);
  MYPOINT *transformP(MYPOINT *input);  

  void setDevice(Device &theDevice);  
  VECTOR portwindow;

 private:
  MATRIX TMat;
  Device *theDevice;
};

#endif
