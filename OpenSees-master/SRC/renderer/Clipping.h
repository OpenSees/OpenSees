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
// $Date: 2000-09-15 08:23:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Clipping.h,v $
                                                                        
                                                                        
#ifndef Clipping_H
#define Clipping_H

#include "db.H"

class Clipping {
 public:
  Clipping();
  ~Clipping();  

  int update(void);
  FACE &transform(FACE &);
  MYPOINT *transformP(MYPOINT *input);
      
 private:  
  float X, Y, Zfar, Znear;
};

#endif
