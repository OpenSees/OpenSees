/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
 

#if !defined(AFX_SURFACEDESIGN_H)
#define AFX_SURFACEDESIGN_H 

#include  "Vector.h"
#include "PrincipalAxis.h"
#include "GridPlane.h" 

class SurfaceDesign  
{
public:
  virtual char * getType()=0;
  SurfaceDesign();
  virtual ~SurfaceDesign();
  
  virtual void setPincipalAxesPtr(PrincipalAxis ** pAxis){};
  virtual void setGridPlanesPtr( GridPlane ** pGridPlane){};
  
  //	virtual void setGradient( Vector * gradient){};
  
  virtual int fitCurve() ;
  virtual double getFunctionValue(Vector * Point);
  virtual double getFunctionValue2(Vector * point, Vector * dp2prime, Vector * gradG2);
  virtual double debug(Vector * point, Vector * dp2prime, Vector * gradG2);
  
 protected: 
  Matrix * theGridValues;
};

#endif // !defined(AFX_SURFACEDESIGN_H__0948F6C0_0502_490A_B32B_A46570768425__INCLUDED_)
