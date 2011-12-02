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
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/GaussLobattoQuadRule1d01.h,v $
                                                                        
                                                                        
// File: ~/QuadRule/GaussLobattoQuadRule1d01.h
//
// Written: rms
// Created: 12/98
// Revision: 
//
// Description: This file contains the class definition for 
// GaussLobattoQuadRule1d01 (Quadrature Rule).
//
// What: "@(#) GaussLobattoQuadRule1d01.h, revA"


#ifndef GaussLobattoQuadRule1d01_h
#define GaussLobattoQuadRule1d01_h

#include <QuadRule1d01.h>

class Vector;
class Matrix;

class GaussLobattoQuadRule1d01: public QuadRule1d01
{
  public:
    GaussLobattoQuadRule1d01 (int quadOrder);
    ~GaussLobattoQuadRule1d01();

    int            setOrder              (int quadOrder);
    int            getOrder              (void) const;
    int            getNumIntegrPoints    (void) const;
    const Matrix & getIntegrPointCoords  (void) const;
    const Vector & getIntegrPointWeights (void) const; 
    
  protected:
    
  private:
    int order;
    Matrix *coord;
    Vector *weight;
};


#endif



