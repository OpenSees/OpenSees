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
                                                                        
// $Revision: 1.2 $
// $Date: 2001-10-02 20:20:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/GaussQuadRule1d.h,v $
                                                                        
// Written: rms
// Created: 12/98
//
// Description: This file contains the class definition for 
// GaussQuadRule1d (Quadrature Rule).

#ifndef GaussQuadRule1d_h
#define GaussQuadRule1d_h

#include <QuadRule1d.h>

class Vector;
class Matrix;

class GaussQuadRule1d: public QuadRule1d
{
  public:
    GaussQuadRule1d ();
    ~GaussQuadRule1d();

    int            setOrder              (int quadOrder);
    int            getOrder              (void) const;
    int            getNumIntegrPoints    (void) const;
    const Matrix & getIntegrPointCoords  (void) const;
    const Vector & getIntegrPointWeights (void) const; 
    const Matrix & getIntegrPointCoords  (int quadOrder);
    const Vector & getIntegrPointWeights (int quadOrder); 
    
  protected:
    
  private:
    int order;

    Matrix *myPts;
    Vector *myWts;

    enum {maxOrder = 10};

    static bool dataSet;

    static double ptsArray[];
    static double wtsArray[];

    static Matrix pts1;
    static Matrix pts2;
    static Matrix pts3;
    static Matrix pts4;
    static Matrix pts5;
    static Matrix pts6;
    static Matrix pts7;
    static Matrix pts8;
    static Matrix pts9;
    static Matrix pts10;

    static Vector wts1;
    static Vector wts2;
    static Vector wts3;
    static Vector wts4;
    static Vector wts5;
    static Vector wts6;
    static Vector wts7;
    static Vector wts8;
    static Vector wts9;
    static Vector wts10;
};

#endif
