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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-05-19 06:00:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/R3vectors.h,v $

// Ed "C++" Love

#include <Vector.h>
#include <Matrix.h>
#include <math.h>

double  LovelyInnerProduct( const Vector &v, const Vector &w ) ;

double  LovelyNorm( const Vector &v ) ; 

Vector  LovelyCrossProduct( const Vector &v, const Vector &w ) ;

Vector  LovelyEig( const Matrix &M ) ;
