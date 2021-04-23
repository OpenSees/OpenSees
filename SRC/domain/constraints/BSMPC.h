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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/BeamSolidMPC2D.h,v $
                                                                        
                                                                        
// File: ~/model/constraints/BeamSolidMPC2D.h
//
// Written: Ziping Zhu,zzpxbkarl@gmail.com
// Revised:
//
// Purpose: This file contains the class definition for BS3d which interacts element
// of beam and solid in 3D.
// BSMPC  $rNodeTag1  $rNodeTag2  $cNodeTag1  $cNodeTag2…$cNodeTag_n $v1 $v2 $v3
// $rNodeTag1 is a node of Beam-column element in multiscale interface.
// $rNodeTag2 is a node of Beam-column element which is not in multiscale interface.
// These two nodes can determine the local coordination x'.
// $cNodeTag1  $cNodeTag2…$cNodeTag_n are nodes of solid element in multiscale interface.
// $v1 $v2 $v3 are the verctor to determine local coordination, which is not parallel to x'.
// local coordination: x'= $rNodeTag2 - $rNodeTag1, y'=v * x', z'=x' * y'


#ifndef BSMPC_h
#define BSMPC_h
#include <Matrix.h>
class Domain;
class ID;
class Matrix;

class BSMPC
{
  public:
    BSMPC(Domain &theDomain, int ndm, int nR1, int nR2, ID &nodeC,const Vector &vectorforlocal);
    virtual ~BSMPC();
    
  protected:
    int computeR2d(Domain &theDomain, int nR, ID &nC);
	int computeR3d(Domain &theDomain, int nR, ID &nC);
  private:
	double cosaerfa1, cosaerfa2, cosaerfa3, cosbeta1, cosbeta2, cosbeta3, cosgama1, cosgama2, cosgama3;
	static Matrix R;
};

#endif