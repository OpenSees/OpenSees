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
// $Date: 2007-02-16 00:24:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/PlasticDamageMaterial.h,v $
                                                                      
// Written: Jeeho Lee
// Created: Feb 2007
//
// Description: This file contains the class definition for 
// PlasticDamageMaterial. PlasticDamageMaterial wraps the FEDEAS
// 1d Material subroutine: PD_1.f

/* *******************************************************************************

References:

1. Lee, J., and Fenves, G. L. (2007). "Plastic-Damage Concrete Model for Beam-Column
Section Elements Subjected to Cyclic Loading", J. of Structural Engineering,
submitted.

2. Lee, J., and Fenves, G. L. (1998). "Plastic-damage model for cyclic loading of
concrete structures.Ó J. of Engineering Mechanics, 124(8), 892-900.

3. Lee, J., and Fenves, G. L. (2001). "A return-mapping algorithm for plastic-damage
models: 3-D and Plane Stress Formulation" Int. J. Numerical Methods in Eng., 50,
487-506.

******************************************************************************** */


#ifndef PlasticDamageMaterial_h
#define PlasticDamageMaterial_h

#include <FedeasMaterial.h>

class PlasticDamageMaterial : public FedeasMaterial
{
  public:
    PlasticDamageMaterial(int tag, double E, double Ft, double Fc, double ft_max,
			  double fcy, double fc_max, double kt_crit, double Relax);
    PlasticDamageMaterial(int tag, const Vector &d);
    PlasticDamageMaterial(void);
    virtual ~PlasticDamageMaterial();
    
    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

    int invokeSubroutine(int ist);
    
 protected:
    
 private:
    
};


#endif

