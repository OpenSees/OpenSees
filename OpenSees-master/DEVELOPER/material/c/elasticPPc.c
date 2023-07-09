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

/*                                                                        
** $Revision: 1.8 $
** $Date: 2009/01/13 01:00:52 $
** $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/c/elasticPPC.c,v $
                                                                        
** Written: fmk 
**
** Description: This file contains the implementation of elasticPP material
*/

#include "elementAPI.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define DBL_EPSILON 1e-18

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C" 
#endif

OPS_Export void
elasticPPc (matObj *thisObj, 
	    modelState *model, 
	    double *strain, 
	    double *tang, 
	    double *stress, 
	    int *isw, 
	    int *result) 
{
  *result = 0;

  if (*isw == ISW_INIT) {

    double dData[2];
    int    iData[1];

    /* get the input data  - tag? E? eyp? */
    int numData = 1;
    OPS_GetIntInput(&numData, iData);
    numData = 2;
    OPS_GetDoubleInput(&numData, dData); 

    /* Allocate the element state */
    thisObj->tag = iData[0];
    thisObj->nParam = 2;  /* E, eyp */
    thisObj->nState = 2;  /* strain, ep */
    OPS_AllocateMaterial(thisObj);

    thisObj->theParam[0] = dData[0];
    thisObj->theParam[1] = dData[1];

    *result = OPS_UNIAXIAL_MATERIAL_TYPE;

  } else if (*isw == ISW_COMMIT) {

    double trialStrain = thisObj->tState[0];
    double f, fYieldSurface;		// yield function

    double E = thisObj->theParam[0];
    double eyp = thisObj->theParam[1];
    double ep = thisObj->cState[1];

    double fyp = E*eyp;
    double fyn = -fyp;

    // compute trial stress
    double sigtrial = E * ( trialStrain - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    fYieldSurface = - E * DBL_EPSILON;
    if ( f > fYieldSurface ) {
      // plastic
      if ( sigtrial > 0.0 ) {
	ep += f / E;
      } else {
	ep -= f / E;
      }
    }

    thisObj->cState[0] = trialStrain;    
    thisObj->cState[1] = ep;

  } else if (*isw == ISW_REVERT) {
    int i;
    for (i=0; i<3; i++)
      thisObj->tState[i] = thisObj->cState[i];

  } else if (*isw == ISW_REVERT_TO_START) {
	int i;
    for (i=0; i<2; i++) {
      thisObj->cState[i] = 0.0;
      thisObj->tState[i] = 0.0;
    }
    
  } else if (*isw == ISW_FORM_TANG_AND_RESID) {

    double trialStrain = *strain;
    double f;		// yield function
    double trialStress, trialTangent, fYieldSurface;

    double E = thisObj->theParam[0]; 
    double eyp = thisObj->theParam[1];
    double ep = thisObj->cState[1];

    double fyp = E*eyp;
    double fyn = -fyp;
    
    // compute trial stress
    double sigtrial = E * ( trialStrain - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    fYieldSurface = - E * DBL_EPSILON;
    if ( f <= fYieldSurface ) {

      // elastic
      trialStress = sigtrial;
      trialTangent = E;

    } else {

      // plastic
      if ( sigtrial > 0.0 ) {
	trialStress = fyp;
      } else {
	trialStress = fyn;
      }

      trialTangent = 0.0;
    }

    thisObj->tState[0] = trialStrain;
    *stress = trialStress;
    *tang = trialTangent;
  }

}


OPS_Export void
localInit() 
{
  OPS_Error("elasticPPC uniaxial material - Written by fmk UC Berkeley Copyright 2008 - Use at your Own Peril\n", 1);
}
