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
** $Revision$
** $Date$
** $Source$
**                                                                        
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

static int initFlag = 0;

OPS_Export void
elasticPlaneStressC(matObj *thisObj, 
		    modelState *model, 
		    double *strain, 
		    double *tang, 
		    double *stress, 
		    int *isw, 
		    int *result) 
{
  *result = 0;

  if (*isw == ISW_INIT) {

    if (initFlag == 0) {
      OPS_Error("elasticPlaneStressC nd material - Written by fmk UC Berkeley Copyright 2008 - Use at your Own Peril\n", 1);
      initFlag = 1;
    }

    double dData[2];
    int    iData[1];

    /* get the input data  - tag? E? v? */
    int numData = 1;
    OPS_GetIntInput(&numData, iData);
    numData = 2;
    OPS_GetDoubleInput(&numData, dData); 

    /* Allocate the element state */
    thisObj->tag = iData[0];
    thisObj->nParam = 2;  /* E, v */
    thisObj->nState = 3;  /* strain */
    OPS_AllocateMaterial(thisObj);

    thisObj->theParam[0] = dData[0];
    thisObj->theParam[1] = dData[1];

    *result = OPS_PLANESTRESS_TYPE;

  } else if (*isw == ISW_COMMIT) {
    int i;
    for (i=0; i<3; i++)
      thisObj->cState[i] = thisObj->tState[i];

  } else if (*isw == ISW_REVERT) {
    int i;
    for (i=0; i<3; i++)
      thisObj->tState[i] = thisObj->cState[i];

  } else if (*isw == ISW_REVERT_TO_START) {
    int i;
    for (i=0; i<3; i++) {
      thisObj->cState[i] = 0.0;
      thisObj->tState[i] = 0.0;
    }
    
  } else if (*isw == ISW_FORM_TANG_AND_RESID) {

    int i;
    for (i=0; i<3; i++) 
      thisObj->tState[i] = strain[i];

    // form tangent
    double E = thisObj->theParam[0];
    double v = thisObj->theParam[1];
    double fact = E/(1-v*v);
    tang[0] = fact;           // tang[1,1]
    tang[1] = fact*v;         // tang[2,1]
    tang[2] = 0.0;            // tang[3,1]
    tang[3] = fact*v;         // tang[1,2]
    tang[4] = fact;           // tang[2,2]
    tang[5] = 0.0;            // tang[3,2]
    tang[6] = 0.0;            // tang[1,3]
    tang[7] = 0.0;            // tang[2,3]
    tang[8] = fact*(1-v)/2.0; // tang[3,3]

    stress[0] = tang[0]*strain[0] + tang[3]*strain[1];
    stress[1] = tang[1]*strain[0] + tang[4]*strain[1];
    stress[2] = 2.0 * tang[8] * strain[2];
  }

}


