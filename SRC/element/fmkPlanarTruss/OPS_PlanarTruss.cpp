/********************************************************************************
(C) Copyright 2001-2022, The Regents of the University of California    
All Rights Reserved.                                               

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************* */

#include "PlanarTruss.h"
#include <elementAPI.h>

#include <UniaxialMaterial.h>

void *OPS_PlanarTruss(void) {

  //
  // create 2 arrays; one each for integer and double command line arguments
  //
  
  int iData[4];    // integer args on command line of size 1, change to suit
  double dData[1]; // double args on command line of size 1, change to suit

  //
  // read values from command line of script
  //
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 5) {
    opserr << "WARNING incorrect # args, want: element PlanarTruss $id $iNode $jNode $matTag $A\n";
    return 0;
  }
  
  int numData = 4; // num integer args
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING failed to read integers, command element PlanarTruss\n";
    return 0;
  }
  
  numData = 1; // reset to num double args
  if (OPS_GetDoubleInput(&numData, &dData[0]) < 0) {
    opserr << "WARNING failed to read doubles, command element PlanarTruss\n";
    return 0;
  }

  //
  // return pointer to new element
  //

  int matTag = iData[3];
  UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matTag);
  
  return new PlanarTruss(iData[0], iData[1], iData[2], theMat, dData[0]);
}
