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
// $Date: 2003-06-10 00:36:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/LobattoBeamIntegration.cpp,v $

#include <LobattoBeamIntegration.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_LobattoBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 3) {
	opserr<<"insufficient arguments:integrationTag,secTag,N\n";
	return 0;
    }

    // inputs: integrationTag,secTag,N
    int iData[3];
    int numData = 3;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) return 0;

    integrationTag = iData[0];
    if(iData[2] > 0) {
	secTags.resize(iData[2]);
    } else {
	secTags = ID();
    }
    for(int i=0; i<secTags.Size(); i++) {
	secTags(i) = iData[1];
    }

    return new LobattoBeamIntegration;
}

LobattoBeamIntegration::LobattoBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_Lobatto)
{
  // Nothing to do
}

LobattoBeamIntegration::~LobattoBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
LobattoBeamIntegration::getCopy(void)
{
  return new LobattoBeamIntegration();
}

void
LobattoBeamIntegration::getSectionLocations(int numSections, 
					    double L,
					    double *xi)
{
  switch(numSections) {
    
  case 2:
    xi[0] = -1.0;
    xi[1] =  1.0;
    break;
    
  case 3:
    xi[0] = -1.0;
    xi[1] =  0.0;
    xi[2] =  1.0;
    break;
    
  case 4:
    xi[0] = -1.0;
    xi[1] = -0.44721360;
    xi[2] =  0.44721360;
    xi[3] =  1.0;
    break;
    
  case 5:
    xi[0] = -1.0;
    xi[1] = -0.65465367;
    xi[2] =  0.0;
    xi[3] =  0.65465367;
    xi[4] =  1.0;
    break;
    
  case 6:
    xi[0] = -1.0;
    xi[1] = -0.7650553239;
    xi[2] = -0.2852315164;
    xi[3] =  0.2852315164;
    xi[4] =  0.7650553239;
    xi[5] =  1.0;
    break;
    
  case 7:
    xi[0] = -1.0;
    xi[1] = -0.8302238962;
    xi[2] = -0.4688487934;
    xi[3] =  0.0;
    xi[4] =  0.4688487934;
    xi[5] =  0.8302238962;
    xi[6] =  1.0;
    break;

  case 8:
    xi[0] = -1.0;
    xi[1] = -0.8717401485;
    xi[2] = -0.5917001814;
    xi[3] = -0.2092992179;
    xi[4] =  0.2092992179;
    xi[5] =  0.5917001814;
    xi[6] =  0.8717401485;
    xi[7] =  1.0;
    break;
    
  case 9:
    xi[0] = -1.0;
    xi[1] = -0.8997579954;
    xi[2] = -0.6771862795;
    xi[3] = -0.3631174638;
    xi[4] =  0.0;
    xi[5] =  0.3631174638;
    xi[6] =  0.6771862795;
    xi[7] =  0.8997579954;
    xi[8] =  1.0;
    break;

  case 10:
    xi[0] = -1.0;
    xi[1] = -0.9195339082;
    xi[2] = -0.7387738651;
    xi[3] = -0.4779249498;
    xi[4] = -0.1652789577;
    xi[5] =  0.1652789577;
    xi[6] =  0.4779249498;
    xi[7] =  0.7387738651;
    xi[8] =  0.9195339082;
    xi[9] =  1.0;
    break;

  default:
    opserr << "LobattoBeamIntegration -- max # integration points is 10\n";
    break;
  
  }
  
  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
LobattoBeamIntegration::getSectionWeights(int numSections, double L,
					  double *wt)
{
  switch (numSections) {
    
  case 2:
    wt[0] = 1.0;
    wt[1] = 1.0;
    break;
    
  case 3:
    wt[0] = 0.333333333333333;
    wt[1] = 1.333333333333333;
    wt[2] = 0.333333333333333;
    break;
    
  case 4:    
    wt[0] = 0.166666666666667;
    wt[1] = 0.833333333333333;
    wt[2] = 0.833333333333333;
    wt[3] = 0.166666666666667;
    break;
    
  case 5:
    wt[0] = 0.1;
    wt[1] = 0.5444444444;
    wt[2] = 0.7111111111;
    wt[3] = 0.5444444444;
    wt[4] = 0.1;
    break;
    
  case 6:    
    wt[0] = 0.06666666667;
    wt[1] = 0.3784749562;
    wt[2] = 0.5548583770;
    wt[3] = 0.5548583770;
    wt[4] = 0.3784749562;
    wt[5] = 0.06666666667;
    break;
    
  case 7:    
    wt[0] = 0.04761904762;
    wt[1] = 0.2768260473;
    wt[2] = 0.4317453812;
    wt[3] = 0.4876190476;
    wt[4] = 0.4317453812;
    wt[5] = 0.2768260473;
    wt[6] = 0.04761904762;
    break;

  case 8:    
    wt[0] = 0.03571428571;
    wt[1] = 0.2107042271;
    wt[2] = 0.3411226924;
    wt[3] = 0.4124587946;
    wt[4] = 0.4124587946;
    wt[5] = 0.3411226924;
    wt[6] = 0.2107042271;
    wt[7] = 0.03571428571;
    break;

  case 9:    
    wt[0] = 0.02777777778;
    wt[1] = 0.1654953615;
    wt[2] = 0.2745387125;
    wt[3] = 0.3464285109;
    wt[4] = 0.3715192743;
    wt[5] = 0.3464285109;
    wt[6] = 0.2745387125;
    wt[7] = 0.1654953615;
    wt[8] = 0.02777777778;
    break;

  case 10:    
    wt[0] = 0.02222222222;
    wt[1] = 0.1333059908;
    wt[2] = 0.2248893421;
    wt[3] = 0.2920426836;
    wt[4] = 0.3275397611;
    wt[5] = 0.3275397611;
    wt[6] = 0.2920426836;
    wt[7] = 0.2248893421;
    wt[8] = 0.1333059908;
    wt[9] = 0.02222222222;
    break;

  default:
    opserr << "LobattoBeamIntegration -- max # integration points is 10\n";
    break;

  }
  
  for (int i = 0; i < numSections; i++)
    wt[i] *= 0.5;
}

void
LobattoBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"Lobatto\"}";
	}
    
	else {
		s << "Lobatto" << endln;
	}
}
