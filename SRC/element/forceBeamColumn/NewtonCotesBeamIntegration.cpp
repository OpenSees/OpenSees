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
// $Date: 2006-01-17 21:12:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/NewtonCotesBeamIntegration.cpp,v $

#include <NewtonCotesBeamIntegration.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_NewtonCotesBeamIntegration(int& integrationTag, ID& secTags)
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

    return new NewtonCotesBeamIntegration;
}

NewtonCotesBeamIntegration::NewtonCotesBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_NewtonCotes)
{
  // Nothing to do
}

NewtonCotesBeamIntegration::~NewtonCotesBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
NewtonCotesBeamIntegration::getCopy(void)
{
  return new NewtonCotesBeamIntegration();
}

void
NewtonCotesBeamIntegration::getSectionLocations(int numSections, double L,
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
    xi[1] = -0.3333333333;
    xi[2] =  0.3333333333;
    xi[3] =  1.0;
    break;
    
  case 5:
    xi[0] = -1.0;
    xi[1] = -0.5;
    xi[2] =  0.0;
    xi[3] =  0.5;
    xi[4] =  1.0;
    break;
    
  case 6:
    xi[0] = -1.0;
    xi[1] = -0.6;
    xi[2] = -0.2;
    xi[3] =  0.2;
    xi[4] =  0.6;
    xi[5] =  1.0;
    break;
    
  case 7:
    xi[0] = -1.0;
    xi[1] = -0.6666666667;
    xi[2] = -0.3333333333;
    xi[3] =  0.0;
    xi[4] =  0.3333333333;
    xi[5] =  0.6666666667;
    xi[6] =  1.0;
    break;

  case 8:
    xi[0] = -1.0;
    xi[1] = -0.7142857143;
    xi[2] = -0.4285714286;
    xi[3] = -0.1428571429;
    xi[4] =  0.1428571429;
    xi[5] =  0.4285714286;
    xi[6] =  0.7142857143;
    xi[7] =  1.0;
    break;
    
  case 9:
    xi[0] = -1.0;
    xi[1] = -0.75;
    xi[2] = -0.5;
    xi[3] = -0.25;
    xi[4] =  0.0;
    xi[5] =  0.25;
    xi[6] =  0.5;
    xi[7] =  0.75;
    xi[8] =  1.0;
    break;

  case 10:
    xi[0] = -1.0;
    xi[1] = -0.7777777778;
    xi[2] = -0.5555555556;
    xi[3] = -0.3333333333;
    xi[4] = -0.1111111111;
    xi[5] =  0.1111111111;
    xi[6] =  0.3333333333;
    xi[7] =  0.5555555556;
    xi[8] =  0.77777777778;
    xi[9] =  1.0;
    break;
  }
  
  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
NewtonCotesBeamIntegration::getSectionWeights(int numSections, double L,
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
    wt[0] = 0.25;
    wt[1] = 0.75;
    wt[2] = 0.75;
    wt[3] = 0.25;
    break;
    
  case 5:
    wt[0] = 0.1555555556;
    wt[1] = 0.7111111111;
    wt[2] = 0.2666666667;
    wt[3] = 0.7111111111;
    wt[4] = 0.1555555556;
    break;
    
  case 6:    
    wt[0] = 0.1319444444;
    wt[1] = 0.5208333333;
    wt[2] = 0.3472222222;
    wt[3] = 0.3472222222;
    wt[4] = 0.5208333333;
    wt[5] = 0.1319444444;
    break;
    
  case 7:    
    wt[0] = 0.09761904762;
    wt[1] = 0.5142857143;
    wt[2] = 0.06428571429;
    wt[3] = 0.6476190476;
    wt[4] = 0.06428571429;
    wt[5] = 0.5142857143;
    wt[6] = 0.09761904762;
    break;

  case 8:    
    wt[0] = 0.0869212963;
    wt[1] = 0.4140046296;
    wt[2] = 0.153125;
    wt[3] = 0.3459490741;
    wt[4] = 0.3459490741;
    wt[5] = 0.153125;
    wt[6] = 0.4140046296;
    wt[7] = 0.0869212963;
    break;

  case 9:    
    wt[0] =  0.0697707231;
    wt[1] =  0.4153791887;
    wt[2] = -0.06546737213;
    wt[3] =  0.7404585538;
    wt[4] = -0.3202821869;
    wt[5] =  0.7404585538;
    wt[6] = -0.06546737213;
    wt[7] =  0.4153791887;
    wt[8] =  0.0697707231;
    break;

  case 10:    
    wt[0] = 0.06377232143;
    wt[1] = 0.3513616071;
    wt[2] = 0.02410714286;
    wt[3] = 0.4317857143;
    wt[4] = 0.1289732143;
    wt[5] = 0.1289732143;
    wt[6] = 0.4317857143;
    wt[7] = 0.02410714286;
    wt[8] = 0.3513616071;
    wt[9] = 0.06377232143;
    break;
  }
  
  for (int i = 0; i < numSections; i++)
    wt[i] *= 0.5;
}

void
NewtonCotesBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"NewtonCotes\"}";
	}
	
	else {
		s << "NewtonCotes" << endln;
	}
}
