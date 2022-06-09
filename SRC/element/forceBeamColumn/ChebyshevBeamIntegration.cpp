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
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ChebyshevBeamIntegration.cpp,v $

#include <ChebyshevBeamIntegration.h>
#include <elementAPI.h>
#include <ID.h>
#include <math.h>

void* OPS_ChebyshevBeamIntegration(int& integrationTag, ID& secTags)
{
  int nArgs = OPS_GetNumRemainingInputArgs();

  if (nArgs < 3) {
    opserr<<"insufficient arguments:integrationTag,secTag,N -or- N,*secTagList\n";
    return 0;
  }
  
  // Read tag
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData,&iData[0]) < 0) {
    opserr << "ChebyshevBeamIntegration - unable to read int data" << endln;
    return 0;
  }
  integrationTag = iData[0];
  
  if (nArgs == 3) {
    // inputs: integrationTag,secTag,N
    numData = 1;
    int Nsections;
    if (OPS_GetIntInput(&numData,&Nsections) < 0) {
      opserr << "ChebyshevBeamIntegration - Unable to read number of sections" << endln;
      return 0;
    }
    if (Nsections < 0)
      return 0;
    
    if (Nsections > 0) {
      secTags.resize(Nsections);
    } else {
      secTags = ID();
    }
    for (int i=0; i<secTags.Size(); i++) {
      secTags(i) = iData[1];
    }
  }
  else {
    // inputs: integrationTag,N,*secTagList
    int Nsections = iData[1];
    if (Nsections < 0)
      return 0;
    int *sections = new int[Nsections];
    if (OPS_GetIntInput(&Nsections,sections) < 0) {
      opserr << "ChebyshevBeamIntegration - Unable to read section tags" << endln;
      return 0;
    }
    if (Nsections > 0) {
      secTags.resize(Nsections);
    } else {
      secTags = ID();
    }
    for (int i=0; i<secTags.Size(); i++) {
      secTags(i) = sections[i];
    }      
    delete [] sections;
  }
  
  return new ChebyshevBeamIntegration();
}

ChebyshevBeamIntegration::ChebyshevBeamIntegration(int t):
  BeamIntegration(BEAM_INTEGRATION_TAG_Chebyshev), type(t)
{
  if (type < 0 || type > 2)
    type = 0;
}

ChebyshevBeamIntegration::~ChebyshevBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
ChebyshevBeamIntegration::getCopy(void)
{
  return new ChebyshevBeamIntegration(type);
}

void
ChebyshevBeamIntegration::getSectionLocations(int numSections, 
					     double L,
					     double *xi)
{
  double pi = acos(-1.0);

  if (type == 0) {
    for (int i = 0; i < numSections; i++)
      xi[i] = 0.0;
    switch (numSections) {
    case 2:
      xi[0] = -0.5773502692;
      xi[1] =  0.5773502692;
      break;
    case 3:
      xi[0] = -0.7071067812;
      xi[1] =  0.0;
      xi[2] =  0.7071067812;
      break;
    case 4:
      xi[0] = -0.7946544723;
      xi[1] = -0.1875924741;
      xi[2] =  0.1875924741;      
      xi[3] =  0.7946544723;
      break;
    case 5:
      xi[0] = -0.8324974870;
      xi[1] = -0.3745414096;
      xi[2] =  0.0;
      xi[3] =  0.3745414096;      
      xi[4] =  0.8324974870;
      break;
    case 6:
      xi[0] = -0.8662468181;
      xi[1] = -0.4225186538;
      xi[2] = -0.2666354015;
      xi[3] =  0.2666354015;
      xi[4] =  0.4225186538;
      xi[5] =  0.8662468181;      
      break;
    case 7:
      xi[0] = -0.8838617008;
      xi[1] = -0.5296567753;
      xi[2] = -0.3239118105;
      xi[3] =  0.0;
      xi[4] =  0.3239118105;      
      xi[5] =  0.5296567753;      
      xi[6] =  0.8838617008;
      break;
    case 8:
      // complex values
      break;
    case 9:
      xi[0] = -0.9115893077;
      xi[1] = -0.6010186554;
      xi[2] = -0.5287617831;
      xi[3] = -0.1679061842;
      xi[4] =  0.0;
      xi[5] =  0.1679061842;      
      xi[6] =  0.5287617831;      
      xi[7] =  0.6010186554;      
      xi[8] =  0.9115893077;
      break;
    case 10:
      // complex values
      break;	
    default:
      break;
    }
  }
  
  if (type == 1) {
    for (int i = 0; i < numSections; i++)
      xi[i] = cos((2*i+1)*pi/(2*numSections));
  }

  if (type == 2) {
    for (int i = 0; i < numSections; i++)
      xi[i] = cos((i+1)*pi/(numSections+1));      
  }
  
  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
ChebyshevBeamIntegration::getSectionWeights(int numSections, double L,
					   double *wt)
{
  double pi = acos(-1.0);

  if (type == 0) {
    double w = 2.0/numSections;
    for (int i = 0; i < numSections; i++)
      wt[i] = w;    
  }
  if (type == 1) {
    double w = pi/numSections;
    //w = 2.0/numSections;
    for (int i = 0; i < numSections; i++)
      wt[i] = w;
  }

  if (type == 2) {
    for (int i = 0; i < numSections; i++) {
      double s = sin((i+1)*pi/(numSections+1));
      wt[i] = pi/(numSections+1)*s*s;
      //wt[i] *= 4/pi;
    }
  }

  for (int i = 0; i < numSections; i++) {
    wt[i] *= 0.5;
  }
}

void
ChebyshevBeamIntegration::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "{\"type\": \"Chebyshev\"}";
  }
  
  else {
    s << "Chebyshev, type=" << type << endln;
  }
}
