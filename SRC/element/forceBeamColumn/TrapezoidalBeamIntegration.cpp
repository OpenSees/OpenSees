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
// $Date: 2007-10-12 20:57:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/TrapezoidalBeamIntegration.cpp,v $

#include <TrapezoidalBeamIntegration.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_TrapezoidalBeamIntegration(int& integrationTag, ID& secTags)
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
    opserr << "TrapezoidalBeamIntegration - unable to read int data" << endln;
    return 0;
  }
  integrationTag = iData[0];
  
  if (nArgs == 3) {
    // inputs: integrationTag,secTag,N
    numData = 1;
    int Nsections;
    if (OPS_GetIntInput(&numData,&Nsections) < 0) {
      opserr << "TrapezoidalBeamIntegration - Unable to read number of sections" << endln;
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
      opserr << "TrapezoidalBeamIntegration - Unable to read section tags" << endln;
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
  
  return new TrapezoidalBeamIntegration;
}

TrapezoidalBeamIntegration::TrapezoidalBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_Trapezoidal)
{
  // Nothing to do
}

TrapezoidalBeamIntegration::~TrapezoidalBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
TrapezoidalBeamIntegration::getCopy(void)
{
  return new TrapezoidalBeamIntegration();
}

void
TrapezoidalBeamIntegration::getSectionLocations(int numSections, double L,
						double *xi)
{
  if (numSections > 1) {
    xi[0] = -1.0;
    xi[numSections-1] = 1.0;

    double dxi = 2.0/(numSections-1);

    for (int i = 1; i < numSections-1; i++)
      xi[i] = -1.0 + dxi*i;
  }

  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
TrapezoidalBeamIntegration::getSectionWeights(int numSections, double L,
					      double *wt)
{
  if (numSections > 1) {

    double wti = 2.0/(numSections-1);

    for (int i = 1; i < numSections-1; i++)
      wt[i] = wti;

    wt[0] = wt[numSections-1] = 0.5*wti;
  }

  for (int i = 0; i < numSections; i++)
    wt[i] *= 0.5;
}

void
TrapezoidalBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"Trapezoidal\"}";
	}
	
	else {
		s << "Trapezoidal" << endln;
	}
}
