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
  int nArgs = OPS_GetNumRemainingInputArgs();

  if (nArgs < 3) {
    opserr<<"insufficient arguments:integrationTag,secTag,N -or- N,*secTagList\n";
    return 0;
  }
  
  // Read tag
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData,&iData[0]) < 0) {
    opserr << "LobattoBeamIntegration - unable to read int data" << endln;
    return 0;
  }
  integrationTag = iData[0];
  
  if (nArgs == 3) {
    // inputs: integrationTag,secTag,N
    numData = 1;
    int Nsections;
    if (OPS_GetIntInput(&numData,&Nsections) < 0) {
      opserr << "LobattoBeamIntegration - Unable to read number of sections" << endln;
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
      opserr << "LobattoBeamIntegration - Unable to read section tags" << endln;
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
    xi[0] = -1.00000000000000e+00;
    xi[1] = 1.00000000000000e+00;
    break;
    
  case 3:
    xi[0] = -1.00000000000000e+00;
    xi[1] = 0.0;
    xi[2] = 1.00000000000000e+00;
    break;
    
  case 4:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -4.47213595499958e-01;
    xi[2] = 4.47213595499958e-01;
    xi[3] = 1.00000000000000e+00;
    break;
    
  case 5:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -6.54653670707977e-01;
    xi[2] = 0.0;
    xi[3] = 6.54653670707977e-01;
    xi[4] = 1.00000000000000e+00;
    break;
    
  case 6:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -7.65055323929465e-01;
    xi[2] = -2.85231516480645e-01;
    xi[3] = 2.85231516480645e-01;
    xi[4] = 7.65055323929465e-01;
    xi[5] = 1.00000000000000e+00;
    break;
    
  case 7:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.30223896278567e-01;
    xi[2] = -4.68848793470714e-01;
    xi[3] = 0.0;
    xi[4] = 4.68848793470714e-01;
    xi[5] = 8.30223896278567e-01;
    xi[6] = 1.00000000000000e+00;
    break;
    
  case 8:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.71740148509606e-01;
    xi[2] = -5.91700181433143e-01;
    xi[3] = -2.09299217902479e-01;
    xi[4] = 2.09299217902479e-01;
    xi[5] = 5.91700181433143e-01;
    xi[6] = 8.71740148509607e-01;
    xi[7] = 1.00000000000000e+00;
    break;
    
  case 9:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.99757995411461e-01;
    xi[2] = -6.77186279510738e-01;
    xi[3] = -3.63117463826178e-01;
    xi[4] = 0.0;
    xi[5] = 3.63117463826178e-01;
    xi[6] = 6.77186279510738e-01;
    xi[7] = 8.99757995411461e-01;
    xi[8] = 1.00000000000000e+00;
    break;
    
  case 10:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -9.19533908166459e-01;
    xi[2] = -7.38773865105505e-01;
    xi[3] = -4.77924949810445e-01;
    xi[4] = -1.65278957666387e-01;
    xi[5] = 1.65278957666387e-01;
    xi[6] = 4.77924949810445e-01;
    xi[7] = 7.38773865105505e-01;
    xi[8] = 9.19533908166459e-01;
    xi[9] = 1.00000000000000e+00;
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
    wt[0] = 1.00000000000000e+00;
    wt[1] = 1.00000000000000e+00;
    break; // sum = 2.00000000000000e+00
    
  case 3:
    wt[0] = 3.33333333333333e-01;
    wt[1] = 1.33333333333333e+00;
    wt[2] = 3.33333333333334e-01;
    break; // sum = 2.00000000000000e+00
    
  case 4:
    wt[0] = 1.66666666666667e-01;
    wt[1] = 8.33333333333333e-01;
    wt[2] = 8.33333333333334e-01;
    wt[3] = 1.66666666666666e-01;
    break; // sum = 2.00000000000000e+00
    
  case 5:
    wt[0] = 1.00000000000000e-01;
    wt[1] = 5.44444444444445e-01;
    wt[2] = 7.11111111111111e-01;
    wt[3] = 5.44444444444444e-01;
    wt[4] = 1.00000000000000e-01;
    break; // sum = 2.00000000000000e+00
    
  case 6:
    wt[0] = 6.66666666666666e-02;
    wt[1] = 3.78474956297848e-01;
    wt[2] = 5.54858377035486e-01;
    wt[3] = 5.54858377035487e-01;
    wt[4] = 3.78474956297848e-01;
    wt[5] = 6.66666666666668e-02;
    break; // sum = 2.00000000000000e+00
    
  case 7:
    wt[0] = 4.76190476190473e-02;
    wt[1] = 2.76826047361566e-01;
    wt[2] = 4.31745381209863e-01;
    wt[3] = 4.87619047619048e-01;
    wt[4] = 4.31745381209863e-01;
    wt[5] = 2.76826047361567e-01;
    wt[6] = 4.76190476190480e-02;
    break; // sum = 2.00000000000000e+00
    
  case 8:
    wt[0] = 3.57142857142856e-02;
    wt[1] = 2.10704227143506e-01;
    wt[2] = 3.41122692483505e-01;
    wt[3] = 4.12458794658704e-01;
    wt[4] = 4.12458794658703e-01;
    wt[5] = 3.41122692483505e-01;
    wt[6] = 2.10704227143506e-01;
    wt[7] = 3.57142857142862e-02;
    break; // sum = 2.00000000000000e+00
    
  case 9:
    wt[0] = 2.77777777777778e-02;
    wt[1] = 1.65495361560806e-01;
    wt[2] = 2.74538712500162e-01;
    wt[3] = 3.46428510973046e-01;
    wt[4] = 3.71519274376417e-01;
    wt[5] = 3.46428510973047e-01;
    wt[6] = 2.74538712500161e-01;
    wt[7] = 1.65495361560806e-01;
    wt[8] = 2.77777777777781e-02;
    break; // sum = 2.00000000000000e+00
    
  case 10:
    wt[0] = 2.22222222222225e-02;
    wt[1] = 1.33305990851070e-01;
    wt[2] = 2.24889342063126e-01;
    wt[3] = 2.92042683679684e-01;
    wt[4] = 3.27539761183898e-01;
    wt[5] = 3.27539761183897e-01;
    wt[6] = 2.92042683679685e-01;
    wt[7] = 2.24889342063126e-01;
    wt[8] = 1.33305990851070e-01;
    wt[9] = 2.22222222222225e-02;
    break; // sum = 2.00000000000000e+00
    
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
