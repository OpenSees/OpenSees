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
// $Date: 2006-08-22 19:05:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/RadauBeamIntegration.cpp,v $

#include <RadauBeamIntegration.h>
#include <math.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_RadauBeamIntegration(int& integrationTag, ID& secTags)
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
    opserr << "RadauBeamIntegration - unable to read int data" << endln;
    return 0;
  }
  integrationTag = iData[0];
  
  if (nArgs == 3) {
    // inputs: integrationTag,secTag,N
    numData = 1;
    int Nsections;
    if (OPS_GetIntInput(&numData,&Nsections) < 0) {
      opserr << "RadauBeamIntegration - Unable to read number of sections" << endln;
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
      opserr << "RadauBeamIntegration - Unable to read section tags" << endln;
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
  
  return new RadauBeamIntegration;
}

RadauBeamIntegration::RadauBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_Radau)
{
  // Nothing to do
}

RadauBeamIntegration::~RadauBeamIntegration()
{
  // Nothing to do
}

BeamIntegration*
RadauBeamIntegration::getCopy(void)
{
  return new RadauBeamIntegration();
}

#ifdef _WIN32

extern "C" int GAUSSQ(int *kind, int *n, double *alpha, double *beta,
                               int *kpts, double *endpts, double *b,
			       double *t, double *w);

#define gaussq_ GAUSSQ

#else

extern "C" int gaussq_(int *kind, int *n, double *alpha, double *beta,
		       int *kpts, double *endpts, double *b,
		       double *t, double *w);

#endif

void
RadauBeamIntegration::getSectionLocations(int numSections, double L,
					  double *xi)
{
  switch(numSections) {

  case 1:
    xi[0] = -1.0;
    break;

  case 2:
    xi[0] = -1.00000000000000e+00;
    xi[1] = 3.33333333333333e-01;
    break;
    
  case 3:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -2.89897948556636e-01;
    xi[2] = 6.89897948556636e-01;
    break;
    
  case 4:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -5.75318923521694e-01;
    xi[2] = 1.81066271118531e-01;
    xi[3] = 8.22824080974592e-01;
    break;
    
  case 5:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -7.20480271312439e-01;
    xi[2] = -1.67180864737834e-01;
    xi[3] = 4.46313972723752e-01;
    xi[4] = 8.85791607770964e-01;
    break;
    
  case 6:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.02929828402347e-01;
    xi[2] = -3.90928546707272e-01;
    xi[3] = 1.24050379505228e-01;
    xi[4] = 6.03973164252783e-01;
    xi[5] = 9.20380285897063e-01;
    break;
    
  case 7:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.53891342639482e-01;
    xi[2] = -5.38467724060109e-01;
    xi[3] = -1.17343037543100e-01;
    xi[4] = 3.26030619437692e-01;
    xi[5] = 7.03842800663031e-01;
    xi[6] = 9.41367145680430e-01;
    break;
    
  case 8:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -8.87474878926156e-01;
    xi[2] = -6.39518616526215e-01;
    xi[3] = -2.94750565773661e-01;
    xi[4] = 9.43072526611107e-02;
    xi[5] = 4.68420354430821e-01;
    xi[6] = 7.70641893678192e-01;
    xi[7] = 9.55041227122575e-01;
    break;
    
  case 9:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -9.10732089420060e-01;
    xi[2] = -7.11267485915709e-01;
    xi[3] = -4.26350485711139e-01;
    xi[4] = -9.03733696068532e-02;
    xi[5] = 2.56135670833455e-01;
    xi[6] = 5.71383041208739e-01;
    xi[7] = 8.17352784200412e-01;
    xi[8] = 9.64440169705273e-01;
    break;
    
  case 10:
    xi[0] = -1.00000000000000e+00;
    xi[1] = -9.27484374233581e-01;
    xi[2] = -7.63842042420003e-01;
    xi[3] = -5.25646030370079e-01;
    xi[4] = -2.36234469390588e-01;
    xi[5] = 7.60591978379783e-02;
    xi[6] = 3.80664840144725e-01;
    xi[7] = 6.47766687674009e-01;
    xi[8] = 8.51225220581608e-01;
    xi[9] = 9.71175180702247e-01;
    break;    

  default:
    break;
  }

  for (int i = 0; i < numSections; i++)
    xi[i]  = 0.5*(xi[i] + 1.0);
}

void
RadauBeamIntegration::getSectionWeights(int numSections, double L,
					double *wt)
{
  switch (numSections) {
    
  case 1:
    wt[0] = 2.0;
    break;

  case 2:
    wt[0] = 5.00000000000000e-01;
    wt[1] = 1.50000000000000e+00;
    break; // sum = 2.00000000000000e+00
    
  case 3:
    wt[0] = 2.22222222222222e-01;
    wt[1] = 1.02497165237684e+00;
    wt[2] = 7.52806125400935e-01;
    break; // sum = 2.00000000000000e+00
    
  case 4:
    wt[0] = 1.25000000000000e-01;
    wt[1] = 6.57688639960120e-01;
    wt[2] = 7.76386937686343e-01;
    wt[3] = 4.40924422353537e-01;
    break; // sum = 2.00000000000000e+00
    
  case 5:
    wt[0] = 8.00000000000000e-02;
    wt[1] = 4.46207802167141e-01;
    wt[2] = 6.23653045951482e-01;
    wt[3] = 5.62712030298925e-01;
    wt[4] = 2.87427121582452e-01;
    break; // sum = 2.00000000000000e+00
    
  case 6:
    wt[0] = 5.55555555555557e-02;
    wt[1] = 3.19640753220511e-01;
    wt[2] = 4.85387188468970e-01;
    wt[3] = 5.20926783189575e-01;
    wt[4] = 4.16901334311908e-01;
    wt[5] = 2.01588385253481e-01;
    break; // sum = 2.00000000000000e+00
    
  case 7:
    wt[0] = 4.08163265306122e-02;
    wt[1] = 2.39227489225312e-01;
    wt[2] = 3.80949873644231e-01;
    wt[3] = 4.47109829014566e-01;
    wt[4] = 4.24703779005955e-01;
    wt[5] = 3.18204231467302e-01;
    wt[6] = 1.48988471112020e-01;
    break; // sum = 2.00000000000000e+00
    
  case 8:
    wt[0] = 3.12500000000000e-02;
    wt[1] = 1.85358154802980e-01;
    wt[2] = 3.04130620646785e-01;
    wt[3] = 3.76517545389118e-01;
    wt[4] = 3.91572167452494e-01;
    wt[5] = 3.47014795634501e-01;
    wt[6] = 2.49647901329865e-01;
    wt[7] = 1.14508814744258e-01;
    break; // sum = 2.00000000000000e+00
    
  case 9:
    wt[0] = 2.46913580246915e-02;
    wt[1] = 1.47654019046315e-01;
    wt[2] = 2.47189378204593e-01;
    wt[3] = 3.16843775670438e-01;
    wt[4] = 3.48273002772966e-01;
    wt[5] = 3.37693966975930e-01;
    wt[6] = 2.86386696357232e-01;
    wt[7] = 2.00553298024552e-01;
    wt[8] = 9.07145049232830e-02;
    break; // sum = 2.00000000000000e+00
    
  case 10:
    wt[0] = 2.00000000000001e-02;
    wt[1] = 1.20296670557482e-01;
    wt[2] = 2.04270131879001e-01;
    wt[3] = 2.68194837841179e-01;
    wt[4] = 3.05859287724423e-01;
    wt[5] = 3.13582457226938e-01;
    wt[6] = 2.90610164832919e-01;
    wt[7] = 2.39193431714380e-01;
    wt[8] = 1.64376012736922e-01;
    wt[9] = 7.36170054867581e-02;
    break; // sum = 2.00000000000000e+00
    
  default:
    break;
  }
  
  for (int i = 0; i < numSections; i++)
    wt[i] *= 0.5;
}

void
RadauBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"Radau\"}";
	}
	
	else {
		s << "Radau" << endln;
	}
}
