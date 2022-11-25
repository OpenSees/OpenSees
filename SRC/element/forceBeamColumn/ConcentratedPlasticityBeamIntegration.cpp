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
**                                                                    **
** ****************************************************************** *
** This Beam Integration Module was Developed by:                     **
**   Silvia Mazzoni (silviamazzoni@yahoo.com)                         **
**   Michael H. Scott (michael@portwooddigital.com)                   **
**   pushed to OpenSees November 2022                                 **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2022-11-20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ConcentratedPlasticityBeamIntegration.cpp,v $

#include <ConcentratedPlasticityBeamIntegration.h>

#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <elementAPI.h>

void* OPS_ConcentratedPlasticityBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments:integrationTag,secTagI,secTagJ,secTagE\n";
	return 0;
    }

    // inputs: 
    int iData[4];
    int numData = 2;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {   
	opserr << "WARNING: failed to get tag and secTagI\n";
	return 0;
    }
    numData = 1;
    if(OPS_GetIntInput(&numData,&iData[2]) < 0) {
	opserr << "WARNING: failed to get secTagJ\n";
	return 0;
    }
    if(OPS_GetIntInput(&numData,&iData[3]) < 0) {
	opserr << "WARNING: failed to get secTagE\n";
	return 0;
    }

    
    int N = 5;
    integrationTag = iData[0];
    secTags.resize(N);
    secTags(0) = iData[1];
    secTags(1) = iData[3];
    secTags(2) = iData[3];
    secTags(3) = iData[3];
    secTags(4) = iData[2];


    
    
    
    return new ConcentratedPlasticityBeamIntegration();
}

ConcentratedPlasticityBeamIntegration::ConcentratedPlasticityBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_ConcentratedPlasticity)
{
  // Nothing to do


}

ConcentratedPlasticityBeamIntegration::~ConcentratedPlasticityBeamIntegration()
{
  // Nothing to do
}

void
ConcentratedPlasticityBeamIntegration::getSectionLocations(int numSections, double L,
					       double *xi)
{
	
    // section locations
    xi[0] = 0.0; // node I
    xi[1] = 0.0; // elastic
    xi[2] = 0.5; // elastic
    xi[3] = 1.0; // elastic
    xi[4] = 1.0; // node J
}

void
ConcentratedPlasticityBeamIntegration::getSectionWeights(int numSections, double L,
					     double *wt)
{
	double oneOverL = 1.0/L;
	int N = 5 ;
    
    // start with the two end nodes at the first two indices, then put them back
    Vector pt(N);
    
    pt[0] = 0.0; // node I
    pt[1] = 1.0; // node J    
    pt[2] = 0.0; // elastic
    pt[3] = 0.5; // elastic
    pt[4] = 1.0; // elastic
    
    int nc = 2;
    Vector wc(nc);
    wc[0] = 1.0*oneOverL; // node I
    wc[1] = 1.0*oneOverL; // node J

    // fill in the rest using low-order integration:
	int nf = N-nc;
    Vector R(nf);
    for (int i = 0; i < nf; i++) {
    	double sum = 0.0;
    	for (int j = 0; j < nc; j++)
		sum += pow(pt(j),i)*wc(j);
	  	R(i) = 1.0/(i+1) - sum;
    }
    
    Matrix J(nf,nf);
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < nf; j++)
	J(i,j) = pow(pt(nc+j),i);
    
    Vector wf(nf);
    
    J.Solve(R, wf);
    
     
    wt[0] = wc(0); // node I
    wt[1] = wf(0); // elastic
    wt[2] = wf(1); // elastic
    wt[3] = wf(2); // elastic
    wt[4] = wc(1); // node J
    
}

BeamIntegration*
ConcentratedPlasticityBeamIntegration::getCopy(void)
{

  return new ConcentratedPlasticityBeamIntegration();
}

int
ConcentratedPlasticityBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = 0.;
  data(1) = 0.;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "ConcentratedPlasticityBeamIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}


int
ConcentratedPlasticityBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				  FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "ConcentratedPlasticityBeamIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  

  return 0;
}

int
ConcentratedPlasticityBeamIntegration::setParameter(const char **argv, int argc,
				      Parameter &param)
{
	
	// there really aren't any parameters that you can change!

	return 0;
}

int
ConcentratedPlasticityBeamIntegration::updateParameter(int parameterID,
					 Information &info)
{

	return 0;
}

int
ConcentratedPlasticityBeamIntegration::activateParameter(int paramID)
{
  

  return 0;
}

void
ConcentratedPlasticityBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"ConcentratedPlasticity\", ";
	}
	
	else {
		s << "ConcentratedPlasticity" << endln;
	}
	return;
}

void 
ConcentratedPlasticityBeamIntegration::getLocationsDeriv(int numSections, double L,
					   double dLdh, double *dptsdh)
{

    dptsdh[0] =  0.0;
    dptsdh[1] =  0.0;
    dptsdh[2] =  0.0;
    dptsdh[3] =  0.0;
    dptsdh[4] =  0.0;
    
  return;

  if (dLdh != 0.0) {
    // STILL TO DO
    //opserr << "getPointsDeriv -- to do" << endln;
  }
      
    return;
}

void
ConcentratedPlasticityBeamIntegration::getWeightsDeriv(int numSections, double L,
					 double dLdh, double *dwtsdh)
{


    dwtsdh[0] = 0.;
    dwtsdh[1] = 0.;
    dwtsdh[2] = 0.;
    dwtsdh[3] = 0.;
    dwtsdh[4] = 0.;

    
  return;
}
