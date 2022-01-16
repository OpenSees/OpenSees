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
                                                                        
// $Revision: 1.8 $
// $Date: 2010-02-16 18:57:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PeerMotion.cpp,v $                                                                        
// Written: fmk 
// Created: 10/06
//
// Purpose: This file contains the class definition for PeerMotion.
// PeerMotion is a concrete class. A PeerMotion object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor. 
//
// What: "@(#) PeerMotion.C, revA"


#include <PeerMotion.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

#include <stdio.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

//#include <Socket.h>

#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGet(char const *URL, char const *page, unsigned int port, char **dataPtr);



#include <elementAPI.h>
#define OPS_Export 

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_PeerMotion)
{
  // Pointer to a uniaxial material that will be returned
  TimeSeries *theSeries = 0;
  
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  
  if (numRemainingArgs < 4) {
    opserr << "WARNING: invalid num args PeerMotion <tag?> $eqMotion $station $type $factor\n";
    return 0;
  }

  int tag = 0;     // default tag = 0
  double factor = 0.0; 
  int numData = 0;

  char *eqMotion = 0;
  char *station = 0;
  char *type = 0;

  // get tag if provided
  if (numRemainingArgs == 5 || numRemainingArgs == 7 || numRemainingArgs == 9) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
      opserr << "WARNING invalid series tag in Constant tag?" << endln;
      return 0;
    }
    numRemainingArgs -= 1;
  }
  
  if ((OPS_GetStringCopy(&eqMotion) != 0) || eqMotion == 0) {
    opserr << "WARNING invalid eqMotion for PeerMotion with tag: " << tag << endln;
    return 0;
  }    

    if ((OPS_GetStringCopy(&station) != 0) || station == 0) {
    opserr << "WARNING invalid station for PeerMotion with tag: " << tag << endln;
    return 0;
  }    

    if ((OPS_GetStringCopy(&type) != 0) || type == 0) {
    opserr << "WARNING invalid type  for PeerMotion with tag: " << tag << endln;
    return 0;
  }    


  if (OPS_GetDouble(&numData, &factor) != 0) {
    opserr << "WARNING invalid facor in PeerMotion Series with tag?" << tag << endln;
    return 0;
  }
  
  theSeries = new PeerMotion(tag, eqMotion, station, type, factor);

  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating PeerMotion with tag: " << tag << "\n";
    return 0;
  }

  delete [] eqMotion;
  delete [] station;
  delete [] type;

  return theSeries;
}


PeerMotion::PeerMotion()	
  :TimeSeries(TSERIES_TAG_PeerMotion),
   thePath(0), dT(0.0), 
   cFactor(0.0), dbTag1(0), dbTag2(0), lastSendCommitTag(-1)
{
  // does nothing
}

		   
PeerMotion::PeerMotion(int tag,
		       const char *earthquake,
		       const char *station,
		       const char *type,
		       double theFactor)
  :TimeSeries(tag, TSERIES_TAG_PeerMotion),
   thePath(0), dT(0.0), 
   cFactor(theFactor), dbTag1(0), dbTag2(0), lastSendCommitTag(-1), lastChannel(0)
{
  char peerPage[124];
  char *nextData, *eqData;
  int nPts,i;
  char tmp1[100];

  if (earthquake != 0 && station != 0 && type != 0) {
    
    if ((strcmp(type,"ACCEL") == 0) || (strcmp(type,"-accel") == 0) || (strcmp(type,"-ACCEL") == 0)
	|| (strcmp(type,"accel") == 0) || (strcmp(type,"ATH") == 0) || (strcmp(type,"-ATH") == 0)) {
      sprintf(peerPage, "/smcat/data/ath/%s/%s.AT2",earthquake,station);
    } else if ((strcmp(type,"DISP") == 0) || (strcmp(type,"-disp") == 0) || (strcmp(type,"-DISP") == 0)
	       || (strcmp(type,"adisp") == 0) || (strcmp(type,"DTH") == 0) || (strcmp(type,"-DTH") == 0)) {
      sprintf(peerPage, "/smcat/data/dth/%s/%s.DT2",earthquake,station);
    } else {
      opserr << "PeerMotion::PeerMotion() - not a valid type:" << type << " (-DISP or -ACCEL requiured)\n";
      return;
    }

    if (httpGet("peer.berkeley.edu",peerPage,80,&eqData) != 0) {
      opserr << "PeerMotion::PeerMotion() - could not connect to PEER Database, ";
      return; 
    }

    if (eqData == 0) {
      opserr << "PeerMotion::PeerMotion() - NO data returned ";
      return; 
    }

    nextData = strstr(eqData,"Page Not Found");
    if (nextData != 0) {
      opserr << "PeerMotion::PeerMotion() - could not get Data for record from Database, ";
      opserr << "page: " << peerPage << " missing \n";
      free(eqData);
      return;
    }

    nextData = strstr(eqData,"NPTS");
    if (nextData == NULL) {
      opserr << "PeerMotion::PeerMotion() - could not find nPts in record, send email opensees-support@berkeley.edu";
      free(eqData);
      return;
    }

    nextData+=5; // NPTS=
    nPts = atoi(nextData);

    nextData = strstr(eqData, "DT");
    if (nextData == NULL) {
      nextData = strstr(eqData, "dt");
      if (nextData == NULL) {
	opserr << "PeerMotion::PeerMotion() - could not find dt in record, send email opensees-support@berkeley.edu";
	free(eqData);
	return;
      }
    }

    nextData+=4; //DT= dT UNIT
    dT = strtod(nextData, &nextData);
    
    sscanf(nextData, "%s", tmp1);
    nextData += strlen(tmp1)+1;
    
    sscanf(nextData, "%s", tmp1);

    thePath = new Vector(nPts);
    //  data = (double *)malloc(nPts*sizeof(double));

    for (i=0; i<nPts; i++) {
      double value = strtod(nextData, &nextData);
      (*thePath)(i) = value;
    }
    
    free(eqData);
  }
}


PeerMotion::PeerMotion(int tag,
		       Vector *theDataPoints,
		       double theTimeStep, 
		       double theFactor)
  :TimeSeries(tag, TSERIES_TAG_PeerMotion),
   thePath(0), dT(theTimeStep), 
   cFactor(theFactor), dbTag1(0), dbTag2(0), lastSendCommitTag(-1), lastChannel(0)
{
  if (theDataPoints != 0)
    thePath = new Vector(*theDataPoints);
}

TimeSeries *
PeerMotion::getCopy(void) 
{
  return new PeerMotion(this->getTag(), thePath, dT, cFactor);
}


PeerMotion::~PeerMotion()
{
  if (thePath != 0)
    delete thePath;
}

double
PeerMotion::getTimeIncr (double pseudoTime)
{
  // NEED TO FILL IN, FOR NOW return 1.0
  return 1.0;
}

double
PeerMotion::getFactor(double pseudoTime)
{
  // check for a quick return
  if (pseudoTime < 0.0 || thePath == 0)
    return 0.0;

  // determine indexes into the data array whose boundary holds the time
  double incr = pseudoTime/dT; 
  int incr1 = floor(incr);
  int incr2 = incr1+1;

  if (incr2 >= thePath->Size())
    return 0.0;

  double value1 = (*thePath)[incr1];
  double value2 = (*thePath)[incr2];
  return cFactor*(value1 + (value2-value1)*(pseudoTime/dT - incr1));
}

double
PeerMotion::getDuration()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PeerMotion::getDuration() on empty Vector" << endln;
	return 0.0;
  }
  return (thePath->Size() * dT);
}

double
PeerMotion::getPeakFactor()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PeerMotion::getPeakFactor() on empty Vector" << endln;
    return 0.0;
  }

  double peak = fabs((*thePath)[0]);
  int num = thePath->Size();
  double temp;
  
  for (int i = 1; i < num; i++)
  {
    temp = fabs((*thePath)[i]);
    if (temp > peak)
      peak = temp;
  }
  
  return (peak*cFactor);
}


double 
PeerMotion::getDt()
{
  return dT;
}
int
PeerMotion::getNPts()
{
  if (thePath == 0)
    return 0;
  else
    return thePath->Size();
}


int
PeerMotion::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  Vector data(5);
  data(0) = cFactor;
  data(1) = dT;
  data(2) = -1;
  
  if (thePath != 0) {
    int size = thePath->Size();
    data(2) = size;
    if (otherDbTag == 0)
      otherDbTag = theChannel.getDbTag();
    data(3) = otherDbTag;
  }

  if ((lastSendCommitTag == -1) && (theChannel.isDatastore() == 1)) {
    lastSendCommitTag = commitTag;
  }

  data(4) = lastSendCommitTag;

  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "PeerMotion::sendSelf() - channel failed to send data\n";
    return result;
  }

  // we only send the vector data if this is the first time it is sent to the database
  // or the channel is for sending the data to a remote process

  if ((lastSendCommitTag == commitTag) || (theChannel.isDatastore() == 0)) {
    if (thePath != 0) {
      result = theChannel.sendVector(otherDbTag, commitTag, *thePath);
      if (result < 0) {
	opserr << "PeerMotion::sendSelf() - ";
	opserr << "channel failed to send the Path Vector\n";
	return result;  
      }
    }
  }

  return 0;
}


int 
PeerMotion::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  Vector data(5);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "PeerMotion::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }

  cFactor = data(0);
  dT = data(1);
  int size = data(2);
  otherDbTag = data(3);
  lastSendCommitTag = data(4);
  
  // get the path vector, only receive it once as it can't change
  if (thePath == 0 && size > 0) {
    thePath = new Vector(size);
    if (thePath == 0 || thePath->Size() == 0) {
      opserr << "PeerMotion::recvSelf() - ran out of memory";
      opserr << " a Vector of size: " <<  size << endln;  
      if (thePath != 0)
	delete thePath;
      thePath = 0;
      return -1;
    }

    result = theChannel.recvVector(otherDbTag, lastSendCommitTag, *thePath);    
    if (result < 0) {
      opserr << "PeerMotion::recvSelf() - ";
      opserr << "channel failed to receive the Path Vector\n";
      return result;  
    }
  }

  return 0;    
}


void
PeerMotion::Print(OPS_Stream &s, int flag)
{
    s << "Path Time Series: constant factor: " << cFactor;
    s << " dT: " << dT << endln;
    if (flag == 1 && thePath != 0) {
      s << " specified path: " << *thePath;

    }
}
