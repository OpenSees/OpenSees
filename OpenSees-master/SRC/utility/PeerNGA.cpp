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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-10-09 21:33:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/PeerNGA.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/07
//
// Description: This file contains the class definition for File used in SimulationINformation.
//
// What: "@(#) SimulationInformation.h, revA"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <StringContainer.h>


#ifdef _WIN32
extern int __cdecl
#else
extern int
#endif
httpGet(char const *URL, char const *page, unsigned int port, char **dataPtr) ;


#ifdef _WIN32
int __cdecl
#else
int
#endif
peerGET(const char *page,
	char **res) {
  
  char *URL="peer.berkeley.edu";
  
  int result = httpGet(URL, page, 80, res);

  if (result == -1) 
    result = httpGet(URL, page, 80, res);

  if (result == -1) {
    fprintf(stderr, "ERROR: peerGET - no data\n");
    return -1;

  }

  char *failure = strstr(*res,"400 Bad Request");
  
  if (failure != 0) {
    fprintf(stderr, "ERROR: peerGET\n");
    return -1;
  }
  
  failure = strstr(*res,"401 Unauthorized");
  
  if (failure != 0) {
    fprintf(stderr, "ERROR: peerGET\n");
    return -2;
  }
  
  failure = strstr(*res,"404 Not Found");

  if (failure != 0) {
    fprintf(stderr, "ERROR: peerGET\n");
    return -3;
  }

  failure = strstr(*res,"410 Gone");

  if (failure != 0) {
    fprintf(stderr, "ERROR: peerGET\n");
    return -4;
  }

  return 0;  
}


#ifdef _WIN32
int __cdecl
#else
int
#endif

peerSearchNGA(const char *Eq,
	      const char *SoilType,
	      const char *Fault,
	      const char *MagLo,
	      const char *MagHi,
	      const char *DistLo,
	      const char *DistHi,
	      const char *VsLo,
	      const char *VsHi,
	      const char *PgaLo,
	      const char *PgaHi,
	      const char *LatSW,
	      const char *LatNE,
	      const char *LngSW,
	      const char *LngNE,
	      StringContainer &ngaRecordNames)
{
  char *resHTML =0;
  char *ngaHTML = 0;

  int trialID;

  char peerPage[256];
  char noData[1] =  "";

  const char *eq = noData;
  const char *soilType = noData;
  const char *fault =noData;
  const char *magLo =noData;
  const char *magHi =noData;
  const char *distLo =noData;
  const char *distHi =noData;
  const char *vsLo =noData;
  const char *vsHi =noData;
  const char *pgaLo =noData;
  const char *pgaHi =noData;
  const char *latSW =noData;
  const char *latNE =noData;
  const char *lngSW =noData;
  const char *lngNE =noData;

  if (Eq != 0)
    eq = Eq;
  if (SoilType != 0)
    soilType = SoilType;
  if (Fault != 0)
    fault = Fault;
  if (MagLo != 0)
    magLo = MagLo;
  if (MagHi != 0)
    magHi = MagHi;
  if (DistLo != 0)
    distLo = DistLo;
  if (DistHi != 0)
    distHi = DistHi;
  if (VsLo != 0)
    vsLo = VsLo;;
  if (VsHi != 0)
    vsHi = VsHi;
  if (PgaLo != 0)
    pgaLo = PgaLo;
  if (PgaHi != 0)
    pgaHi = PgaHi;
  if (LatSW != 0)
    latSW =LatSW;
  if (LatNE != 0)
    latNE =LatNE;
  if (LngSW != 0)
    lngSW =LngSW;
  if (LngNE != 0)
    lngNE =LngNE;

  if (strlen(latNE) == 0)
    sprintf(peerPage,"/nga/search?qid=&fault=%s&mag_lo=%s&mag_hi=%s&dist_lo=%s&dist_hi=%s&vs30_lo=%s&vs30_hi=%s&pga_lo=%s&pga_hi=%s&loc=&format=&latSW=&latNE=&lngSW=&lngNE=",
	    fault, magLo, magHi, distLo, distHi, vsLo, vsHi, pgaLo, pgaHi);
  else 
    sprintf(peerPage,"/nga/search?qid=&fault=%s&mag_lo=%s&mag_hi=%s&dist_lo=%s&dist_hi=%s&vs30_lo=%s&vs30_hi=%s&pga_lo=%s&pga_hi=%s&loc=&format=&latSW=%s&latNE=%s&lngSW=%s&lngNE=%s",
	    fault, magLo, magHi, distLo, distHi, vsLo, vsHi, pgaLo, pgaHi, latSW, latNE, lngSW, lngNE);


  fprintf(stderr,"PeerNGA - 1\n");


  trialID = peerGET(peerPage, &resHTML);

  if (trialID == 0 && resHTML != 0) {    

    char *loc = strstr(resHTML,"No matching records were found");
    if (loc != 0) {
      // no records found
      fprintf(stderr, "PeerNGA - no records found\n");
      delete [] resHTML;
      return 0;
    }




    loc = resHTML;
    
    //
    // now get each record from search & get record locations
    //

    while (loc != 0) {
      loc = strstr(loc,"http://peer.berkeley.edu/nga/data?doi=");
      if (loc != 0) {
	loc = strstr(loc,"NGA"); 
	char *endName = strstr(loc,"\n");
	int  ngaNameLength = endName-loc;
	char *ngaName = new char [ngaNameLength+1];
	strncpy(ngaName, loc, ngaNameLength);
	strcpy(&ngaName[ngaNameLength],"");
	sprintf(peerPage,"/nga/data?doi=%s",ngaName);

	fprintf(stderr,"PeerNGA - %s\n", ngaName);

	trialID = peerGET(peerPage, &ngaHTML);
	
	if (ngaHTML != 0) {

	  int ok = 0;

	  //
	  // if soil specified, get soil type & check
	  //

	  if (soilType != noData) {
	    char *siteSoilType = strstr(ngaHTML,"Geomatrix 3:");
	    if (siteSoilType != 0) {
	      siteSoilType += 15;

	      if (siteSoilType[0] != soilType[0]) {
		ok = 1;
	      }
	    }
	  }

	  //
	  // get first 2 recordings
	  // ignoring vertical as some given and some not
	  //
	   
	  if (ok == 0) {

	    char *recordLoc = strstr(ngaHTML,"/nga_files/ath/");
	    
	    //
	    // get the first 2 records
	    //
	    
	    if (recordLoc != 0) {
	      recordLoc+=14;
	      char *recordEndLoc = strstr(recordLoc,".AT2");
	      int recordLength = recordEndLoc - recordLoc;
	      char *recordName = new char[recordLength+1];
	      strncpy(recordName, recordLoc, recordLength);
	      strcpy(&recordName[recordLength],"");
	      ngaRecordNames.addString(recordName);
	      delete [] recordName;
	    } 
	    
	    recordLoc = strstr(recordLoc,"/nga_files/ath/");
	    if (recordLoc != 0) {
	      recordLoc+=14;
	      char *recordEndLoc = strstr(recordLoc,".AT2");
	      int recordLength = recordEndLoc - recordLoc;
	      char *recordName = new char[recordLength+1];
	      strncpy(recordName, recordLoc, recordLength);
	      strcpy(&recordName[recordLength],"");
	      ngaRecordNames.addString(recordName);
	      delete [] recordName;
	    } 
	  }
	  delete [] ngaName;
	  free(ngaHTML);
	}
      }
    }


  //
  // clean up memory
  //
  
  if (resHTML != 0)
    free(resHTML);
  }
  fprintf(stderr,"PeerNGA - DONE\n");

  
  return trialID;
}

