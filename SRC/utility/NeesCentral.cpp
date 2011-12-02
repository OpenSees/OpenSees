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
// $Date: 2007-10-05 22:07:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/NeesCentral.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/07
//
// Description: This file contains the class definition for File used in SimulationINformation.
//
// What: "@(#) SimulationInformation.h, revA"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef _HTTPS

#ifdef _WIN32
extern int __cdecl
#else
extern int
#endif
httpsGet(char const *URL, char const *page, char const *cookie, unsigned int port, char **dataPtr);


#ifdef _WIN32
extern int __cdecl
#else
extern int
#endif
httpsSEND(const char *URL, const char *page, const char *cookie, const char *contentType, 
	  const char *dataToPost, unsigned int port, bool returnHeader, bool doPOST, char **resPtr);
	  

#ifdef _WIN32
extern int __cdecl
#else
extern int
#endif
httpsSEND_File(const char *URL, const char *page, const char *cookie, const char *contentType,
	       const char *filename, unsigned int port, bool returnHeader, bool doPOST, char **resPtr);
	       
#ifdef _WIN32
int __cdecl
#else
int
#endif
neesLogin(const char *user,
	  const char *pass,
	  char **cookieRes)
{
  char *URL="central.nees.org";
  char *loginPage = "/login.php";
  char *loginType="application/x-www-form-urlencoded";
  char *cookie = 0;
  char loginData[512];
  char *postRes = 0;
  
  sprintf(loginData,"redirect=https%s3A%s2F%s2Fcentral.nees.org%s2Findex.php%s3F&user=%s&pass=%s","%","%","%","%","%",user,pass);  

  /*
   * login & get GAsession cookie
   */

  httpsSEND(URL, loginPage, cookie, loginType, loginData,  443, true, true, &postRes);

  // make sure GAsession in cookie
  char *GAsession = strstr(postRes,"GAsession");  
  if (GAsession == 0) {
    fprintf(stderr, "neesLogin: could not login, check username, password or central.nees.org availability\n");
    free(postRes);
    return -1;
  }

  /*
   * set the coookie data for next call
   */

  cookie = (char *)malloc(8*sizeof(char));
  int sizeCookie = 8;
  strcpy(cookie, "Cookie:");

  char *startNextCookie = strstr(postRes,"Set-Cookie");
  while (startNextCookie != 0) {

    startNextCookie += 11; // "Set-Cookie:");
    char *endNextCookie = strstr(startNextCookie, ";");
    int sizeNextCookie = endNextCookie-startNextCookie+1;

    char *nextCookie = (char *)malloc((sizeNextCookie+sizeCookie+1)*sizeof(char));
    
    strncpy(nextCookie, cookie, sizeCookie-1);
    strncpy(&nextCookie[sizeCookie-1], startNextCookie, sizeNextCookie);

    free(cookie);
    cookie = nextCookie;
    sizeCookie += sizeNextCookie;

    strcpy(&cookie[sizeCookie],"");

    startNextCookie++;
    startNextCookie = strstr(startNextCookie,"Set-Cookie");
  }
  strcpy(&cookie[sizeCookie-2],"\n");

  free(postRes);
  *cookieRes = cookie;

  return 0;
}




#ifdef _WIN32
int __cdecl
#else
int
#endif
neesSEND(const char *page,
	 const char *cookie,
	 const char *postData, 
	 const char *contentType,
	 bool doPOST,
	 char **res) {

  char *URL="central.nees.org";

  httpsSEND(URL, page, cookie, contentType, postData,  443, true, doPOST, res);

  char *failure = strstr(*res,"400 Bad Request");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND\n");
    return -1;
  }

  failure = strstr(*res,"401 Unauthorized");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND\n");
    return -2;
  }

  failure = strstr(*res,"404 Not Found");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND\n");
    return -3;
  }

  failure = strstr(*res,"410 Gone");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND\n");
    return -4;
  }

  return 0;
}



#ifdef _WIN32
int __cdecl
#else
int
#endif
neesGET(const char *page,
	const char *cookie,
	char **res) {

  char *URL="central.nees.org";

  httpsGet(URL, page, cookie, 443, res);

  char *failure = strstr(*res,"400 Bad Request");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesGET\n");
    return -1;
  }

  failure = strstr(*res,"401 Unauthorized");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesGET\n");
    return -2;
  }

  failure = strstr(*res,"404 Not Found");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesGET\n");
    return -3;
  }

  failure = strstr(*res,"410 Gone");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesGET\n");
    return -4;
  }

  return 0;  
}

#ifdef _WIN32
int __cdecl
#else
int
#endif
neesSEND_File(const char *page,
	      const char *cookie,
	      const char *filename, 
	      const char *contentType,
	      bool doPOST,
	      char **res) {
  

  char *URL="central.nees.org";

  httpsSEND_File(URL, page, cookie, contentType, filename,  443, true, doPOST, res);

  char *failure = strstr(*res,"400 Bad Request");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND_File\n");
    return -1;
  }

  failure = strstr(*res,"401 Unauthorized");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND_File\n");
    return -2;
  }

  failure = strstr(*res,"404 Not Found");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND_File\n");
    return -3;
  }

  failure = strstr(*res,"410 Gone");

  if (failure != 0) {
    fprintf(stderr, "ERROR: neesSEND_File\n");
    return -4;
  }  

  return 0;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
neesADD_Trial(const char *cookie,
	      int projID,
	      int expID,
	      const char *name,
	      const char *title,
	      const char *objective,
	      const char *description)
{
  char *resHTML =0;
  int trialID;

  char neesPage[96];
  sprintf(neesPage,"/REST/Project/%d/Experiment/%d/Trial%c",projID,expID,'\0');

  char *xmlData = (char *)malloc(512+strlen(name)+strlen(title)+strlen(objective)+strlen(description));

  sprintf(xmlData,"<?xml version=\"1.0\" encoding=\"UTF-8\"?> <central xmlns=\"http://central.nees.org/api\" xmlns:type=\"central\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:tns=\"http://central.nees.org/api\"> <Trial> <name>%s</name> <title>%s</title> <objective>%s</objective> <description>%s</description> </Trial> </central> %c",name, title, objective, description,'\0');

  trialID = neesSEND(neesPage, cookie, xmlData, 0, true, &resHTML);



  if (trialID == 0 && resHTML != 0) {    

      char *startNextCookie = strstr(resHTML, neesPage);
      
      if (startNextCookie != NULL) {
	startNextCookie += strlen(neesPage);
	sscanf(startNextCookie, "/%d", &trialID);
	
      } else {
	trialID = -1;
      }
      
  } else
    trialID = -1;

  //
  // clean up memory
  //

  free(xmlData);
  if (resHTML != 0)
    free(resHTML);

  
  return trialID;
}




#ifdef _WIN32
int __cdecl
#else
int
#endif
neesADD_Trial_XmlFromFromFile(const char *cookie,
			      int projID,
			      int expID,
			      const char *filename)
{
  char *resHTML =0;
  int trialID;
  
  char neesPage[96];
  sprintf(neesPage,"/REST/Project/%d/Experiment/%d/Trial%c",projID,expID,'\0');
  
  trialID = neesSEND_File(neesPage, cookie, filename, 0, true, &resHTML);

  if (trialID == 0 && resHTML != 0) {    
    
      char *startNextCookie = strstr(resHTML, neesPage);
      
      if (startNextCookie != NULL) {
	startNextCookie += strlen(neesPage);
	sscanf(startNextCookie, "/%d", &trialID);
	
      } else {
	trialID = -1;
      }
      
  } else
    trialID = -1;
  
  if (resHTML != 0)
    free(resHTML);
  
  return trialID;
}



#ifdef _WIN32
int __cdecl
#else
int
#endif
neesADD_TrialAnalysisFile(const char *cookie,
			  int projID,
			  int expID,
			  int trialID,
			  const char *path,
			  int basePathLength,
			  const char *filename,
			  const char *description)
{
  char *resHTML =0;
  int res;

  //
  // we first need to get the projects experiments trial details
  // this includes information  about where files are store
  //

  char neesPage[96];
  sprintf(neesPage,"/REST/Project/%d/Experiment/%d/Trial/%d%c",projID,expID,trialID,'\0');

  res = neesGET(neesPage, cookie, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to find Trial\n");    
    free(resHTML);
    return -1;
  }

  //
  // we get information about where the path for analysis part of trial is
  // and we can create the page address which contains this path, the provided path and filename
  //

  char *startData = strstr(resHTML, "DataFile");
  startData+=15; // DataFile link="

  char *endData = strstr(resHTML,"Analysis");
  endData+=8;

  if (startData == 0 || endData == 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to get Analysis dir in Trial\n");    
    free(resHTML);
    return -1;
  }

  char *nextStartData = startData+1;
  while (nextStartData != 0 && nextStartData < endData) {
    nextStartData = strstr(nextStartData, "DataFile");
    nextStartData+=15;
    if (nextStartData < endData)
      startData = nextStartData;
  }

  char *neesFilePage = (char *)malloc(endData-startData+10+strlen(path)+strlen(filename)); // 10 for /content
  strncpy(neesFilePage, startData, endData-startData);
  strncpy(&neesFilePage[endData-startData], "\0", 1);
  strcat(neesFilePage,&path[basePathLength]);
  if (neesFilePage[strlen(neesFilePage)-1] != '/')
    strcat(neesFilePage,"/");
  strcat(neesFilePage,filename);

  free(resHTML);

  //
  // we can create meta data for file, to be sent to page in a POST
  //

  char *xmlData = (char *)malloc(512+strlen(filename)+strlen(description)+strlen(path));  

  sprintf(xmlData,"<?xml version=\"1.0\" encoding=\"UTF-8\"?> <central xmlns=\"http://central.nees.org/api\" xmlns:type=\"central\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:tns=\"http://central.nees.org/api\"> <DataFile> <name>%s</name> <Description>%s</Description> </DataFile> </central> %c",filename, description,'\0');


  res = neesSEND(neesFilePage, cookie, xmlData, 0, true, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to send file meta data\n");    
    free(resHTML);
    return -1;
  }

  free(resHTML);

  //
  // once we have POSTED file xml data we can do a PUT of actual data, add /content to page
  // placing local location of file in xmlData array
  //

  strcat(neesFilePage,"/content");
  strcpy(xmlData, path);
  strcat(xmlData,filename);

  res = neesSEND_File(neesFilePage, cookie, xmlData, 0, false, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to send file contents\n");    
    free(resHTML);
    return -1;
  }

  //
  // clean up memory
  //

  free(neesFilePage);
  free(xmlData);

  if (resHTML != 0)
    free(resHTML);

  return 0;
}


int neesADD_TrialAnalysisDir(const char *cookie,
			     int projID,
			     int expID,
			     int trialID,
			     const char *path,
			     const char *dirName)
{

  char *resHTML =0;
  int res;

  //
  // we first need to get the projects experiments trial details
  // this includes information  about where files are stored
  //

  char neesPage[96];
  sprintf(neesPage,"/REST/Project/%d/Experiment/%d/Trial/%d%c",projID,expID,trialID,'\0');

  res = neesGET(neesPage, cookie, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to find Trial\n");    
    free(resHTML);
    return -1;
  }

  //
  // we get information about where the path for analysis part of trial is
  // and we can create the page address which contains this path, the provided path and filename
  //

  char *startData = strstr(resHTML, "DataFile");
  startData+=15; // DataFile link="

  char *endData = strstr(resHTML,"Analysis");
  endData+=8;

  if (startData == 0 || endData == 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to get Analysis dir in Trial\n");    
    free(resHTML);
    return -1;
  }

  char *nextStartData = startData+1;
  while (nextStartData != 0 && nextStartData < endData) {
    nextStartData = strstr(nextStartData, "DataFile");
    nextStartData+=15;
    if (nextStartData < endData)
      startData = nextStartData;
  }

  char *neesFilePage = (char *)malloc(endData-startData+10+strlen(path)+strlen(dirName)); // 10 for /content
  strncpy(neesFilePage, startData, endData-startData);
  strncpy(&neesFilePage[endData-startData], "\0", 1);
  strcat(neesFilePage,path);
  if (neesFilePage[strlen(neesFilePage)-1] != '/')
    strcat(neesFilePage,"/");

  char *restPath = (char *)malloc(strlen(neesFilePage));
  strcpy(restPath,&neesFilePage[11]);
  strcat(&restPath[strlen(restPath)-2],"");  

  strcat(neesFilePage,dirName);

  free(resHTML);

  //
  // we can create meta data for file, to be sent to page in a POST
  //

  char *contentLink = (char *)malloc(strlen(neesFilePage)+10); 
  strcpy(contentLink, neesFilePage);
  strcat(contentLink,"/content");

  char *xmlData = (char *)malloc(512+strlen(dirName)+strlen(neesFilePage)+strlen(contentLink));  

  sprintf(xmlData,"<?xml version=\"1.0\" encoding=\"UTF-8\"?> <central xmlns=\"http://central.nees.org/api\" xmlns:type=\"central\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:tns=\"http://central.nees.org/api\"> <DataFile isDirectory=\"true\"> <name>%s</name> <path>%s</path> <contentLink>%s</contentLink></DataFile> </central> %c",dirName, restPath, contentLink, '\0');

  res = neesSEND(neesFilePage, cookie, xmlData, 0, true, &resHTML);


  if (res < 0) {
    fprintf(stderr, "IGNORE neesSEND error - currently a bug in central.nees.org nees\n");    
    //    free(resHTML);
    //    return -1;
  }


  free(resHTML);

  //
  // once we have file name we can do a PUT with blank data
  //


  res = neesSEND(neesFilePage, cookie, "", 0, false, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisFile - failed to send file contents\n");    
    free(resHTML);
    return -1;
  }

  //
  // clean up memory
  //

  free(restPath);
  free(neesFilePage);
  free(xmlData);
  free(contentLink);
  if (resHTML != 0)
    free(resHTML);

  return 0;
}

/*
int neesADD_TrialAnalysisDir(const char *cookie,
			     int projID,
			     int expID,
			     int trialID,
			     const char *path,
			     const char *dirName)
{

  // 
  // as of this moment there is nothing in the REST interface for adding a directory
  // we will do it using there browser interface .. takes a couple of gets to get some info
  // and not too stable .. however what can you do!
  //

  char *resHTML =0;
  int res;
  char *startData, *endData;
  char neesPage[1024];

  //
  // first we need to do a GET to get Analysis directory
  //
  sprintf(neesPage,"/REST/Project/%d/Experiment/%d/Trial/%d%c",projID,expID,trialID,'\0');

  res = neesGET(neesPage, cookie, &resHTML);
  
  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisDir - failed to get trial information\n");    
    free(resHTML);
    return -1;
  }

  startData = strstr(resHTML, "DataFile");
  startData+=20; // <DataFile link="/REST

  endData = strstr(resHTML,"Analysis");
  endData+=8;
  
  char *neesFilePage = new char[endData-startData+1];


  strncpy(neesFilePage, startData, endData-startData);
  strncpy(&neesFilePage[endData-startData], "\0", 1);


  //char *expData = strstr(startData,"Experiment");
  //  char *trialData = strstr(startData,"/Trial");
  //strncpy(neesFilePage, startData, expData-startData);
  //strncpy(&neesFilePage[expData-startData], trialData, endData-trialData);
  //strncpy(&neesFilePage[expData-startData+endData-trialData+1], "\0", 1);

  free(resHTML);

  int dataLength = 512+strlen(neesFilePage)+strlen(path);

  //
  // now we can create the directory
  //

  char *neesBrowserPage =  new char[dataLength];
  
  //  sprintf(neesBrowserPage,"?projid=%d&expid=%d&trialid=%d&tloc=Analysis&section=&view=browser&basepath=%s&path=%s&floc=Mkdir",
  sprintf(neesBrowserPage,"?projid=%d&expid=%d&trialid=%d&action=DisplayTrialAnalysis&view=browser&basepath=%s&path=%s&floc=Mkdir",
	  projID, expID, trialID, neesFilePage, path);


  dataLength+= strlen(dirName);

  char *data = new char[dataLength];
  sprintf(data,"projid=%d&expid=%d&trialid=%d&action=DisplayTrialAnalysis&view=browser&basepath=%s&path=%s&doit=1&newdir=%s",
	  projID, expID, trialID, &neesFilePage[5], path, dirName); //&neesFilePage[5] to remove /File 


//  sprintf(data,"projid=%d&expid=%d&trialid=%d&action=DisplayTrialAnalysis&view=browser&basepath=%s&path=%s&doit=1&newdir=%s",
//	  projID, expID, trialID, &neesFilePage[5], path, dirName); //&neesFilePage[5] to remove /File 

  char *contentType="application/x-www-form-urlencoded";

  res = neesSEND(neesBrowserPage, cookie, data, contentType, true, &resHTML);

  if (res < 0) {
    fprintf(stderr, "ERROR neesADD_TrialAnalysisDir - failed to add dir\n");    
    return -1;
  }

  // clean up memory
  delete [] data;
  delete [] neesBrowserPage;
  delete [] neesFilePage;
  free(resHTML);

  return 0;
}
*/

#endif

