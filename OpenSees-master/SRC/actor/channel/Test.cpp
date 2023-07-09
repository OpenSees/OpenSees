
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern int httpGet(char const *URL, char const *page, unsigned int port, char **dataPtr);


extern int httpsGet(char const *URL, char const *page, unsigned int port, char **dataPtr);


extern int neesLogin(const char *user,
		     const char *pass,
		     char **cookieRes);

extern int neesSEND(const char *page,
		    const char *cookie,
		    const char *postData, 
		    const char *contentType,
		    bool, 
		    char **resData);


extern int neesSENDTrial(const char *cookie,
			 int projID,
			 int expID,
			 const char *name,
			 const char *title,
			 const char *objective,
			 const char *description);


extern int neesPOSTTrial_File(const char *cookie,
			      int projID,
			      int expID,
			      const char *name);

extern int neesADD_TrialAnalysisFile(const char *cookie,
				     int projID,
				     int expID,
				     int trialID,
				     const char *path,
				     int length,
				     const char *name,
				     const char *description);

extern int neesADD_TrialAnalysisDir(const char *cookie,
				    int projID,
				    int expID,
				    int trialID,
				    const char *path,
				    const char *dirName);


int main(int argc, char **argv) {
  int res;
  char *resData;
  char *user = "fmckenna";
  char *pass = "amk875002";
  char *cookie = 0;

  res = neesLogin(user,pass, &cookie);

  if (res != 0) {
    fprintf(stderr, "Failed to logon\n");
  } else {
    fprintf(stderr, "Success Logging on\nCookie: %s\n", cookie);    
  }

  char *page = "/REST/Project/150/Experiment/197/Trial";
  char *data = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> <central xmlns=\"http://central.nees.org/api\" xmlns:type=\"central\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:tns=\"http://central.nees.org/api\"> <Trial> <name>123456</name> <title>opensees upload test</title> <objective>the test objective</objective> <description>the test description</description> </Trial> </central> ";

  char *contentType = 0;

  //  res = neesPOST(page, cookie, data, contentType, &resData);

  //  res = neesSENDTrial(cookie, 150, 197, "AnotheTrial", "A New Title", "A New Objective", "A New Description");

  /*
  res = neesPOSTTrial_File(cookie, 150, 197, "TrialPost.xml");
  
  if (res > 0) {
    fprintf(stderr, "Success Adding Trial %d\n", res);  
  
    res = neesPOSTTrialFile(cookie, 150, 197, res, "Test.cpp","A File Description");

  
  }

  */
  /*n


  page = "?projid=150&expid=198&trialid=769&tloc=Analysis&section=&view=browser&basepath=/File/NEES-2006-0150/Experiment-2/Trial-1/Analysis&/path=hello/hello2&floc=Mkdir";
  */

  //  neesADD_TrialAnalysisDir(cookie, 150, 710, 3084, "/output/", "output10");
  /*
  page = "?projid=150&expid=710&trialid=3084&action=DisplayTrialAnalysis&view=browser&basepath=/File/NEES-2006-0150//Trial-1/Analysis&path=/&floc=Mkdir";
  data = "projid=150&expid=710&trialid=3084&action=DisplayTrialAnalysis&view=browser&basepath=/NEES-2006-0150//Trial-1/Analysis&path=/&doit=1&newdir=anothernewdirfmk5";
  contentType="application/x-www-form-urlencoded";
  res = neesSEND(page, cookie, data, contentType, true, &resData);
  */  


  //  neesADD_TrialAnalysisDir(cookie, 150, 710, 3096, "/", "output4");
  //  neesADD_TrialAnalysisDir(cookie, 150, 710, 3096, "/output4", "output5");
  char *path = "/home/fmk/OpenSees/SRC/utility";
  
  //  neesADD_TrialAnalysisFile(cookie, 150, 710, 3096, path, strlen(path)-17, "Test.cpp","test.coo description");    

  neesADD_TrialAnalysisFile(cookie, 150, 710, 3096, path, strlen(path)-17, "Test.cpp","test.coo description");    

  /*
neesADD_TrialAnalysisFile(cookie, 150, 198, 769, "/", "curl.txt","test.coo description");    
		   
  page = "?projid=150&expid=198&trialid=769&tloc=Analysis&section=&view=browser&basepath=/File/NEES-2006-0150/Experiment-2/Trial-1/Analysis&path=/&floc=Mkdir";
  data = "projid=150&expid=198&trialid=769&tloc=Analysis&section=&view=browser&basepath=/NEES-2006-0150/Experiment-2/Trial-1/Analysis&path=/&doit=1&newdir=anothernewdirfmk4";
  contentType="application/x-www-form-urlencoded";
  //    res = neesSEND(page, cookie, data, contentType, true, &resData);
  */


  free(cookie);
  //  free(resData);

  fprintf(stdout,"BYE\n");
  
  exit(0);
}
