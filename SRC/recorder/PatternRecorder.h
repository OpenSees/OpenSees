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
                                                                        
// $Revision: 1.5 $
// $Date: 2004-11-24 22:45:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/PatternRecorder.h,v $
                                                                        
#ifndef PatternRecorder_h
#define PatternRecorder_h

// Written: MHS
// Created: 2002
//
// Description: This file contains the class definition for 
// PatternRecorder. A PatternRecorder records loads values from
// a LoadPattern.


#include <Recorder.h>
#include <ID.h>
#include <Vector.h>

#include <fstream>
using std::ofstream;
class Domain;
class FE_Datastore;

class PatternRecorder: public Recorder
{
  public:
    PatternRecorder(int thePattern,
		 Domain &theDomain,
		 const char *fileName,
		 double deltaT = 0.0,
		 double relDeltaTTol = 0.00001,
		 int startFlag = 0); 

    ~PatternRecorder();
    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    int restart(void);    
    int flush(void);    
    
  protected:
    
  private:	
	int thePattern;
    Domain *theDomain;

    int flag;   // flag indicating whether time, load factor or nothing printed
	        // at start of each line in file
    char *fileName;
    ofstream theFile;     

    double deltaT;
    double relDeltaTTol;
    double nextTimeStampToRecord;
};

#endif
