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
// $Date: 2000-12-18 10:03:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.h,v $
                                                                        
                                                                        
#ifndef ElementRecorder_h
#define ElementRecorder_h

// File: ~/recorder/ElementRecorder.h
//
// Written: fmk 
// Created: 09/99
// Revision: A
//
// Description: This file contains the class definition for ElementRecorder.
// A ElementRecorder is used to obtain a response from an element during 
// the analysis.
//
// What: "@(#) ElementRecorder.h, revA"

#include <Recorder.h>
#include <fstream.h>

#include <Information.h>
#include <G3Globals.h>
#include <ID.h>


class Domain;
class Vector;
class Matrix;
class Element;
class Response;

class ElementRecorder: public Recorder
{
  public:
    ElementRecorder(const ID &eleID, Domain &theDomain, char **argv, int argc,
		    bool echoTime, char *fileName =0);

    ~ElementRecorder();
    int record(int commitTag);
    int playback(int commitTag);

    void restart(void);    
    
  protected:
    
  private:	
    int numEle;
    ID responseID;                 // integer element returns in setResponse
    Information	*eleInfoObjects;  // object into whcih element places the response
    Element **theElements;         // pointer to the elements

	Response **theResponses;

    Domain *theDomain;
    bool echoTimeFlag;     // flag indicating if pseudo time also printed
    char theFileName[MAX_FILENAMELENGTH];  // file name  
    ofstream theFile; 	   // output stream
};


#endif
