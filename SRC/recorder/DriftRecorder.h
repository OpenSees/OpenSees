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
// $Date: 2004-01-29 23:30:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/DriftRecorder.h,v $
                                                                        
#ifndef DriftRecorder_h
#define DriftRecorder_h

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for 
// DriftRecorder. 

#include <Recorder.h>
#include <ID.h>

#include <fstream>
using std::ofstream;
class Domain;

class DriftRecorder: public Recorder
{
 public:
  DriftRecorder(int ndI, int ndJ, int dof, int perpDirn,
		Domain &theDomain, const char *fileName, int startFlag = 0); 
  
  ~DriftRecorder();
  int record(int commitTag, double timeStamp);
  int playback(int commitTag);
  void restart(void);    
  
 protected:
  
 private:	
  int ndI;
  int ndJ;
  int dof;
  int perpDirn;
  double oneOverL;
  Domain *theDomain;
  int flag;   // flag indicating whether time, load factor or nothing printed
  // at start of each line in file
  char *theFileName;
  ofstream theFile;     
};

#endif
