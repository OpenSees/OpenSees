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

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/TimePoints.h,v $

#ifndef TimePoints_h
#define TimePoints_h

#include <Vector.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class TimePoints
{
  public:
	TimePoints(double,double,double,double,bool);
	~TimePoints();

	int getNumTimePoints(void);
	int MakeOrder(double MirTime);
	int getAnalysisStep(int ii);
	double getAnalysisTime(int ii)
	{ return (*TmPts)(ii);}
	int getOrder(int ii);
	void printTimePts(ofstream&);
	void initializeOrder();
  protected: 
    
  private:
	  void MakeTimePoints(void);
      int StartPoint(void){return anaStps[0];}

 	  int StepsToStart;
	  int StepsToEnd;
	  int SamplingFrequency;
	  int NumPoints;
	  Vector* TmPts;
	  int* anaStps;
	  int* OrdAna;
	  double delta;
	  bool print;
	  ofstream output;

};

#endif
