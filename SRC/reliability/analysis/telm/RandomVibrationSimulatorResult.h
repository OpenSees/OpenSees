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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomVibrationSimulatorResult.h,v $

#ifndef RandomVibrationSimulatorResult_h
#define RandomVibrationSimulatorResult_h

#include <Matrix.h>
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

class RandomVibrationSimulatorResult
{
  public:
	RandomVibrationSimulatorResult(int passedlsf, int passednumtimepoints, 
								   int passenumfragility, double passedeps);

	~RandomVibrationSimulatorResult();
	
	int updateq(int itimepoint, int ifragility, double qvalue);
	int checkconvergence(int itime, int ifrag,int numsimulation);

	void print1(ofstream& output);
	void print2(ofstream& output);
	
  protected: 
    
  private:
	int lsf;
	int numTimePoints;
	int numFragility;
	Matrix* sum_of_q;
	Matrix* sum_of_qsquared;
	Matrix* probability;
	Matrix* cov;
	Matrix* probability_at_convergence;
	int**   convergenceFlag;
	int**   num_of_simulation_at_convergence;
	double eps;
};

#endif
