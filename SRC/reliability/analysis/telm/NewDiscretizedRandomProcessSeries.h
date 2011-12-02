/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.3 $
// $Date: 2010-02-04 18:31:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewDiscretizedRandomProcessSeries.h,v $
                                                                       

#ifndef NewDiscretizedRandomProcessSeries_h
#define NewDiscretizedRandomProcessSeries_h

#include <TimeSeries.h>
#include <ModulatingFunction.h>
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

class NewDiscretizedRandomProcessSeries : public TimeSeries
{
public:
    NewDiscretizedRandomProcessSeries(int num, 
				      ModulatingFunction **theModFuncs,
				      double p_mean,
				      double targetStdv);
    ~NewDiscretizedRandomProcessSeries();

    TimeSeries *getCopy(void);
    // method to get load factor
    double getFactor(double pseudoTime);

    // None of the following functions should be invoked on this type
    // of object
    double getDuration () {return 0.0;} // dummy function
    double getPeakFactor () {return 0.0;} // dummy function
    double getTimeIncr (double pseudoTime) {return 1.0;} // dummy function
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0); 
    int getparameterID() const{ return parameterID;}
    void setparameterID(int param) { parameterID = param ;}
    
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter     (const char **argv, int argc, Parameter &param);
    int updateParameter  (int parameterID, Information &info);
    int activateParameter(int parameterID);
    double getFactorSensitivity(double time);
    double getFactorSensitivity(double time, double ktime);
    double getkickInTimes(int rvNum) const;
    int getPulseSequentialID( int rvNum ) const;
    int getNumPulses() const { return numRandVar; }
    int updateRV (int nrv, double value);
    //	void sortPulse (void);
    // AddingSensitivity:END ////////////////////////////////////////////
    
 protected:
    
 private:
    int numModFuncs;
    double c;
    double mean;
    double maxStdv;
    double maxStdvTime;
    ModulatingFunction **theModulatingFunctions;
    Vector *randomVariables;	// current random pulses in this class;
    Vector *kickInTimes;		// launching time of each random pulses;
    int numRandVar;	  // number of random variable 
    int parameterID;	
    bool* active;
    ofstream output;
};

#endif
