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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/DiscretizedRandomProcessSeries.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), February 2002
// Revised: 
//

#ifndef DiscretizedRandomProcessSeries_h
#define DiscretizedRandomProcessSeries_h

#include <TimeSeries.h>
#include <ModulatingFunction.h>

class DiscretizedRandomProcessSeries : public TimeSeries
{
public:
  DiscretizedRandomProcessSeries(int tag,
				 int num, 
				 ModulatingFunction **theModFuncs,
				 double p_mean,
				 double targetStdv);
    ~DiscretizedRandomProcessSeries();

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
    
// AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter  (int parameterID, Information &info);
    int activateParameter(int parameterID);
    double getFactorSensitivity(double time);
// AddingSensitivity:END ////////////////////////////////////////////

protected:
	
private:
    int numModFuncs;
	double c;
	double mean;
	double maxStdv;
	double maxStdvTime;
    ModulatingFunction **theModulatingFunctions;
	Vector *randomVariables;
	Vector *kickInTimes;
	int parameterID;
};

#endif
