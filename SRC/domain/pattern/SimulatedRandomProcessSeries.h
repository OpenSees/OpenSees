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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/SimulatedRandomProcessSeries.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), February 2002
// Revised: 
//

#ifndef SimulatedRandomProcessSeries_h
#define SimulatedRandomProcessSeries_h

#include <TimeSeries.h>
#include <Spectrum.h>
#include <RandomNumberGenerator.h>

class SimulatedRandomProcessSeries : public TimeSeries
{
public:
  SimulatedRandomProcessSeries(int tag,
			       RandomNumberGenerator *theRandNumGenerator,
			       Spectrum *theSpectrum,
			       int numFreqIntervals,
			       double mean);

    ~SimulatedRandomProcessSeries();
    
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
    
protected:
	
private:
	RandomNumberGenerator *theRandomNumberGenerator;
	Spectrum *theSpectrum;
	int numFreqIntervals;
	double mean;
	double deltaW;
	Vector *theta;
	Vector *A;
};

#endif
