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
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/MPAccSeries.h,v $


#ifndef MPAccSeries_h
#define MPAccSeries_h

// Written: Tang.S 
// Created: 2023-01
// Revision: A

#include <TimeSeries.h>

class MPAccSeries : public TimeSeries
{
public:
    // constructors
    MPAccSeries(int tag,
        double tStart, 
        double tFinish,
        double period, 
        double gammaMP, 
        double nuMP = 1.0,
        double AFactor = 1.0);

    MPAccSeries();

    // destructor
    ~MPAccSeries();

    TimeSeries *getCopy(); 

    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return tFinish-tStart;}
    double getPeakFactor () {return AFactor;}
    double getTimeIncr (double pseudoTime) {return tFinish-tStart;}

    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    double tStart;      // start time of time series (sec)
    double tFinish;     // end time of time series (sec)
    double period;      // period of MPAcc series (sec)
    double gammaMP;  // ¦Ã factor in M&P pulse model 
    double nuMP;     // nu in degree in M&P pulse model
    double AFactor;   // the M&P velocity amplificarion factor(optional,default=1.0)
};

#endif
