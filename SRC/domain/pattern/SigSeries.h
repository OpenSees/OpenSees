
#ifndef SigSeries_h
#define SigSeries_h

// Written: Codi McKee (Texas A&M University) 
// Created: 14-Jun-2021


#include <TimeSeries.h>

class SigSeries : public TimeSeries
{
public:
    // constructors
    SigSeries(int tag,
        double tStart, 
        double tFinish, double cFactor, double c1, double c2);

    SigSeries();

    // destructor
    ~SigSeries();

    TimeSeries *getCopy(); 

    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return tFinish-tStart;}
    double getPeakFactor () {return cFactor;}
    double getTimeIncr (double pseudoTime) {return tFinish-tStart;}
    double getStartTime() { return 0.0; } // dummy function
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    double tStart{};      // start time of time series (sec)
    double tFinish{};     // end time of time series (sec)
    double c1{};          // c1 factor
    double c2{};          // c2 factor
    double cFactor{};     // amplitude of sig series
};

#endif
