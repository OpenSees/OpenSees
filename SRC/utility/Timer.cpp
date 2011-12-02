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
// $Date: 2003-02-14 23:02:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/Timer.cpp,v $
                                                                        
                                                                        
// File: ~/utility/Timer.C
//
// Written: fmk 
// Created: Mar 1997
// Revision: A
//
// Description: This file contains the class definition for Timer.
// Timer is a stopwatch.
//
// What: "@(#) Timer.h, revA"

#include<Timer.h>
#include <bool.h>

#ifdef NOW
extern "C" int getrusage(int who, struct rusage *rusage);
#endif

#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

Timer::Timer() 
{
#ifdef _WIN32
    // fill in later
#else    
    r1us = &r1usage;
    r2us = &r2usage;
#endif
}

Timer::~Timer()
{ 

}

void 
Timer::start(void)
{
#ifdef _WIN32
    // fill in later
#else        
    t1 = times(&tmsstart);
    getrusage(0,r1us);
#endif
}

void 
Timer::pause(void)
{
#ifdef _WIN32
    // fill in later
#else        
    t2 = times(&tmsend);
    getrusage(0,r2us);    
#endif
} 


double
Timer::getReal(void) const
{
#ifdef _WIN32
    // fill in later
    return 0.0;
#else        
    long clktck = CLK_TCK;    
    double Real = (t2-t1)/(double) clktck;
    return Real;
#endif
}    

double
Timer::getCPU(void) const
{
#ifdef _WIN32
    // fill in later
    return 0.0;
#else        
    long clktck = CLK_TCK;
    double CPU  = (tmsend.tms_utime - tmsstart.tms_utime)/(double) clktck;    
    return CPU;
#endif    
}    

int
Timer::getNumPageFaults(void) const
{
#ifdef _WIN32
    // fill in later
    return 0;
#else        
    int r2yes = r2us->ru_majflt;
    int r1yes = r1us->ru_majflt;
    return r2yes-r1yes;
#endif
}    



void 
Timer::Print(OPS_Stream &s) const
{
#ifdef _WIN32
    // fill in later
#else        
    long clktck = CLK_TCK;
    double Real = (t2-t1)/(double) clktck;
    double CPU  = (tmsend.tms_utime - tmsstart.tms_utime)/(double) clktck;
    double System  = (tmsend.tms_stime - tmsstart.tms_stime)/(double) clktck;
    s << endln;
    s << "TIME(sec) Real: " << Real << "  CPU: " << CPU;
    s << "   System: " << System << endln;

    int r2no = r2us->ru_minflt;
    int r2yes = r2us->ru_majflt;
    int r1no = r1us->ru_minflt;
    int r1yes = r1us->ru_majflt;
    int r1page = r1no + r1yes;
    int r2page = r2no + r2yes;
    
    s << "PAGE FAULTS: " << r2page-r1page << " (NO i/o: ";
    s << r2no-r1no << " YES i/o " << r2yes-r1yes << ") ";

    r2no = r2us->ru_nivcsw;
    r2yes = r2us->ru_nvcsw;
    r1no = r1us->ru_nivcsw;
    r1yes = r1us->ru_nvcsw;
    r1page = r1no + r1yes;
    r2page = r2no + r2yes;

    s << "CONTEXT SWITCHES " << r2page-r1page << " (Invol: ";
    s << r2no-r1no << " Voluntary " << r2yes-r1yes << ") ";

    r2no = r2us->ru_nswap;
    r1no = r1us->ru_nswap;
    r2yes = r2us->ru_maxrss;
    
    s << "Swapped: " << r2no-r1no << " Max Res Set Size: " << r2yes << endln;
    s << endln;
#endif    
}    

OPS_Stream &operator<<(OPS_Stream &s, const Timer &E)
{
    E.Print(s);
    return s;
}


