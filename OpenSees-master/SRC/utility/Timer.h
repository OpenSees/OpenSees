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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:02:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/Timer.h,v $
                                                                        
                                                                        
// File: ~/utility/Timer.h
//
// Written: fmk 
// Created: Mar 1997
// Revision: A
//
// Description: This file contains the class definition for Timer.
// Timer is a stopwatch.
//
// What: "@(#) Timer.h, revA"

#ifndef Timer_h
#define Timer_h

#ifdef _WIN32

#else
#ifdef _MAC
#include <sys/time.h>
#else

#include <time.h>
#include <unistd.h>

#endif
#include <sys/times.h>
#include <sys/resource.h>
#endif

#include <OPS_Globals.h>

class Timer
{
  public:
    Timer();    
    virtual ~Timer();

    void start(void);
    void pause(void);
    double getReal(void) const;
    double getCPU(void) const;
    int getNumPageFaults(void) const;
    
    virtual void Print(OPS_Stream &s) const;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const Timer &E);    

  protected:
    
  private:

#ifdef _WIN32

// fill in later

#else

#ifdef TIMER_USE_MPIWTIME
    double t1, t2;
    double tmsstart, tmsend;
    // struct rusage r1usage, r2usage;
    // struct rusage *r1us, *r2us;
#else 
    clock_t t1, t2;
    struct tms tmsstart, tmsend;
    struct rusage r1usage, r2usage;
    struct rusage *r1us, *r2us;
#endif //TIMER_USE_MPIWTIME

#endif
};


#endif

