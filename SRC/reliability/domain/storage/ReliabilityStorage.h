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
                                                                        
                                                                        
#ifndef ReliabilityStorage_h
#define ReliabilityStorage_h

// Written: Kevin Mackie

#include <Information.h>

class ReliabilityStorage
{
  public:
    ReliabilityStorage();    
    virtual ~ReliabilityStorage();

    virtual const char *getClassType(void) const;

    virtual int setVariable(const char *variable, Information &);
    virtual int getVariable(const char *variable, Information &);
    
  protected:
    
  private:

    
};

#endif
