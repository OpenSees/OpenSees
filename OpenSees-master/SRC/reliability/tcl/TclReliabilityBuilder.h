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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-03-04 00:46:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/tcl/TclReliabilityBuilder.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef TclReliabilityBuilder_h
#define TclReliabilityBuilder_h

#include <tcl.h>
#include <ReliabilityDomain.h>

int inputCheck();

class TclReliabilityBuilder
{
public:
	TclReliabilityBuilder(Domain &theDomain, Tcl_Interp *interp);
	~TclReliabilityBuilder();
    
	ReliabilityDomain *getReliabilityDomain();
    
protected:

private:
    Tcl_Interp *theInterp;
};

#endif







