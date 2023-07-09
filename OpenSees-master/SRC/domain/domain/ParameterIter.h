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
// $Date: 2006-12-13 18:17:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/ParameterIter.h,v $

// Description: This file contains the class definition for ParameterIter.
// ParameterIter is an abstract base class. An ParameterIter is an iter for 
// returning the elements of an object of class  Domain. ParameterIters 
// must be written for each subclass of Domain (this is done for 
// efficiency reasons), hence the abstract base class.

#ifndef ParameterIter_h
#define ParameterIter_h

class Parameter;

class ParameterIter
{
  public:
    ParameterIter() {};
    virtual ~ParameterIter() {};

    virtual Parameter *operator()(void) =0;

  protected:
    
  private:

};

#endif





