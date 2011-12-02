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
// $Date: 2006-11-08 20:06:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/StringContainer.h,v $
                                                                        
// Written: fmk 
// Created: 11/06
//
// Description: This file contains the class definition for StringContainer.
// StringContainer is used to store information about a simulation; this
// includes start and end times; program version, files opened for reading, files
// opened for writing, and all parameters used (specified with pset or -par option
// to program)
//
// What: "@(#) StringContainer.h, revA"

#ifndef StringContainer_h
#define StringContainer_h

class StringContainer
{
 public:
  StringContainer();
  ~StringContainer();
  int addString(const char *);
  const char *getString(int) const;
  const char *operator()(void);
  int getNumStrings() const;
  void clear(void);

 private:
  char **strings;
  int numStrings;
};

#endif
