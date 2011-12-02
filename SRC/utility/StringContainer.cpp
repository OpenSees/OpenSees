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
// $Source: /usr/local/cvs/OpenSees/SRC/utility/StringContainer.cpp,v $
//
// Description: This file contains the class definition for StringContainer.
// SimulationInformation is a stopwatch.
//
// What: "@(#) SimulationContainer.h, revA"



#include<StringContainer.h>
#include <OPS_Globals.h>
#include <string.h>
#include <time.h>
#include <math.h>

StringContainer::StringContainer()
  :strings(0), numStrings(0)
{
}

StringContainer::~StringContainer()
{
  this->clear();
}


const char *
StringContainer::getString(int num) const
{
  if (num >= 0 && num < numStrings)
    return strings[num];

  return 0;
}


int
StringContainer::getNumStrings(void) const
{
  return numStrings;
}

void
StringContainer::clear(void)
{
  if (strings != 0) {
    for (int i=0; i<numStrings; i++)
      delete [] strings[i];
    delete [] strings;
  }

  numStrings = 0;
  strings = 0;
}

int 
StringContainer::addString(const char *theString)
{
  // if valid theString
  if (theString == 0) 
    return 0;

  // create new array to hold pointers and copy pointers there
  char **nextStrings = new char *[numStrings+1];
  if (nextStrings == 0) 
    return -2;

  for (int i=0; i<numStrings; i++)
    nextStrings[i] = strings[i];
  
  // create new pointer for new file name and add to end of array
  char *copyString = new char[strlen(theString)+1];
  if (copyString == 0)
    return -3;

  strcpy(copyString, theString);

  nextStrings[numStrings] = copyString;


  // delete old array and reset the array pointer
  if (strings != 0)
    delete [] strings;
  strings = nextStrings;
  
  // increment number of files
  numStrings++;

  // return ok
  return 0;
}
