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
// $Date: 2009-03-27 19:18:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/FileIter.cpp,v $
                                                                        
// Written: fmk 
// Created: 09/07


#include <File.h>
#include <FileIter.h>

FileIter::FileIter(File &theFilePtr)
  
{
  theDirFiles = &(theFilePtr.dirFiles);
  iter = theDirFiles->begin();
}

FileIter::~FileIter()
  
{

}


File *
FileIter::operator()(void)
{

  if (iter != theDirFiles->end()) {
    File *result = iter->second;
    iter++;
    return result;
  } else
    return 0;
}

