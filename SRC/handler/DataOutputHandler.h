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
// $Date: 2004-11-24 22:39:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputHandler.h,v $

#ifndef _DataOutputHandler
#define _DataOutputHandler

#include <MovableObject.h>

class Vector;

class DataOutputHandler: public MovableObject
{
 public:
  DataOutputHandler(int classTag);
  virtual ~DataOutputHandler();

  virtual int open(char **dataDescription, int numData) =0;
  virtual int write(Vector &data) =0;
};

#endif
