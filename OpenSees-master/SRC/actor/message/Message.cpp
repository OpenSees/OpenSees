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
                                                                        
// $Revision: 1.6 $
// $Date: 2004-11-25 00:54:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/message/Message.cpp,v $
                                                                        
                                                                        
// File: ~/actor/Message.C
//

#include "Message.h"
#include <OPS_Globals.h>

#ifndef _WIN32
#include <strings.h>
#endif

Message::Message()
:length(0),data(0)
{

}

Message::Message(double *Ptr, int size)
{
    length = size*sizeof(double);
    data = (char *)Ptr;
}

Message::Message(int *Ptr, int size)
{
    length = size*sizeof(int);
    data = (char *)Ptr;
}

Message::Message(char *Ptr, int size)
{
    length = size*sizeof(char);
    data = Ptr;
}

Message::~Message()
{

}

#ifdef _WIN32
extern "C" void bcopy(const char *scource, char *dest, int length);
#include <winsock.h>
#define bcopy(s,d,l) memmove(d,s,l)
#endif

int
Message::putData(char *theData, int startLoc, int endLoc)
{
    if (startLoc > 0 && startLoc <= length &&
	endLoc <= length && endLoc > startLoc) {
	int theLength = endLoc - startLoc;
	char *dataPos = &data[startLoc];
	bcopy(theData, dataPos, theLength);
	return 0;
    } else {
	opserr << "Message::putData() - invalid length of data given\n";
	return -1;
    }
}

void
Message::setData(char *theData, int theLength)
{
	data = theData;
	length = theLength;
}

const char *
Message::getData(void)
{
    return data;
}

int 
Message::getSize()
{
    return length;
}


