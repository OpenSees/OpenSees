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
// $Date: 2010-03-04 19:05:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/TCP_Stream.h,v $


// Written: fmk

#ifndef _TCP_Stream
#define _TCP_Stream

#include <OPS_Stream.h>
#include <Vector.h>
class TCP_Socket;

class TCP_Stream : public OPS_Stream
{
 public:
    TCP_Stream();        
    TCP_Stream(unsigned int other_Port, 
	       const char *other_InetAddr,
	       bool checkEndianness = false); 

    ~TCP_Stream();

    int setFile(const char *fileName, openMode mode = OVERWRITE);
    int open(void);
    int close(void);
    
    // xml stuff
    int tag(const char *);
    int tag(const char *, const char *);
    int endTag();
    int attr(const char *name, int value);
    int attr(const char *name, double value);
    int attr(const char *name, const char *value);
    int write(Vector &data);
    
    // regular stuff
    OPS_Stream& write(const char *s, int n);
    OPS_Stream& write(const unsigned char *s, int n);
    OPS_Stream& write(const signed char *s, int n);
    OPS_Stream& write(const void *s, int n);
    OPS_Stream& operator<<(char c);
    OPS_Stream& operator<<(unsigned char c);
    OPS_Stream& operator<<(signed char c);
    OPS_Stream& operator<<(const char *s);
    OPS_Stream& operator<<(const unsigned char *s);
    OPS_Stream& operator<<(const signed char *s);
    OPS_Stream& operator<<(const void *p);
    OPS_Stream& operator<<(int n);
    OPS_Stream& operator<<(unsigned int n);
    OPS_Stream& operator<<(long n);
    OPS_Stream& operator<<(unsigned long n);
    OPS_Stream& operator<<(short n);
    OPS_Stream& operator<<(unsigned short n);
    OPS_Stream& operator<<(bool b);
    OPS_Stream& operator<<(double n);
    OPS_Stream& operator<<(float n);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
 private:
    int sendSize;
    Vector data;
    TCP_Socket *theChannel;
};

#endif
