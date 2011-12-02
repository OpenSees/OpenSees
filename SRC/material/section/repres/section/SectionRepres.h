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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/section/SectionRepres.h,v $
                                                                        
                                                                        
// File: SectionRepres.h
//
// Written by Remo M. de Souza
// November 1998


// Purpose: This file contains the class definition for SectionRepres.  
// SectionRepres is an abstract base class and thus no objects of it's
// type can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.


#ifndef SectionRepres_h 
#define SectionRepres_h 

#include <TaggedObject.h>


class SectionRepres: public TaggedObject
{
  public:

    // Section creation functions

    SectionRepres(int tag);    
        
    // Section edition functions

    virtual ~SectionRepres();
   
    // Section inquiring functions
     
    virtual int  getType(void) const = 0;
    friend ostream &operator<<(ostream &s, const SectionRepres &sectionRepres);    
    
  protected:
    
  private:
};


#endif

