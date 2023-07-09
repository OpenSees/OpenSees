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
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/feap/fElmt05.h,v $
                                                                        
                                                                        
#ifndef fElmt05_h
#define fElmt05_h

// File: ~/element/fortran/fElmt05.h
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// Description: This file contains the class definition for fElmt05. fElmt05
// is a wrapper used to call fortran element subroutine elmt05. elmt05 is a 
// linear elastic 2d element
//
// What: "@(#) fElmt05.h, revA"

#include <fElement.h>

class fElmt05 : public fElement
{
  public:
    // constructor
    fElmt05(int tag,
	    int Nd1, int Nd2,
	    double E, double A, double rho = 0.0);
    
    fElmt05();    
    
    // destructor
    ~fElmt05();

    virtual int getNumExternalNodes(void) const;
    virtual const ID &getExternalNodes(void);

    virtual int invokefRoutine(double *d, double *ul, double *xl, int *ix, 
			       double *tl, double *s, double *r, int ndf, 
			       int ndm, int nst, int isw, double dm, int nen, int n, 
			       int nh1, int nh2, int nh3, double *h,  
			       double *ctan, int ior, int iow);

  protected:
	     
  private:

};

#endif



