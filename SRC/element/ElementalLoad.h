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
                                                                        
// $Revision: 1.3 $
// $Date: 2002-06-11 20:48:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ElementalLoad.h,v $
                                                                        
                                                                        
#ifndef ElementalLoad_h
#define ElementalLoad_h

// Written: fmk 
//
// Purpose: This file contains the class definition for ElementalLoad.
// ElementalLoad is an abstract class.

#include <Load.h>
class Element;

class ElementalLoad : public Load
{
  public:
    ElementalLoad(int tag, int classTag, const ID &theElementTags);
    ElementalLoad(int tag, int classTag);
    ElementalLoad(int classTag);    
    ~ElementalLoad();

    virtual void setDomain(Domain *theDomain);
    virtual void applyLoad(double loadfactor);
    virtual const Vector &getData(int &type, double loadFactor) = 0;

    virtual const ID &getElementTags(void);
    virtual int removeElement(int tag); // returns -1 if fails, 
                                        // numElements left if removed
  protected:
    int setElementTags(const ID &theEleTags);
    ID *theElementTags;     // copy of element tags, removed in setDomain

  private:
    Element **theElements;  // pointer to associated elements
    int numElements;        // number of associated elements
};

#endif

