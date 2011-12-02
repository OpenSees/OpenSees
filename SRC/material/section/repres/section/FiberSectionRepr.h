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
// $Date: 2000-12-13 05:56:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/section/FiberSectionRepr.h,v $
                                                                        
                                                                        
// File: FiberSectionRepr.h
// Written by Remo M. de Souza (November 1998)
// modified by rms (July 1999) - doesn't prespecify number of patches and reinf. layers
//             rms (August 1999) - add fibers to section for the analysis 
		    
#ifndef FiberSectionRepr_h 
#define FiberSectionRepr_h 

#include <SectionRepres.h>

class Patch;
class ReinfLayer;
class Fiber;

#define SEC_TAG_FiberSection 1

class FiberSectionRepr: public SectionRepres
{
  public:

    // constructor and destructor

    FiberSectionRepr(int sectionID);
    FiberSectionRepr(int sectionID, int maxNumPatches, int maxNumReinfLayers);
    FiberSectionRepr(int sectionID, int numPatches, Patch **patches, 
                     int numReinfLayers, ReinfLayer **reinfLayers);

    ~FiberSectionRepr();
        
    // edition functions
    void setNumPatches     (int numPatches);
    int  setPatches        (Patch **patches);    
    void setNumReinfLayers (int numReinfLayers);
    int  setReinfLayers    (ReinfLayer **reinfLayers);
 
    int  addPatch          (const Patch & aPatch);
    int  addReinfLayer     (const ReinfLayer & aReinfLayer);
    int  addFiber(Fiber &theFiber);    

    // inquiring functions
    int  getType(void) const;
      
    int          getNumPatches     (void) const;
    int          getNumReinfLayers (void) const;
    int          getNumFibers(void) const;    
    Patch      **getPatches        (void) const;     
    ReinfLayer **getReinfLayers    (void) const;
    Fiber     **getFibers(void) const;
    
    void Print(ostream &s, int flag =0);
    friend ostream &operator<<(ostream &s, FiberSectionRepr &fiberSectionRepr);    
  

  protected:
    
  private:
    
    int        sectID;
    int        maxNPatches;       // maximum number of patches that can be stored
    int        maxNReinfLayers;   // maximum number of reinforcing layers  that can be stored
    Patch      **patch;           // patch array (change to linked list later!!!!) 
    ReinfLayer **reinfLayer;      // reinforcing bar array (change to linked
                                  // list later!!!!!!!!)
    int        nPatches;          // current number of patches
    int        nReinfLayers;      // current number of reinforcing layers

    int numFibers;       // number of fibers in the section
    Fiber **theFibers;   // array of pointers to fibers
                         // that form the section
    int sizeFibers;      // size of the fibers array
};


#endif

