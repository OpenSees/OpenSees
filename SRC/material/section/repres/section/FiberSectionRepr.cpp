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
// $Date: 2005-03-25 00:32:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/section/FiberSectionRepr.cpp,v $
                                                                        
                                                                        
// File: FiberSectionRepr.C
// Written by Remo M. de Souza
// November 1998

#include <stdlib.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber2d.h>
#include <FiberSection.h>

#include <ReinfBar.h>
#include <Cell.h>
#include <QuadPatch.h>
#include <ReinfLayer.h>
#include <FiberSectionRepr.h>


FiberSectionRepr::FiberSectionRepr(int sectionID, int maxNumPatches, int maxNumReinfLayers)
  :SectionRepres(sectionID), 
  sectID(sectionID), 
  maxNPatches(maxNumPatches),
  maxNReinfLayers(maxNumReinfLayers),
  patch(0), 
  reinfLayer(0), 
  nPatches(0), 
  nReinfLayers(0), 
  numFibers(0), 
  theFibers(0), 
  sizeFibers(32)
{

    theFibers = new Fiber *[sizeFibers];

    if (theFibers == 0) {
      opserr << "FiberSectionRepr::FiberSectionRepr -- failed to allocate Fiber pointers\n";
			      
	sizeFibers = 0;
    }

   patch      = new Patch*[maxNumPatches];
   if (!patch)
   {
      opserr << "FATAL ERROR: FiberSectionRepr - not enough memory to allocate " << maxNumPatches << " patches";
      exit (-1);
   }
   int i;
   
   for (i=0; i< maxNumPatches; i++)
       patch[i] = 0;

   reinfLayer = new ReinfLayer*[maxNumReinfLayers];
   if (!reinfLayer)
   {
      opserr << "FATAL ERROR: FiberSectionRepr - not enough memory to allocate " << maxNumReinfLayers << "reinforcing layers";
      exit (-1);
   }

   for (i=0; i< maxNumReinfLayers; i++)
       reinfLayer[i] = 0;
      
}


FiberSectionRepr::FiberSectionRepr(int sectionID)
  :SectionRepres(sectionID), 
  sectID(sectionID), 
  maxNPatches(0),
  maxNReinfLayers(0),
  patch(0), reinfLayer(0),
  nPatches(0), 
  nReinfLayers(0)
{
   opserr << "Function FiberSectionRepr::FiberSectionRepr not implemented yet";
   exit (-1);
}


FiberSectionRepr::FiberSectionRepr(int sectionID, 
                                   int numPatches, Patch **patches,
                                   int numReinfLayers,
                                   ReinfLayer **reinfLayers)
 :SectionRepres(sectionID), 
  sectID(sectionID), 
  maxNPatches(numPatches),
  maxNReinfLayers(numReinfLayers),
  patch(patches), reinfLayer(reinfLayers),
  nPatches(numPatches), 
  nReinfLayers(numReinfLayers)
{
   opserr << "Function FiberSectionRepr::FiberSectionRepr not implemented yet";
   exit (-1);
}



FiberSectionRepr::~FiberSectionRepr(void)
{
   int i;

   if (patch)
   {
      for (i = 0; i < maxNPatches; i++)
         if (patch[i] != 0)
              delete patch[i];

      delete [] patch;
   }

   if (reinfLayer)
   { 
      for (i = 0; i < maxNReinfLayers; i++)
	if (reinfLayer[i] != 0) 
            delete reinfLayer[i];

      delete [] reinfLayer;
   } 
   
   if (theFibers != 0)
       delete [] theFibers;  // NOTE: don't delete fiber objects themselves
                             //       leave this to FiberSection destructor
}
        
 
void FiberSectionRepr::setNumPatches (int numPatches)
{
   opserr << "Function FiberSectionRepr::setnumPatches not implemented yet";
   exit (-1);
}

int FiberSectionRepr::setPatches (Patch **patches)
{
   opserr << "Function FiberSectionRepr::setPatches not implemented yet";
   exit (-1);

   return 1;
}

void FiberSectionRepr::setNumReinfLayers (int numReinfLayers)
{
   opserr << "Function FiberSectionRepr::setnumReinfLayers not implemented yet";
   exit (-1);
}

int FiberSectionRepr::setReinfLayers (ReinfLayer **reinfLayers)
{
   opserr << "Function FiberSectionRepr::setReinfLayers not implemented yet";
   exit (-1);

   return 1;
}
 
int FiberSectionRepr::addPatch (const Patch & aPatch)
{
   int error = 0;

   if (nPatches < maxNPatches)
      patch[nPatches++] = aPatch.getCopy();

   // out of room: make a new bigger array, copy old over, delete old array and then add new patch
   else {
     maxNPatches *= 2;
     Patch **patches = new Patch*[maxNPatches];     
     if (patches == 0) {
       opserr << "FiberSectionRepr::addPatch() - out of memory\n";
       return 1;
     }
     for (int i=0; i<nPatches; i++)
       patches[i] = patch[i];
     for (int j=nPatches; j<maxNPatches; j++)
       patches[j] = 0;
     
     delete [] patch;
     patch = patches;
     patch[nPatches++] = aPatch.getCopy();
   }

   return error;   
}

int FiberSectionRepr::addReinfLayer (const ReinfLayer & aReinfLayer)
{
   int error = 0;

   if (nReinfLayers < maxNReinfLayers)
      reinfLayer[nReinfLayers++] = aReinfLayer.getCopy();


   // out of room: make a new bigger array, copy old over, delete old array and then add new layer
   else {
     maxNReinfLayers *= 2;
     ReinfLayer **reinfLayers = new ReinfLayer*[maxNReinfLayers];     
     if (reinfLayers == 0) {
       opserr << "FiberSectionRepr::addReinLayer() - out of memory\n";
       return 1;
     }
     for (int i=0; i<nReinfLayers; i++)
       reinfLayers[i] = reinfLayer[i];
     for (int j=nReinfLayers; j<maxNReinfLayers; j++)
       reinfLayers[j] = 0;

     delete [] reinfLayer;
     reinfLayer = reinfLayers;
     reinfLayer[nReinfLayers++] = aReinfLayer.getCopy();
   }

   return error;   
}


int FiberSectionRepr::getType (void) const
{
   return SEC_TAG_FiberSection;
}

int FiberSectionRepr::getNumPatches (void) const
{
   return nPatches;
}

int FiberSectionRepr::getNumReinfLayers (void) const
{
   return nReinfLayers;
}

Patch ** 
FiberSectionRepr::getPatches (void) const
{
   return patch;
}

ReinfLayer ** 
FiberSectionRepr::getReinfLayers (void) const
{
   return reinfLayer;
}


void FiberSectionRepr::Print(OPS_Stream &s, int flag)
{
//   int i;
   
   s << "\nSection representation type: Fiber Section";
   s << "\nMaximum Number of patches: "       << maxNPatches;
   s << "\nMaximum Number of reinf. layers: " << maxNReinfLayers;
   s << "\nCurrent Number of patches: "       << nPatches;
    s << "\nCurrent Number of reinf. layers: " << nReinfLayers;

//   for (i=0; i<nPatches; i++)
//      s << "\nPatch "<<i<<" :" << *patch[i];
   
//   for (i=0; i<nReinfLayers; i++)
//     s << "\nReinfLayer "<<i<<" :" << *reinfLayer[i];

}
   
    
OPS_Stream &operator<<(OPS_Stream &s, FiberSectionRepr &fiberSectionRepr)
{
   fiberSectionRepr.Print(s);
   return s;
}
  


int
FiberSectionRepr::addFiber(Fiber &newFiber)
{
    if (numFibers < sizeFibers) {
	// space available in array .. set new pointer and increment number
	theFibers[numFibers] = &newFiber;
	numFibers++;
    }
    else {
	// need to create a larger array
	int newSize = 2*numFibers;
	if (newSize == 0) 
	    newSize = 2; // in case failed in constructor
	
	Fiber **newArray = new Fiber *[newSize]; 
	
	if (newArray == 0) {
	  opserr << "FiberSection::addFiber -- failed to allocate Fiber pointers\n";
	    return -1;
	}
	    
	// set the new size of the array
	sizeFibers = newSize;
	
	// copy the old pointers
	for (int i = 0; i < numFibers; i++)
	    newArray[i] = theFibers[i];
	
	// add the new pointer
	newArray[numFibers] = &newFiber;
	numFibers++;
	
	// zero the last elements of the array
	for (int j = numFibers; j < newSize; j++) 
	    newArray[j] = 0;
	
	delete [] theFibers;
	
	theFibers = newArray;
    }
  
    return 0;
}


int
FiberSectionRepr::getNumFibers(void) const
{
    return numFibers;
}

Fiber **
FiberSectionRepr::getFibers(void) const
{
    return theFibers;
}
