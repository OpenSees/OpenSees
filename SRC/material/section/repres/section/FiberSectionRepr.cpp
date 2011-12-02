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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/section/FiberSectionRepr.cpp,v $
                                                                        
                                                                        
// File: FiberSectionRepr.C
// Written by Remo M. de Souza
// November 1998

#include <iostream.h>
#include <stdlib.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber2d.h>
#include <FiberSection.h>

#include <ReinfBar.h>
#include <Cell.h>
#include <QuadPatch.h>
#include <ReinfLayer.h>
#include <FiberSectionRepr.h>


FiberSectionRepr::FiberSectionRepr(int sectionID, int maxNumPatches, int maxNumReinfLayers):
                                   SectionRepres(sectionID), 
                                   sectID(sectionID), 
				   maxNPatches (maxNumPatches),
				   maxNReinfLayers (maxNumReinfLayers),
                                   nPatches(0), 
                                   nReinfLayers(0), 
                                   patch(0), reinfLayer(0)
{
   patch      = new Patch*[maxNumPatches];
   if (!patch)
   {
      cerr << "FATAL ERROR: FiberSectionRepr - not enough memory to allocate " << maxNumPatches << " patches";
      exit (-1);
   }
   int i;
   
   for (i=0; i< maxNumPatches; i++)
       patch[i] = 0;
      
   reinfLayer = new ReinfLayer*[maxNumReinfLayers];
   if (!reinfLayer)
   {
      cerr << "FATAL ERROR: FiberSectionRepr - not enough memory to allocate " << maxNumReinfLayers << "reinforcing layers";
      exit (-1);
   }

   for (i=0; i< maxNumReinfLayers; i++)
       reinfLayer[i] = 0;
      
}


FiberSectionRepr::FiberSectionRepr(int sectionID):
                                   SectionRepres(sectionID), 
                                   sectID(sectionID), 
				   maxNPatches(0),
				   maxNReinfLayers(0),
                                   nPatches(0), 
                                   nReinfLayers(0), 
                                   patch(0), reinfLayer(0)
{
   cerr << "Function FiberSectionRepr::FiberSectionRepr not implemented yet";
   exit (-1);
}


FiberSectionRepr::FiberSectionRepr(int sectionID, 
                                   int numPatches, Patch **patches,
                                   int numReinfLayers,
                                   ReinfLayer **reinfLayers):
                                   SectionRepres(sectionID), 
                                   sectID(sectionID), 
				   maxNPatches(numPatches),
				   maxNReinfLayers(numReinfLayers),
                                   nPatches(numPatches), 
                                   nReinfLayers(numReinfLayers), 
                                   patch(patches), reinfLayer(reinfLayers)
{
   cerr << "Function FiberSectionRepr::FiberSectionRepr not implemented yet";
   exit (-1);
}



FiberSectionRepr::~FiberSectionRepr(void)
{
   int i;

   if (patch)
   {
      for (i = 0; i < maxNPatches; i++)
         if (patch[i])
              delete patch[i];

      delete [] patch;
   }

   if (reinfLayer)
   { 
      for (i = 0; i < maxNReinfLayers; i++)
         if (reinfLayer[i])
            delete reinfLayer[i];
      
      delete [] reinfLayer;
   } 
}
        
 
void FiberSectionRepr::setNumPatches (int numPatches)
{
   cerr << "Function FiberSectionRepr::setnumPatches not implemented yet";
   exit (-1);
}

int FiberSectionRepr::setPatches (Patch **patches)
{
   cerr << "Function FiberSectionRepr::setPatches not implemented yet";
   exit (-1);

   return 1;
}

void FiberSectionRepr::setNumReinfLayers (int numReinfLayers)
{
   cerr << "Function FiberSectionRepr::setnumReinfLayers not implemented yet";
   exit (-1);
}

int FiberSectionRepr::setReinfLayers (ReinfLayer **reinfLayers)
{
   cerr << "Function FiberSectionRepr::setReinfLayers not implemented yet";
   exit (-1);

   return 1;
}
 
int FiberSectionRepr::addPatch (const Patch & aPatch)
{
   int error = 0;

   if (nPatches < maxNPatches)
      patch[nPatches++] = aPatch.getCopy();
   else
      error = 1;

   return error;   
}

int FiberSectionRepr::addReinfLayer (const ReinfLayer & aReinfLayer)
{
   int error = 0;

   if (nReinfLayers < maxNReinfLayers)
      reinfLayer[nReinfLayers++] = aReinfLayer.getCopy();
   else
      error = 1;

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


void FiberSectionRepr::Print(ostream &s, int flag)
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
   
    
ostream &operator<<(ostream &s, FiberSectionRepr &fiberSectionRepr)
{
   fiberSectionRepr.Print(s);
   return s;
}
  
