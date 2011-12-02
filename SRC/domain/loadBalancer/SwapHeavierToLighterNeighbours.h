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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/SwapHeavierToLighterNeighbours.h,v $
                                                                        
                                                                        
// File: ~/domain/loadBalancer/SwapHeavierToLighterNeighbours.h
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for SwapHeavierToLighterNeighbours.
// A SwapHeavierToLighterNeighbours is an object used to balance the partition.
// It will go thought the weigted partition graph and shed load from a heavier
// to its lighter neighbours if the ratio between the two is greater than a 
// certain percent. It repaets this process a number of times.
//
// What: "@(#) SwapHeavierToLighterNeighbours.h, revA"

#ifndef SwapHeavierToLighterNeighbours_h
#define SwapHeavierToLighterNeighbours_h

#include <LoadBalancer.h>

class SwapHeavierToLighterNeighbours: public LoadBalancer
{
  public:
    SwapHeavierToLighterNeighbours();
    SwapHeavierToLighterNeighbours(double factorGreater, 
				   int numReleases);
    
    virtual  ~SwapHeavierToLighterNeighbours();    

    virtual int balance(Graph &theWeightedGraph);

  protected:    
	
  private:
    int numReleases;
    double factorGreater;
    bool disallowDisconnectedGraphs;            
};

#endif


