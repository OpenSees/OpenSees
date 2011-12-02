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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-25 23:34:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/AlgorithmIncrements.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 01/01
// Revision: A
//
// Description: This file contains the class definition for AlgorithmIncrements.
// A AlgorithmIncrements will display the X and B in the SOE associated with the
// algorithm on a record.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef AlgorithmIncrements_h
#define AlgorithmIncrements_h

#include <Recorder.h>

#include <fstream>
using std::ofstream;

class EquiSolnAlgo;
class Renderer;
class ColorMap;
class ID;
class Vector;

class AlgorithmIncrements : public Recorder
{
  public:
    AlgorithmIncrements(EquiSolnAlgo *theAlgo,
			const char *windowTitle, 
			int xLoc, int yLoc, int width, int height,
			bool displayRecord = false,
			const char *fileName = 0);
    
    ~AlgorithmIncrements();    

    int plotData(const Vector &X, const Vector &B);

    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    void restart(void);    

  protected:

  private:
    ColorMap *theMap;
    Renderer *theRenderer;
    EquiSolnAlgo *theAlgo;

    int numRecord;
    bool displayRecord;
    char *fileName;
    ofstream theFile;     
};

#endif







