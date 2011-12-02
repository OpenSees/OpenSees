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
// $Date: 2004-01-29 23:30:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/FilePlotter.h,v $
                                                                        
                                                                        
// File: ~/recorder/tcl/FilePlotter.h.h
// 
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class definition for FilePlotter.
// A FilePlotter will create a line graph using xy points found in a file.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef FilePlotter_h
#define FilePlotter_h

#include <Recorder.h>

class Renderer;
class ColorMap;
class ID;

class FilePlotter : public Recorder
{
  public:
    FilePlotter(const char *fileName,
		const char *windowTitle, 
		int xLoc, int yLoc, int width, int height, double dT);
    
    ~FilePlotter();    

    int plotFile();

    int record(int commitTag, double timeStamp);
    int playback(int commitTag);
    void restart(void);    

    int setFile(char *newFile);
    int setCol(const ID &theCols);

  protected:

  private:
    ColorMap *theMap;
    Renderer *theRenderer;
    ID *cols;
    char *fileName;
    
    double deltaT;
    double nextTimeStampToRecord;    
};

#endif







