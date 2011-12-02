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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/TclVideoPlayer.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/tcl/TclVideoPlayer.h.h
// 
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class definition for TclVideoPlayer.
// A TclVideoPlayer can be used to playback movies.
//
// What: "@(#) ModelBuilder.h, revA"

#ifndef TclVideoPlayer_h
#define TclVideoPlayer_h

#include <G3Globals.h>

class Renderer;
class ColorMap;

extern "C" {
#include <tcl.h>
#include <tk.h>
}

class TclVideoPlayer
{
  public:
    TclVideoPlayer(char *title, char *fileName, char *imageName,
				   Tcl_Interp *interp);
    ~TclVideoPlayer();    

    int play(void);

  protected:

  private:
    ColorMap *theMap;
    Renderer *theRenderer;

    char *theFileName;   // file name  
    ifstream theFile; 	 // output stream    
};

#endif







