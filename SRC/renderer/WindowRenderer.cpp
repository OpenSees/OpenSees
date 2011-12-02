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
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/WindowRenderer.cpp,v $
                                                                        
                                                                        
// File: ~/renderer/WindowRenderer.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for WindowRenderer.
// WindowRenderer is an class which diplays using X11 or openGL.
//
// What: "@(#) WindowRenderer.h, revA"

#include <WindowRenderer.h>
#include <ColorMap.h>

#include <db.H>

#include <Vector.h>
#include <View.h>
#include <Projection.h>
#include <Viewport.h>
#include <Clipping.h>
#include <Device.h>
#include <Scan.h>

WindowRenderer::WindowRenderer(int width, int height,
			       ColorMap &_theMap)
  :Renderer(_theMap)
{

}

WindowRenderer::~WindowRenderer()
{

}

