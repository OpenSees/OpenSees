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

// $Revision$
// $Date$
// $URL$

// Written: Minjie Zhu
//
// Description: This class defines the BackgroundMesh
//

#include "BackgroundMesh.h"

static BackgroundMesh bgmesh;

BackgroundMesh& OPS_getBgMesh()
{
    return bgmesh;
}

// OPS_BgMesh
int OPS_BgMesh()
{
    return 0;
}

BackgroundMesh::BackgroundMesh()
{
}

BackgroundMesh::~BackgroundMesh()
{
}

void
BackgroundMesh::addRecorder(Recorder*)
{
}

int
BackgroundMesh::remesh()
{
    return 0;
}
