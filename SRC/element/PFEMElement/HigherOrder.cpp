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

// $Revision $
// $Date$
// $URL$

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: Track the edge nodes in higher order elements
// 
// 

#include "HigherOrder.h"

static HigherOrder high;

HigherOrder& OPS_getHigherOrder()
{
    return high;
}



HigherOrder::HigherOrder()
    :ho(), emp(), midele()
{
}


HigherOrder::~HigherOrder()
{
}

const HigherOrder::VInt& 
HigherOrder::operator()(const VInt& face) const
{
    SInt fs(face.begin(), face.end());
    HO::const_iterator it = ho.find(fs);
    if (it == ho.end()) return emp;

    return it->second;
}

HigherOrder::VInt& 
HigherOrder::operator()(const VInt& face)
{
    SInt fs(face.begin(), face.end());
    return ho[fs];
}

void
HigherOrder::addEle(const VInt& face, int ele)
{
    SInt fs(face.begin(), face.end());
    midele[fs].insert(ele); 
}

bool
HigherOrder::removeEle(const VInt& face, int ele)
{
    SInt fs(face.begin(), face.end());
    SInt& eles = midele[fs];
    eles.erase(ele);

    if (eles.empty()) {
	// remove the face
	midele.erase(fs);
	ho.erase(fs);
    }
    
    return eles.empty();
}
