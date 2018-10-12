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


#ifndef HigherOrder_H
#define HigherOrder_H

#include <vector>
#include <map>
#include <set>

class HigherOrder
{
public:
    
    typedef std::vector<int> VInt;
    typedef std::set<int> SInt;
    typedef std::map<SInt, VInt> HO;
    typedef std::map<SInt, SInt> MidEle;

public:
    HigherOrder();
    ~HigherOrder();

    const VInt& operator()(const VInt& face) const;
    VInt& operator()(const VInt& face);

    void addEle(const VInt& face, int ele);
    bool removeEle(const VInt& face, int ele);

private:
    
    HO ho;
    VInt emp;
    MidEle midele;
    
};

HigherOrder& OPS_getHigherOrder();

#endif
