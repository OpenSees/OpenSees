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

// $Revision: 1.0 $
// $Date: 2016-1-27  $
                                                                        
// Written: Minjie Zhu
//
// Description: This class defines the BackgroundMesh 
//

#ifndef BackgroundFixData_h
#define BackgroundFixData_h

#include <vector>
#include <Domain.h>
#include <ID.h>

class BackgroundFixData
{
public:
    BackgroundFixData();
    ~BackgroundFixData();

    void addInfo(const Vector& min, const Vector& max, const ID& fix);

    int tryFix(int ndtag, Domain& domain);
    void clear(Domain& domain);

private:
    std::vector<Vector> min;
    std::vector<Vector> max;
    std::vector<ID> fix;
    std::vector<SP_Constraint*> sps;
};

#endif
