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
// $Date: 02/10/06 00:06:00 $
// $Source: /usr/local/cvs/OpenSees/OTHER/DRM/DRMLoadPattern.h $

/* 
* @author: gnp <petropoulos@gmail.com>
 *
 * @Description: DRM Load Pattern Header File
 *
 * @Date: 2/10/06
 */


#ifndef DRMLoadPattern_h
#define DRMLoadPattern_h

#include <LoadPattern.h>
#include <map>
#include <set>
#include <fstream>

class DRMInputHandler;
class Domain;

class DRMLoadPattern : public LoadPattern
{
  public :
    
    DRMLoadPattern(int tag, double cfact, DRMInputHandler* my_handler, Domain* domain);
  
    virtual ~DRMLoadPattern();

    void applyLoad(double time);

    void setMaps();

 private:
   Domain *myDomain;
   DRMInputHandler* myHandler;

   double factor;
   std::map<int,int> eNodes;
   std::map<int,Element*> elem;
   
   std::map<int,Vector*> storage;
   std::map<int,int> storage2;

   std::ofstream ofile;
};
#endif
