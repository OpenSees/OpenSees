/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

#ifndef Flume_h
#define Flume_h

#include <Node.h>

#include <vector>

#include "Mesh.h"

class Flume : public Mesh {
   public:
    Flume(int tag, const std::vector<double>& crds,
          const std::vector<double>& dimensions, bool top);
    ~Flume();
    int mesh();

    virtual Node* create_node(const std::vector<double>& crds,
                              int& nodeTag);
    virtual int create_line(Node* nd1, Node* nd2, int& nodeTag,
                            int dir);
    virtual int create_face(Node* nd1, Node* nd2, int& nodeTag,
                            int dir1, int dir2);

    const std::vector<double>& getCrds() const { return crds; }
    const std::vector<double>& getDimensions() const {
        return dimensions;
    }

   private:
    const std::vector<double> crds;
    const std::vector<double> dimensions;
    bool top;
};

#endif
