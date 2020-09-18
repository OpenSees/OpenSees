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

#ifndef BackgroundMesh_h
#define BackgroundMesh_h

#include <ID.h>
#include <Node.h>
#include <Recorder.h>

#include <fstream>
#include <set>
#include <vector>

#include "BCell.h"
#include "BNode.h"
#include "BackgroundDef.h"

class BackgroundMesh {
   public:
    BackgroundMesh();
    virtual ~BackgroundMesh();

    // background info
    void setTol(double t) { this->tol = t; }
    double getTol() const { return this->tol; }
    void setRange(const VDouble& l, const VDouble& u);
    void setBasicSize(double size) { bsize = size; }
    void addRecorder(Recorder* recorder);
    int record(bool init = false);
    void setLocs(const VDouble& l) { this->locs = l; }
    int setFile(const char* name);
    void setNumSub(int num) { numsub = num; }
    void addStructuralNodes(VInt& snodes, int sid);
    void setContactData(const VDouble& data);
    void addLargeSize(int numbasic, const VDouble& range_low,
                      const VDouble& range_up);
    bool isDispOn() const { return dispon; }
    void setAlphaS(int sid, double alpha) { alphaS[sid] = alpha; }
    void setDispOn(bool on);

    // remesh all
    int remesh(bool init = false);

    // get grids
    void getIndex(const VDouble& crds, double incr, VInt& index) const;
    void lowerIndex(const VDouble& crds, VInt& index) const;
    void upperIndex(const VDouble& crds, VInt& index) const;
    void nearIndex(const VDouble& crds, VInt& index) const;
    void getCrds(const VInt& index, VDouble& crds) const;
    void getCorners(const VInt& index, int num, VVInt& indices) const;

    // particles
    int addParticles();
    void gatherParticles(const VInt& minindex, const VInt& maxindex,
                         VParticle& pts, bool checkfsi = false);
    int moveParticles();
    int convectParticle(Particle* pt, VInt index, int nums);
    int moveFixedParticles();

    // create grid nodes and elements
    int addStructure();
    int gridNodes();
    int gridFluid();
    int gridFSI();
    int gridFSInoDT();
    int gridEles();
    static int createContact(const VInt& ndtags, const VInt& sids,
                             VInt& elends);

    // particle kernel
    static double QuinticKernel(double q, double h, int ndm);
    static int preNForTri(double x1, double y1, double x2, double y2,
                          double x3, double y3, VDouble& coeff);
    void getNForTri(const VDouble& coeff, double x, double y, VDouble& N);
    static void getNForRect(double x0, double y0, double hx, double hy,
                            double x, double y, VDouble& N);

    static int preNForTet(const VDouble& crds1, const VDouble& crds2,
                          const VDouble& crds3, const VDouble& crds4,
                          VVDouble& coeff);
    void getNForTet(const VVDouble& coeff, const VDouble& crds,
                    VDouble& N);
    static void getNForRect(double x0, double y0, double z0, double hx,
                            double hy, double hz, double x, double y,
                            double z, VDouble& N);
    static bool check_area(const VDouble& ndcrds1, const VDouble& ndcrds2,
                           const VDouble& ndcrds3);
    static bool check_vol(const VDouble& ndcrds1, const VDouble& ndcrds2,
                          const VDouble& ndcrds3, const VDouble& ndcrds4);
    //    1
    //    /\ upper
    // 0 /  \ 2
    //  ------
    //  | 4  |
    //  | /\ | lower
    //  |/  \|
    // 3------5
    static void splitPrism(const VInt& prism, VVInt& tets, bool incl1,
                           bool incl4);

    // clear all
    void clearAll();
    int clearBackground();
    static void clearGridEles();
    void clearGrid();

    // interpolate in a cell
    int interpolate(Particle* pt, const VVInt& index, const VVDouble& vels,
                    const VVDouble& dvns, const VDouble& pns,
                    const VDouble& dpns, const VVDouble& crds,
                    const std::vector<BackgroundType>& types,
                    const VVInt& ndtags, double dt);
    static int interpolate(const VVDouble& values, const VDouble& N,
                           VDouble& newvalue);
    static int interpolate(const VDouble& values, const VDouble& N,
                           double& newvalue);
    static int solveLine(const VDouble& p1, const VDouble& dir, int dim,
                         double crd, double& k);
    static bool inEle(const VDouble& N);

   private:
    VInt lower, upper;
    std::map<VInt, BCell> bcells;
    std::map<VInt, BNode> bnodes;
    double tol;
    double bsize;
    int numave, numsub;
    std::vector<Recorder*> recorders;
    VDouble locs;
    double currentTime;
    std::ofstream theFile;
    std::map<int, VInt> structuralNodes;  // >0:structure, <0: fluid,
                                          // 0:invalid, larger:debris
    VDouble contactData;
    VInt contactEles;
    bool dispon;
    std::map<int, double> alphaS;  // alphaS for sids

    static const int contact_tag = -13746;
};

BackgroundMesh& OPS_getBgMesh();

#endif
