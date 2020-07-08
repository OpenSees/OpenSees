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

#include "BackgroundDef.h"
#include <vector>
#include <set>
#include <ID.h>
#include <Node.h>
#include <Recorder.h>
#include <fstream>

class BackgroundMesh
{
private:

    static int FLUID, STRUCTURE, FIXED;

    // FLUID - a grid fluid node
    // STRUCTURE - a structural node
    // FIXED - a fixed grid fluid node
    struct BNode {

        VInt tags;
        VVDouble crdsn;
        VVDouble vn;
        VVDouble incrv;
        VVDouble dvn;
        VDouble pn;
        VDouble dpn;
        int type;
        VInt sid; // structure id, <0:fluid, >0:structure, =0:not in contact

        BNode():tags(),crdsn(),vn(),incrv(), dvn(),pn(),dpn(),type(FLUID){}
        void addNode(int tag, const VDouble& crds, const VDouble& v,
                     const VDouble& dv, double p, double dp, int tp, int id=-1) {
            tags.push_back(tag);
            crdsn.push_back(crds);
            vn.push_back(v);
            incrv.push_back(v);
            dvn.push_back(dv);
            pn.push_back(p);
            dpn.push_back(dp);
            type = tp;
            sid.push_back(id);
        }
        void clear() {
            tags.clear();
            crdsn.clear();
            vn.clear();
            incrv.clear();
            dvn.clear();
            pn.clear();
            dpn.clear();
            type = FLUID;
            sid.clear();
        }
        int size() const {return (int)tags.size();}

    };

    // FLUID - a grid fluid cell
    // STRUCTURE - a structural cell, which should have no particles
    struct BCell {
        VParticle pts;
        int type;
        std::vector<BNode*> bnodes;
        std::vector<VInt> bindex;
        int sizeLevel;

        BCell():pts(),type(FLUID),bnodes(),bindex(),sizeLevel(0) {}
        void add(Particle* pt) {pts.push_back(pt);}
    };

public:
    BackgroundMesh();
    virtual ~BackgroundMesh();

    // background info
    void setTol(double t) {this->tol = t;}
    double getTol() const {return this->tol;}
    void setMeshTol(double t) {meshtol = t;}
    double getMeshTol() const {return meshtol;}
    void setRange(const VDouble& l, const VDouble& u);
    void setBasicSize(double size) {bsize = size; dispon= false;}
    void addRecorder(Recorder* recorder);
    int record(bool init=false);
    void setLocs(const VDouble& l) {this->locs = l;}
    int setFile(const char* name);
    void setNumSub(int num) {numsub = num;}
    void addStructuralNodes(VInt& snodes, int sid);
    void setContactData(const VDouble& data);
    bool isIncrVel() const { return incrVel;}
    void setIncrVel(bool ivel) {incrVel = ivel;}
    bool isFSITri() const {return fsiTri;}
    void setFSITri(bool fsi) {fsiTri = fsi;}
    void setBoundReduceFactor(double factor) {boundReduceFactor = factor;}
    void addLargeSize(int numbasic,
                      const VDouble& range_low,
                      const VDouble& range_up);
    int getSizeLevel(VInt& index);
    void setPressureOnce(bool flag) {pressureonce = flag;}
    bool isPressureOnce() const {return pressureonce;}
    bool isDispOn() const {return dispon;}
    void setFastAssembly(bool flag) {fastAssembly = flag;}
    bool isFastAssembly() const {return fastAssembly;}
    void setKernelClose(bool flag) {kernelClose = flag;}
    bool isKernelClose() const {return kernelClose;}

    // remesh all
    int remesh(bool init=false);

    // get grids
    void getIndex(const VDouble& crds, double incr, VInt& index) const;
    void lowerIndex(const VDouble& crds, VInt& index) const;
    void upperIndex(const VDouble& crds, VInt& index) const;
    void nearIndex(const VDouble& crds, VInt& index) const;
    void getCrds(const VInt& index, VDouble& crds) const;
    void getCorners(const VInt& index, int num, int level, VVInt& indices) const;

    // particles
    int addParticles();
    void gatherParticles(const VInt& minindex, const VInt& maxindex,
                         VParticle& pts, bool checkfsi=false);
    int moveParticles();
    int convectParticle(Particle* pt, VInt index, int level, int nums);
    int moveFixedParticles();
    int inlet();
    void addInlet(const VDouble& crds, const VDouble& vel);
    void setInletNum(const VInt& num) {inletNum = num;}

    // create grid nodes and elements
    int addStructure();
    int gridNodes();
    int gridFluid();
    int gridFSI(ID& freenodes);
    int gridEles();
    static int createContact(const VInt& ndtags, const VInt& sids, VInt& elends);

    // particle kernel
    static double QuinticKernel(double q, double h, int ndm);
    static int preNForTri(double x1, double y1, double x2, double y2,
                          double x3, double y3, VDouble& coeff);
    void getNForTri(const VDouble& coeff, double x, double y,
                    VDouble& N);
    static void getNForRect(double x0, double y0, double hx, double hy,
                            double x, double y, VDouble& N);

    static int preNForTet(const VDouble& crds1, const VDouble& crds2,
                          const VDouble& crds3, const VDouble& crds4,
                          VVDouble& coeff);
    void getNForTet(const VVDouble& coeff, const VDouble& crds,
                    VDouble& N);
    static void getNForRect(double x0, double y0, double z0,
                            double hx, double hy, double hz,
                            double x, double y, double z,
                            VDouble& N);


    // clear all
    void clearAll();
    int clearBackground();
    static void clearGridEles();
    void clearGrid();

    // interpolate in a cell
    int interpolate(Particle* pt, const VVInt& index,
                    const VVDouble& vels, const VVDouble& incrvels,
                    const VVDouble& dvns,
                    const VDouble& pns, const VDouble& dpns,
                    const VVDouble& crds,
                    const VInt& fixed, const VVInt& ndtags,
                    double dt);
    static int interpolate(const VVDouble& values, const VDouble& N, VDouble& newvalue);
    static int interpolate(const VDouble& values, const VDouble& N, double& newvalue);
    static int solveLine(const VDouble& p1, const VDouble& dir,
                         int dim, double crd, double& k);
    static bool inEle(const VDouble& N);

    // find free surface
    int findFreeSurface(const ID& freenodes);
    void setFreeSurface() {freesurface = true;}

    // wall
    void getWall(VDouble& dir, double& dist, const VDouble& xbnd,
                 const VDouble& ybnd, const VDouble& zbnd,
                 const VDouble pcrds);

private:

    VInt lower, upper;
    std::map<VInt,BCell> bcells;
    std::map<VInt,BNode> bnodes;
    double tol, meshtol;
    double bsize;
    int numave, numsub;
    std::vector<Recorder*> recorders;
    VDouble locs;
    double currentTime;
    std::ofstream theFile;
    std::map<int, VInt> structuralNodes; // >0:structure, <0: fluid, 0:invalid, larger:debris
    bool freesurface;
    VDouble contactData;
    VInt contactEles;
    bool incrVel;
    bool fsiTri; // move partiles in fsi area through triangles
    double boundReduceFactor;
    VVInt inletLoc;
    VVDouble inletVel;
    VInt inletNum;
    VVInt largesize;
    bool pressureonce;
    bool dispon;
    bool fastAssembly;
    bool kernelClose;

    static const int contact_tag = -13746;
};

BackgroundMesh& OPS_getBgMesh();

#endif
