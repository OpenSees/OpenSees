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

class Domain;
class Node;

#include <ID.h>
#include <Vector.h>
#include "ParticleGroup.h"
#include <vector>
#include <PVDRecorder.h>
#include "BackgroundGrid.h"
#include "BackgroundFixData.h"
#include "BackgroundStructure.h"

class BackgroundMesh
{
public:
    BackgroundMesh();
    ~BackgroundMesh();

    // set mesh size
    void setMeshsize(double size) {grids.setSize(size);}
    double getMeshsize() const {return grids.getSize();}

    // add new particle group
    ParticleGroup* newParticleGroup();

    // get all groups
    int numParticleGroups() const {return (int)groups.size();}
    ParticleGroup* getParticleGroup(int i) {return (i>=0&&i<numParticleGroups())? groups[i]:0;return 0;}

    // mapping
    int mapParticleToBack();
    int mapBackToParticle();

    // add fix information
    void addFixInfo(const Vector& min, const Vector& max, const ID& fix);

    // add structural node
    void addStructuralNode(int tag);

    // saving
    void save(const char* filename);

    // set recorder
    void setRecorder(PVDRecorder* recorder);
    int record(int ctag, double timestamp);

private:
    int particlesInGrids();
    int nodesInGrids();
    int structureInGrids();
    int mesh();
    int meshGrid(Node** nodes, GridIndex* index, ID& order);
    int reorder(Node** nodes, ID& order);
    int removeEmptyElements();
    int fix();

    int moveParticles();
    // int structureToGrids();
    void clear();
    
    static double QuinticKernel(double q, double h, int ndm);
    static void getNForRect(double x0, double y0, double hx, double hy, double x, double y,
    			    Vector& N);
    static void getNForTri(double x1, double y1, double x2, double y2,
    			   double x3, double y3, double x, double y,
    			   Vector& N);
    static double pi;
    int findNodeTag();
    int findEleTag();

private:

    std::vector<ParticleGroup*> groups;
    BackgroundFixData fixInfo;
    BackgroundGrid grids;
    ID structuralNodes;
    std::vector<Vector> structuralCoord;
    std::map<int, ID> connectedNodes;
    PVDRecorder* theRecorder;
};

BackgroundMesh& OPS_GetBackgroundMesh();


#endif
