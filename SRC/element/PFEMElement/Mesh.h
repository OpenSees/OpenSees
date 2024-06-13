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
// Description: A piece of mesh
//

#ifndef Mesh_h
#define Mesh_h

#include <TaggedObject.h>
#include <TaggedObjectIter.h>
#include <ID.h>

class Node;
class Element;

class Mesh : public TaggedObject
{
public:

    // constructor
    Mesh(int tag, int numelends);
    virtual ~Mesh();

    // set/get
    virtual void setMeshsize(double size) {meshsize = size;}
    virtual double getMeshsize() const {return meshsize;}

    virtual void setID(int i) {id = i;}
    virtual int getID() const {return id;}

    virtual void setNdf(int n) {ndf = n;}
    virtual int getNdf() const {return ndf;}

    virtual const ID& getNodeTags() const {return ndtags;}
    virtual void setNodeTags(const ID& tags) {ndtags=tags;}

    virtual const ID& getNewNodeTags() const {return newndtags;}
    virtual void setNewNodeTags(const ID& tags) {newndtags=tags;}

    virtual const ID& getEleTags() const {return eletags;}
    virtual void setEleTags(const ID& tags) {eletags=tags;}
    virtual void addEleTags(const ID& tags);

    virtual const ID& getEleNodes() const {return elenodes;}
    virtual void setEleNodes(const ID& tags) {elenodes=tags;}
    virtual void addEleNodes(const ID& tags);

    virtual int getNumEleNodes() const {return numelenodes;}

    virtual int setEleArgs();
    virtual void* getEleArgs();
    virtual bool hasEleArgs() {return eleType != 0;}

    virtual bool isFluid() const {return fluid;}
    virtual int getEleType() const {return eleType;}

    // mesh
    virtual bool ismeshed() const {return newndtags.Size()!=0;}
    virtual int mesh() = 0;
    virtual int updateMesh();
    virtual void Print(OPS_Stream &s, int flag =0);

    // clear mesh
    virtual int clearEles();
    virtual int clearNodes();

    // create new element
    virtual int newElements(const ID& elenodes);
    virtual Node* newNode(int tag, const Vector& crds);

    // find the next available tag
    static int nextNodeTag();
    static int nextEleTag();
    static void setStartNodeTag(int tag) {startNodeTag = tag;}

private:

    double meshsize;
    int id;
    int ndf;
    int numelenodes;

    ID ndtags, newndtags, eletags, elenodes;
    int eleType;
    bool fluid;

    static int startNodeTag;
};

bool OPS_addMesh(Mesh* msh);
bool OPS_removeMesh(int tag);
Mesh *OPS_getMesh(int tag);
void OPS_clearAllMesh(void);
TaggedObjectIter& OPS_getAllMesh();

#endif
