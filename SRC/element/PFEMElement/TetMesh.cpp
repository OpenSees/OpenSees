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
// $Date $

// Written: Minjie Zhu
//
// Description: The class TetMesh is for meshing triangles
//

#include "TetMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <string.h>
#include "TetMeshGenerator.h"
#include <elementAPI.h>
#include <Element.h>
#include <set>
#include <map>
#include <vector>
#include <Pressure_Constraint.h>

int OPS_TetMesh()
{
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr<<"WARNING: want tag? nummesh? mtags? id? ndf? size? eleType? eleArgs?\n";
        return -1;
    }

    // get tag and number mesh
    int num = 2;
    int idata[2];
    if (OPS_GetIntInput(&num,idata) < 0) {
        opserr<<"WARNING: failed to read mesh tag and number of 2D boundary mesh\n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < idata[1]+3) {
        opserr<<"WARNING: want mtags? id? ndf? size? <eleType? eleArgs?>\n";
        return -1;
    }

    // create mesh
    TetMesh* mesh = new TetMesh(idata[0]);
    if(OPS_addMesh(mesh) == false) {
        opserr<<"WARNING: failed to add mesh\n";
        return -1;
    }

    // get mesh tags
    num = idata[1];
    ID mtags(num);
    if (OPS_GetIntInput(&num,&mtags(0)) < 0) {
        opserr<<"WARNING: failed to read boundary mesh tags\n";
        return -1;
    }
    mesh->setMeshTags(mtags);

    // get id and ndf
    num = 2;
    int data[2];
    if (OPS_GetIntInput(&num,data) < 0) {
        opserr<<"WARNING: failed to read id and ndf\n";
        return -1;
    }
    mesh->setID(data[0]);
    mesh->setNdf(data[1]);

    // get size
    num = 1;
    double size;
    if (OPS_GetDoubleInput(&num,&size) < 0) {
        opserr<<"WARNING: failed to read mesh size\n";
        return -1;
    }
    mesh->setMeshsize(size);

    // set eleArgs
    if (mesh->setEleArgs() < 0) {
        opserr << "WARNING: failed to set element arguments\n";
        return -1;
    }

    // mesh
    if (mesh->mesh() < 0) {
        opserr<<"WARNING: failed to do triangular mesh\n";
        return -1;
    }

    return 0;
}

TetMesh::TetMesh(int tag)
    :Mesh(tag,4), mtags()
{
}

TetMesh::~TetMesh()
{
}

int
TetMesh::mesh()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // check
    double size = this->getMeshsize();
    if(size <= 0) {
        opserr<<"WARNING: mesh size <= 0\n";
        return -1;
    }
    if (mtags.Size() == 0) return 0;

    // get nodes and elements from boundary mesh
    ID ndtags;
    std::vector<TetMeshGenerator::Facet> facets;
    for (int i=0; i<mtags.Size(); ++i) {
	// get mesh
        Mesh* mesh = OPS_getMesh(mtags(i));
        if (mesh == 0) {
            opserr << "WARNING: mesh "<<mtags(i)<<" does not exist\n";
            return -1;
        }

	// get nodes
        const ID& tags = mesh->getNodeTags();
        for (int j=0; j<tags.Size(); ++j) {
            ndtags.insert(tags(j));
        }
        const ID& newtags = mesh->getNewNodeTags();
        for (int j=0; j<newtags.Size(); ++j) {
            ndtags.insert(newtags(j));
        }

	// get polygons
	int numelenodes = mesh->getNumEleNodes();
        const ID& elends = mesh->getEleNodes();
	TetMeshGenerator::Facet facet;
        for (int j=0; j<elends.Size()/numelenodes; ++j) {
	    TetMeshGenerator::Polygon polygon(numelenodes);
	    for (int k=0; k<numelenodes; ++k) {
		polygon[k] = elends(numelenodes*j+k);
	    }
            facet.push_back(polygon);
        }

	// add the mesh as a facet
	facets.push_back(facet);
    }
    if(ndtags.Size() < 4) {
        opserr<<"WARNING: input number of nodes < 4\n";
        return -1;
    }
    if(facets.size() < 4) {
        opserr<<"WARNING: input number of facets < 4\n";
        return -1;
    }
    this->setNodeTags(ndtags);

    // get correct index for factes
    for (unsigned int i=0; i<facets.size(); ++i) {
	for (unsigned int j=0; j<facets[i].size(); ++j) {
	    for (unsigned int k=0; k<facets[i][j].size(); ++k) {
		facets[i][j][k] = ndtags.getLocationOrdered(facets[i][j][k]);
		if (facets[i][j][k] < 0) {
		    opserr << "WARNING: failed to get a node in a facet -- Tetmesh::mesh\n";
		    return -1;
		}
	    }
	}
    }

    // calling mesh generator
    TetMeshGenerator gen;
    int nodecounter = nextNodeTag();
    for(int i=0; i<ndtags.Size(); i++) {
        // get node
        Node* theNode = domain->getNode(ndtags(i));
        if(theNode == 0) {
            opserr<<"WARNING: node "<<ndtags(i)<<" is not defined\n";
            return -1;
        }
        Vector crds = theNode->getCrds();
        const Vector& disp = theNode->getTrialDisp();
        if(crds.Size() != 3) {
            opserr<<"WARNING: ndm != 3\n";
            return -1;
        }
        if (disp.Size() >= crds.Size()) {
            for (int j=0; j<crds.Size(); ++j) {
                crds(j) += disp(j);
            }
        }

        // add point
        gen.addPoint(crds(0), crds(1), crds(2), 0);

        // add pc
        if (this->isFluid()) {
	    // create pressure constraint
	    Pressure_Constraint* thePC = domain->getPressure_Constraint(ndtags(i));
	    if(thePC != 0) {
		thePC->setDomain(domain);
	    } else {

		// create pressure node
		Node* pnode = 0;
		pnode = new Node(nodecounter++, 1, crds[0], crds[1], crds[2]);
		if (pnode == 0) {
		    opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		    return -1;
		}
		if (domain->addNode(pnode) == false) {
		    opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
		    delete pnode;
		    return -1;
		}

		thePC = new Pressure_Constraint(ndtags(i), pnode->getTag());
		if(thePC == 0) {
		    opserr<<"WARNING: no enough memory for Pressure_Constraint\n";
		    return -1;
		}
		if (domain->addPressure_Constraint(thePC) == false) {
		    opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
		    delete thePC;
		    return -1;
		}
	    }
        }
    }

    for (unsigned int i=0; i<facets.size(); ++i) {
	gen.addFacet(facets[i],0);
    }

    // meshing
    gen.mesh(size*size*size/(6.0*1.414),false);

    // get points and create nodes
    int nump = gen.getNumPoints();
    if (nump == 0) {
        opserr << "WARNING: no nodes is meshed\n";
        return -1;
    }
    ID newndtags(nump-ndtags.Size());
    ID allndtags(nump);
    for (int i=0; i<ndtags.Size(); ++i) {
        allndtags(i) = ndtags(i);
    }

    for (int i=ndtags.Size(); i<nump; ++i) {
        // get point
        Vector crds(3);
        int mark = 0;
        gen.getPoint(i,crds(0),crds(1),crds(2),mark);
        Node* node = newNode(nodecounter++,crds);
        if (node == 0) {
            opserr << "WARING: failed to create node\n";
            return -1;
        }
        if (domain->addNode(node) == false) {
            opserr << "WARNING: failed to add node to domain\n";
            delete node;
            return -1;
        }
        allndtags(i) = node->getTag();
        newndtags(i-ndtags.Size()) = node->getTag();

        // add pc
        if (this->isFluid()) {
	    // create pressure constraint
	    Pressure_Constraint* thePC = domain->getPressure_Constraint(node->getTag());
	    if(thePC != 0) {
		thePC->setDomain(domain);
	    } else {

		// create pressure node
		Node* pnode = 0;
		pnode = new Node(nodecounter++, 1, crds(0), crds(1), crds(2));
		if (pnode == 0) {
		    opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		    return -1;
		}
		if (domain->addNode(pnode) == false) {
		    opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
		    delete pnode;
		    return -1;
		}

		thePC = new Pressure_Constraint(node->getTag(), pnode->getTag());
		if(thePC == 0) {
		    opserr<<"WARNING: no enough memory for Pressure_Constraint\n";
		    return -1;
		}
		if (domain->addPressure_Constraint(thePC) == false) {
		    opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
		    delete thePC;
		    return -1;
		}
	    }
        }
    }
    this->setNewNodeTags(newndtags);

    // get tetrahedrons
    int numtet = gen.getNumTets();
    if (numtet == 0) return 0;
    ID elenodes(numtet*4);

    for(int i=0; i<numtet; i++) {
        int p1,p2,p3,p4;
        gen.getTet(i,p1,p2,p3,p4);
        elenodes(4*i) = allndtags(p1);
        elenodes(4*i+1) = allndtags(p2);
        elenodes(4*i+2) = allndtags(p3);
        elenodes(4*i+3) = allndtags(p4);
    }
    this->setEleNodes(elenodes);

    // create elemnts
    if (this->newElements(elenodes) < 0) {
        opserr << "WARNING: failed to create elements\n";
        return -1;
    }

    return 0;
}

int
TetMesh::remesh(double alpha)
{
    int ndm = OPS_GetNDM();
    if (ndm != 3) {
	opserr << "WARNING: TetMesh::remesh is only for 3D problem\n";
	return -1;
    }
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // get all nodes
    std::map<int, std::vector<int> > ndinfo;
    TaggedObjectIter& meshes = OPS_getAllMesh();
    ID nodetags;
    Mesh* msh = 0;
    while((msh = dynamic_cast<Mesh*>(meshes())) != 0) {

        int id = msh->getID();
        int mtag = msh->getTag();

        if (id == 0) continue;

        // get bound nodes
        const ID& tags = msh->getNodeTags();
        for (int i=0; i<tags.Size(); ++i) {
            std::vector<int>& info = ndinfo[tags(i)];
            info.push_back(mtag);
            info.push_back(id);
            nodetags.insert(tags(i));
        }

        // get internal nodes
        const ID& newtags = msh->getNewNodeTags();
        for (int i=0; i<newtags.Size(); ++i) {
            std::vector<int>& info = ndinfo[newtags(i)];
            info.push_back(mtag);
            info.push_back(id);
            nodetags.insert(newtags(i));
        }
    }

    if (nodetags.Size() == 0) return 0;

    // calling mesh generator
    TetMeshGenerator gen;
    int nodecounter = nextNodeTag();
    for(int i=0; i<nodetags.Size(); ++i) {

        // get node
        Node* theNode = domain->getNode(nodetags(i));
        if(theNode == 0) {
            opserr<<"WARNING: node "<<nodetags(i)<<" is not defined\n";
            return -1;
        }
        const Vector& crds = theNode->getCrds();
        const Vector& disp = theNode->getTrialDisp();
        if(crds.Size() < ndm || disp.Size() < ndm) {
            opserr<<"WARNING: ndm < 3 or ndf < 3\n";
            return -1;
        }

	// create pc if not
	Pressure_Constraint* thePC = domain->getPressure_Constraint(nodetags(i));
	if(thePC != 0) {
	    thePC->setDomain(domain);
	} else {

	    // create pressure node
	    Node* pnode = 0;
	    pnode = new Node(nodecounter++, 1, crds(0), crds(1), crds(2));
	    if (pnode == 0) {
		opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		return -1;
	    }
	    if (domain->addNode(pnode) == false) {
		opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
		delete pnode;
		return -1;
	    }

	    thePC = new Pressure_Constraint(nodetags(i), pnode->getTag());
	    if(thePC == 0) {
		opserr<<"WARNING: no enough memory for Pressure_Constraint\n";
		return -1;
	    }
	    if (domain->addPressure_Constraint(thePC) == false) {
		opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
		delete thePC;
		return -1;
	    }
	}

        // add point
        gen.addPoint(crds(0)+disp(0), crds(1)+disp(1), crds(2)+disp(2), 0);
    }

    // meshing
    gen.remesh(alpha);

    // get elenodes
    std::map<int,ID> meshelenodes;
    for(int i=0; i<gen.getNumTets(); i++) {

        // get points
        int p1,p2,p3,p4;
        gen.getTet(i,p1,p2,p3,p4);

        // get nodes
        int nds[4];
        nds[0] = nodetags(p1);
        nds[1] = nodetags(p2);
        nds[2] = nodetags(p3);
	nds[3] = nodetags(p4);

        // check if all nodes in same mesh
        std::vector<int>& info1 = ndinfo[nds[0]];
        int mtag = 0, id = 0;
        bool same = false;
        for (int k=0; k<(int)info1.size()/2; ++k) {
            // check if any mesh of node 1 is same for another three nodes
            mtag = info1[2*k];
            id = info1[2*k+1];

            int num = 0;
            for (int j=1; j<4; ++j) {
                std::vector<int>& infoj = ndinfo[nds[j]];
                for (int kj=0; kj<(int)infoj.size()/2; ++kj) {
                    int mtagj = infoj[2*kj];
                    if (mtag == mtagj) {
                        ++num;
                        break;
                    }
                }

            }
            if (num == 3) {
                same = true;
                break;
            }
        }

        // nodes in different mesh
        if (!same) {
            mtag = 0;
            id = 0;
            for (int j=0; j<4; ++j) {
                std::vector<int>& info = ndinfo[nds[j]];
                for (int k=0; k<(int)info.size()/2; ++k) {
                    if (info[2*k+1] < id) {
                        if (dynamic_cast<TetMesh*>(OPS_getMesh(info[2*k])) != 0) {
                            mtag = info[2*k];
                            id = info[2*k+1];
                        }
                    }
                }
            }
        }

        // if all connected to structure
        if (id >= 0) continue;

        // add elenodes to its mesh
        ID& elenodes = meshelenodes[mtag];
        for (int j=0; j<4; ++j) {
            elenodes[elenodes.Size()] = nds[j];
        }
    }

    // creat elements
    for (std::map<int,ID>::iterator it=meshelenodes.begin();
	 it!=meshelenodes.end(); ++it) {

        int mtag = it->first;
        ID& elenodes = it->second;

        msh = OPS_getMesh(mtag);
        if (msh != 0) {

            int id = msh->getID();

            // remove mesh for id<0
            if (id < 0) {
                if (msh->clearEles() < 0) {
                    opserr << "WARNING: failed to clear element in mesh"<<mtag<<"\n";
                    return -1;
                }

                if (msh->newElements(elenodes) < 0) {
                    opserr << "WARNING: failed to create new elements in mesh"<<mtag<<"\n";
                    return -1;
                }
            }
        }
    }

    return 0;
}
