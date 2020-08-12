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
// Description: For defining mesh
//

#include "Mesh.h"
#include <MapOfTaggedObjects.h>
#include <Node.h>
#include <Element.h>
#include <Domain.h>
#include <elementAPI.h>
#include <Pressure_Constraint.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <vector>

int Mesh::startNodeTag = 1;

#ifdef _PARALLEL_INTERPRETERS

#include <mpi.h>
#endif

void *OPS_ElasticBeam2d(const ID &info);

void *OPS_ForceBeamColumn2d(const ID &info);

void *OPS_DispBeamColumn2d(const ID &info);

void *OPS_PFEMElement2DCompressible(const ID &info);

void *OPS_PFEMElement2DBubble(const ID &info);

void *OPS_PFEMElement3DBubble(const ID &info);

void *OPS_PFEMElement2Dmini(const ID &info);

void *OPS_Tri31(const ID &info);

void *OPS_FourNodeTetrahedron(const ID &info);

void *OPS_ShellMITC4(const ID &info);

void *OPS_CorotTrussElement(const ID &info);

// msh objects
static MapOfTaggedObjects theMeshObjects;

bool OPS_addMesh(Mesh *msh) {
    return theMeshObjects.addComponent(msh);
}

bool OPS_removeMesh(int tag) {
    TaggedObject *obj = theMeshObjects.removeComponent(tag);
    if (obj != 0) {
        delete obj;
        return true;
    }
    return false;
}

Mesh *OPS_getMesh(int tag) {

    TaggedObject *theResult = theMeshObjects.getComponentPtr(tag);
    if (theResult == 0) {
        return 0;
    }
    Mesh *theMsh = (Mesh *) theResult;

    return theMsh;
}

void OPS_clearAllMesh(void) {
    theMeshObjects.clearAll();
}

TaggedObjectIter &OPS_getAllMesh() {
    return theMeshObjects.getComponents();
}

//  mesh class

Mesh::Mesh(int tag, int numelends)
        : TaggedObject(tag), meshsize(0.0),
          id(0), ndf(OPS_GetNDF()), numelenodes(numelends),
          ndtags(), newndtags(), eletags(), elenodes(),
          eleType(0), fluid(false) {
}

Mesh::~Mesh() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) return;

    clearEles();
    clearNodes();
}

int
Mesh::clearEles() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // clear elements
    for (int i = 0; i < eletags.Size(); ++i) {
        Element *ele = domain->removeElement(eletags(i));
        if (ele != 0) {
            delete ele;
        }
    }
    eletags = ID();
    elenodes = ID();

    return 0;
}

int
Mesh::clearNodes() {
    Domain *domain = OPS_GetDomain();

    // clear nodes
    for (int i = 0; i < newndtags.Size(); ++i) {
        Node *node = domain->removeNode(newndtags(i));
        if (node != 0) {
            delete node;
        }
        Pressure_Constraint *pc = domain->removePressure_Constraint(newndtags(i));
        if (pc != 0) {
            delete pc;
        }
    }

    newndtags = ID();

    return 0;
}

void
Mesh::Print(OPS_Stream &s, int flag) {
    s << "tag: " << this->getTag() << "\n";
    s << "mesh size: " << meshsize << "\n";
    s << "id: " << id << "\n";
    s << "ndf: " << ndf << "\n";
}

Node *
Mesh::newNode(int tag, const Vector &crds) {
    // check ndf
    if (ndf <= 0) return 0;

    // craete new node
    Node *node = 0;
    if (crds.Size() == 1) {
        node = new Node(tag, ndf, crds(0));
    } else if (crds.Size() == 2) {
        node = new Node(tag, ndf, crds(0), crds(1));
    } else if (crds.Size() == 3) {
        node = new Node(tag, ndf, crds(0), crds(1), crds(2));
    }

    return node;
}

int
Mesh::nextNodeTag() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    NodeIter &nodes = domain->getNodes();
    Node *node = 0;
    int ndtag = 0;

    while ((node = nodes()) != 0) {
        ndtag = node->getTag();
    }

    if (startNodeTag > ndtag + 1) {
        return startNodeTag;
    }

    return ndtag + 1;
}

int
Mesh::nextEleTag() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    ElementIter &eles = domain->getElements();
    Element *ele = 0;
    int etag = 0;

    while ((ele = eles()) != 0) {
        etag = ele->getTag();
    }

    return etag + 1;
}

int
Mesh::updateMesh() {
    // remove
    this->clearEles();
    this->clearNodes();

    // remesh
    if (this->mesh() < 0) {
        opserr << "WARNING: failed to update mesh\n";
        return -1;
    }

    return 0;
}

int
Mesh::setEleArgs() {
    // no elements
    if (OPS_GetNumRemainingInputArgs() < 1) {
        eleType = 0;
        return 0;
    }

    // get type and set info
    const char *type = OPS_GetString();
    int ndm = OPS_GetNDM();
    ID info(2);
    info(0) = 1; //save data
    info(1) = this->getTag();

    if (strcmp(type, "elasticBeamColumn") == 0) {
        if (ndm == 2) {
            eleType = ELE_TAG_ElasticBeam2d;
            if (OPS_ElasticBeam2d(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            numelenodes = 2;
        }

    } else if (strcmp(type, "forceBeamColumn") == 0) {
        if (ndm == 2) {
            eleType = ELE_TAG_ForceBeamColumn2d;
            if (OPS_ForceBeamColumn2d(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            numelenodes = 2;
        }

    } else if (strcmp(type, "dispBeamColumn") == 0) {
        if (ndm == 2) {
            eleType = ELE_TAG_DispBeamColumn2d;
            if (OPS_DispBeamColumn2d(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            numelenodes = 2;
        }

    } else if (strcmp(type, "PFEMElementBubble") == 0) {
        if (ndm == 2) {
            eleType = ELE_TAG_PFEMElement2DBubble;
            if (OPS_PFEMElement2DBubble(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            fluid = true;
            numelenodes = 3;

        } else {
            eleType = ELE_TAG_PFEMElement3DBubble;
            if (OPS_PFEMElement3DBubble(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            fluid = true;
            numelenodes = 4;
        }

    } else if (strcmp(type, "MINI") == 0) {
        if (ndm == 2) {
            eleType = ELE_TAG_PFEMElement2Dmini;
            if (OPS_PFEMElement2Dmini(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            fluid = true;
            numelenodes = 3;
        }

    } else if (strcmp(type, "PFEMElementCompressible") == 0) {
        opserr << "WARNING: PFEMElementCompressible needs fix in TriMesh\n";
        return -1;
        if (ndm == 2) {
            eleType = ELE_TAG_PFEMElement2DCompressible;
            if (OPS_PFEMElement2DCompressible(info) == 0) {
                opserr << "WARNING: failed to read eleArgs\n";
                return -1;
            }
            fluid = true;
            numelenodes = 3;
        }

    } else if (strcmp(type, "tri31") == 0) {
        eleType = ELE_TAG_Tri31;
        if (OPS_Tri31(info) == 0) {
            opserr << "WARNING: failed to read eleArgs\n";
            return -1;
        }
        numelenodes = 3;

    } else if (strcmp(type, "FourNodeTetrahedron") == 0) {
        eleType = ELE_TAG_FourNodeTetrahedron;
        if (OPS_FourNodeTetrahedron(info) == 0) {
            opserr << "WARNING: failed to read eleArgs\n";
            return -1;
        }
        numelenodes = 4;
    } else if (strcmp(type, "ShellMITC4") == 0) {
        eleType = ELE_TAG_ShellMITC4;
        if (OPS_ShellMITC4(info) == 0) {
            opserr << "WARNING: failed to read eleArgs\n";
            return -1;
        }
        numelenodes = 4;

    } else if (strcmp(type, "corotTruss") == 0) {
        eleType = ELE_TAG_CorotTruss;
        if (OPS_CorotTrussElement(info) == 0) {
            opserr << "WARNING: failed to read eleArgs\n";
            return -1;
        }
        numelenodes = 2;

    } else {
        opserr << "WARNING: element " << type << " is not currently supported in mesh\n";
        return -1;
    }

    return 0;
}

void*
Mesh::getEleArgs()
{
    // function to call
    void *(*OPS_Func)(const ID &info) = 0;
    switch (eleType) {
        case ELE_TAG_PFEMElement2DBubble:
            OPS_Func = OPS_PFEMElement2DBubble;
            break;
        case ELE_TAG_PFEMElement3DBubble:
            OPS_Func = OPS_PFEMElement3DBubble;
            break;
        default:
            return 0;
    }

    ID info(2);
    info(0) = 3; // get arguments
    info(1) = this->getTag(); // mesh tag
    return OPS_Func(info);
}

void
Mesh::addEleTags(const ID &tags) {
    for (int i = 0; i < tags.Size(); ++i) {
        eletags.insert(tags(i));
    }
}

void
Mesh::addEleNodes(const ID &tags) {
    int oldsize = elenodes.Size();
    elenodes.resize(oldsize + tags.Size());
    for (int i = 0; i < tags.Size(); ++i) {
        elenodes(oldsize + i) = tags(i);
    }
}

int
Mesh::newElements(const ID &elends) {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // if no element, no new elements
    if (eleType == 0) return 0;
    if (elends.Size() < numelenodes) {
        return 0;
    }

    // function to call
    void *(*OPS_Func)(const ID &info) = 0;
    switch (eleType) {
        case ELE_TAG_ElasticBeam2d:
            OPS_Func = OPS_ElasticBeam2d;
            break;
        case ELE_TAG_ForceBeamColumn2d:
            OPS_Func = OPS_ForceBeamColumn2d;
            break;
        case ELE_TAG_DispBeamColumn2d:
            OPS_Func = OPS_DispBeamColumn2d;
            break;
        case ELE_TAG_PFEMElement2DBubble:
            OPS_Func = OPS_PFEMElement2DBubble;
            break;
        case ELE_TAG_PFEMElement3DBubble:
            OPS_Func = OPS_PFEMElement3DBubble;
            break;
        case ELE_TAG_PFEMElement2Dmini:
            OPS_Func = OPS_PFEMElement2Dmini;
            break;
        case ELE_TAG_PFEMElement2DCompressible:
            OPS_Func = OPS_PFEMElement2DCompressible;
            break;
        case ELE_TAG_Tri31:
            OPS_Func = OPS_Tri31;
            break;
        case ELE_TAG_FourNodeTetrahedron:
            OPS_Func = OPS_FourNodeTetrahedron;
            break;
        case ELE_TAG_ShellMITC4:
            OPS_Func = OPS_ShellMITC4;
            break;
        case ELE_TAG_CorotTruss:
            OPS_Func = OPS_CorotTrussElement;
            break;
        default:
            OPS_Func = OPS_ElasticBeam2d;
    }


    int eletag = this->nextEleTag();

    // create elements
    ID neweletags(elends.Size() / numelenodes);
    std::vector<Element *> neweles(neweletags.Size());

#pragma omp parallel for
    for (int i = 0; i < neweletags.Size(); ++i) {

        // element tag
        neweletags(i) = eletag + i;

        // info
        ID info(numelenodes + 3);
        info(0) = 2; // load data
        info(1) = this->getTag(); // mesh tag
        info(2) = neweletags(i); // ele tag
        for (int j = 0; j < numelenodes; ++j) {
            // get elenode
            info(3 + j) = elends(numelenodes * i + j);
        }

        // create element
        neweles[i] = (Element *) OPS_Func(info);
    }

    // add elements to domain
    for (unsigned int i = 0; i < neweles.size(); ++i) {
        if (neweles[i] == 0) {
            opserr << "WARING: run out of memory for creating element\n";
            return -1;
        }
        if (domain->addElement(neweles[i]) == false) {
            opserr << "WARNING: failed to add element to domain\n";
            delete neweles[i];
            return -1;
        }
    }

    this->addEleTags(neweletags);

    return 0;
}
