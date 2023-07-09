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

#ifndef BackgroundGrid_h
#define BackgroundGrid_h

#include <vector>
#include <map>
#include <ID.h>

class Particle;
class Node;
class Element;
class Domain;

class GridIndex
{
public:
    GridIndex():i(-1), j(-1), valid(false) {}
    GridIndex(int ii, int jj):i(ii),j(jj), valid(true) {}
    ~GridIndex(){}

    bool operator<(const GridIndex& index) const {
	return std::make_pair(i,j) < std::make_pair(index.i,index.j);
	// if (i == index.i) return j < index.j;
	// return i < index.i;
    }

    bool isValid() const {return valid;}

    GridIndex north() const {return GridIndex(i,j+1);}
    GridIndex south() const {return GridIndex(i,j-1);}
    GridIndex east() const {return GridIndex(i+1,j);}
    GridIndex west() const {return GridIndex(i-1,j);}
    GridIndex northEast() const {return GridIndex(i+1,j+1);}
    GridIndex northWest() const {return GridIndex(i-1,j+1);}
    GridIndex southEast() const {return GridIndex(i+1,j-1);}
    GridIndex southWest() const {return GridIndex(i-1,j-1);}

    double getX(double size) const {return i*size;}
    double getY(double size) const {return j*size;}

private:
    int i,j;
    bool valid;
};

class BackgroundGrid
{
    struct GridData {
	GridData():particles(),node(0),elements() {}
	
	std::vector<Particle*> particles;
	Node* node;
	std::vector<Element*> elements;
    };
    
public:
    
    BackgroundGrid();
    ~BackgroundGrid();

    // clear all
    void clear(const ID& structuralNodes);

    // size
    void setSize(double sz) {size = sz;}
    double getSize() const {return size;}

    // add and get
    void addGrid(const GridIndex& index);
    void addParticle(const GridIndex& index, Particle* p);
    void setNode(const GridIndex& index, Node* nd);
    void addElement(const GridIndex& index, Element* e);
    
    std::vector<Particle*>* getParticles(const GridIndex& index);
    Node* getNode(const GridIndex& index);
    std::vector<Element*>* getElements(const GridIndex& index);
    bool hasGrid(const GridIndex& index);

    int numGridPoints() const {return (int)data.size();}

    // iterators
    GridIndex getIndex() const;
    std::vector<Particle*>* getParticles();
    std::vector<Element*>* getElements();
    Node* getNode();
    void addParticle(Particle* p);
    void addElement(Element* e);
    void setNode(Node* nd);
    void reset() {iter = data.begin();}
    void reset(const GridIndex& index);
    bool isEnd() const {return iter == data.end();}
    void next();

    // check if corner
    bool isCorner(const GridIndex& center) const;

private:
    std::map<GridIndex,GridData*> data;
    double size;
    std::map<GridIndex,GridData*>::iterator iter;
};

#endif
