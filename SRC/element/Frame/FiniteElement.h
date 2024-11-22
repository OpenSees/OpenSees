#pragma once
#include <array>

#include <Domain.h>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <ID.h>

class Node;
class Domain;
class Response;
class Rotation;

using namespace OpenSees;

#include <State.h>


template <int nen, int ndm, int ndf>
class FiniteElement : public Element {
public:
     FiniteElement(int tag, int classtag)
       : Element(tag, classtag),
         e_state(State::None),
         connectedExternalNodes(nen) {
         for (int i=0; i<nen; i++)
           theNodes[i] = nullptr;
     }

     FiniteElement(int tag, int classtag, std::array<int, nen>& nodes)
       : Element(tag, classtag),
         e_state(State::None),
         connectedExternalNodes(nen) {
         for (int i=0; i<nen; i++) {
           connectedExternalNodes(i) = nodes[i];
           theNodes[i] = nullptr;
         }
     }

//  ~FiniteElement();

    // For Element
    virtual const ID& getExternalNodes() final {return connectedExternalNodes;};
    virtual Node **getNodePtrs() final {return theNodes.data();};
    virtual int  getNumExternalNodes() const final {return nen;};
    virtual int  getNumDOF() final {return nen*ndf;};

    void setDomain(Domain *theDomain) final
    {
        if (theDomain == nullptr) {
         for (int i=0; i<nen; i++) {
           theNodes[i] = nullptr;
         }
         return;
        }

        for (int i=0; i<nen; i++) {
          theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
          if (theNodes[i] == nullptr) {
            opserr << "FiniteElement::setDomain  tag: " << this->getTag() << " -- Node " 
                   << connectedExternalNodes(i) << " does not exist\n";
            return;
          }

          if (theNodes[i]->getNumberDOF() != ndf) {
            opserr << "FiniteElement::setDomain  tag: " << this->getTag() << " -- Node " << connectedExternalNodes(i) 
                   << " has incorrect number of DOF\n";
            return;
          }
        }

        this->DomainComponent::setDomain(theDomain);

        if (this->setState(State::Init) != 0)
          return;

//      if (this->setState(State::Pres) != 0)
//        return;
    }

//  virtual const Vector &getResistingForce();
//  virtual const Vector &getResistingForceIncInertia();
//  virtual const Matrix &getTangentStiff();
//  virtual const Matrix &getInitialStiff();
//  virtual const Matrix &getMass();    

//  //
//  // Implemented by Children   
//  //
//  int  commitState();
//  int  revertToLastCommit();
//  int  revertToStart();
//  void update();

//  void zeroLoad();
//  const char *getClassType() const;
//  int addLoad(ElementalLoad *theLoad, double loadFactor);

//  const Vector &getResistingForce();
//  const Vector &getResistingForceIncInertia();            

//  int sendSelf(int commit, Channel &theChan);
//  int recvSelf(int commit, Channel &theChan, FEM_ObjectBroker &theBroker);

//  void Print(OPS_Stream &s, int flag =0);

//  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
//  int getResponse(int responseID, Information &eleInfo);

protected:

#ifdef FEFT
    // Implemented by children
    virtual MatrixND<ndf*nen,ndf*nen> getTangent(State state, int rate) =0;
    virtual VectorND<ndf*nen>         getForce(State state, int rate) =0;

    // Supplied to children
    const VectorND<ndf>& getNodeUnknowns(int tag, int rate);
    const VectorND<ndm>& getNodePosition(int tag, State state);
    const Rotation&      getNodeRotation(int tag, State state);
    const VectorND<ndm>& getNodeVelocity(int tag);
    const VectorND<ndm>& getNodeLocation(int tag, State state);
#endif
    virtual                       int setNodes() = 0;

    // Supplied for children
    inline int setState(State state) {

      if ((e_state & state) == state)
        return 0;

      else {
        int status = -1;
        switch (state) {
          case State::Init:
            if ((status = this->setNodes()) == 0)
              e_state |= State::Init;
            break;

          case State::Pres:
            if ((status = this->update()) == 0)
              e_state |= State::Pres;
            break;

          default:
            break;
        }
        return status;
      }
    }

    //
    std::array<Node*, nen> theNodes;

    ID  connectedExternalNodes;    

private:
    State  e_state;

};

