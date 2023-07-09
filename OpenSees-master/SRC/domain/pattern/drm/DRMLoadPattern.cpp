/* 
* @author: gnp <petropoulos@gmail.com>
 *
 * @Description: DRM Load Pattern Implementation
 *               Main functionality performed by DRMInputHandler
 *
 * @Date: 2/10/06
 */


#include "DRMLoadPattern.h"
#include <PlaneDRMInputHandler.h>
#include <LoadPattern.h>
#include <Domain.h>
#include "DRMBoundaryLayerDecorator.h"
#include "DRMInputHandler.h"
#include "Mesh3DSubdomain.h"
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <math.h>
#include <iostream>

//using namespace::std;

DRMLoadPattern::DRMLoadPattern(int tag, double cfact, DRMInputHandler* my_handler, Domain* domain)
:LoadPattern(tag, PATTERN_TAG_DRMLoadPattern)
{
  this->factor = cfact;

  // node map
  std::map<int,int> mapn;
  // ele map
  std::map<int,Element* > mape;
  // double* map for eles
  std::map<int,Vector* > maps;
  // int map for eles
  std::map<int,int > maps2;

  this->eNodes = mapn;
  this->elem = mape;
  this->storage = maps;
  this->storage2 = maps2;

  this->myDomain = domain;
  this->myHandler = my_handler;

  this->setMaps();
}


DRMLoadPattern::~DRMLoadPattern()
{
  // clean up maps
  // still need to do!
}

void
DRMLoadPattern::setMaps()
{
  this->myHandler->seteNodeMap(this->eNodes);
  this->myHandler->seteleMap(this->elem,this->storage,this->storage2);
}

void
DRMLoadPattern::applyLoad(double time)
{
  DRMBoundaryLayerDecorator *myDecorator = new DRMBoundaryLayerDecorator();
  myDecorator->setDomain(this->getDomain());
  Vector U(24);
  Vector Ud(24);
  Vector Udd(24);
  Vector load(24);
  U.Zero();
  Ud.Zero();
  Udd.Zero();
  load.Zero();
  myDecorator->setMap(this->eNodes);
  for(std::map<int,Element*>::iterator pos=this->elem.begin(); pos!=this->elem.end(); pos++) {
//    int eleTag = pos->first;
    Element* ele = (Element*) pos->second;
    if (ele != 0) {
      U.Zero();
      Ud.Zero();
      Udd.Zero();
      load.Zero();
  
      myDecorator->setBrick(ele);
      this->myHandler->getMotions(ele, time, U, Ud, Udd);
      myDecorator->applyDRMLoad(this->factor,load, U, Ud, Udd);
    }
  }
  delete myDecorator;
}
