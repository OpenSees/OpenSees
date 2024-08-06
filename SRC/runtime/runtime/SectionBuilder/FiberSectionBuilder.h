//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#ifndef FiberSectionBuilder_h 
#define FiberSectionBuilder_h 

#include <ReinfBar.h>
#include <Patch.h>
#include <ReinfLayer.h>
#include <cell/Cell.h>
#include <TaggedObject.h>

#include <FiberSection2dInt.h>

// Inherit tagged object so that the model builder can delete
class SectionBuilder: public TaggedObject {
public:
  SectionBuilder(int tag) : TaggedObject(tag) {}

  virtual int addFiber(int tag, int mat, double area, const Vector& cPos) =0;
  virtual int addHFiber(int tag, int mat, double area, const Vector& cPos)=0;

  int addPatch(const Patch& patch) {
    Cell**  cells  = patch.getCells();
    const int nc   = patch.getNumCells();
    const int mat  = patch.getMaterialID();
    for(int j=0; j<nc; j++) {
      // get fiber data
      double area        = cells[j]->getArea();
      const Vector& cPos = cells[j]->getCentroidPosition();

      if (this->addFiber(j, mat, area, cPos) != 0)
        return -1;
    }
    return 0;
  }

  int addLayer(const ReinfLayer& layer) {

    int numReinfBars   = layer.getNumReinfBars();
    ReinfBar* reinfBar = layer.getReinfBars();
    int mat            = layer.getMaterialID();

    for(int j=0; j<numReinfBars; j++) {
	// get fiber data
	double area        = reinfBar[j].getArea();
	const Vector& cPos = reinfBar[j].getPosition();
        if (this->addFiber(j, mat, area, cPos) != 0)
          return -1;
    }
    return 0;
  }

};

template <int ndm, class MatT, class SecT>
class FiberSectionBuilder : public SectionBuilder {
public:
  FiberSectionBuilder(BasicModelBuilder& builder_, SecT& section_)
    : SectionBuilder(section_.getTag()), builder(builder_), section(section_) {}

  virtual int addHFiber(int tag, int mat, double area, const Vector& cPos);
  int addFiber(int tag, int mat, double area, const Vector& cPos) {

      MatT * theMaterial = builder.getTypedObject<MatT>(mat);
      if (theMaterial == nullptr) {
        opserr << "no material with tag " << mat << " for fiber " << tag << "\n";
        return -1;
      }

      if constexpr (ndm==2) {
          section.addFiber(*theMaterial, area, cPos(0));

      } else {
          section.addFiber(*theMaterial, area, cPos(0), cPos(1));
      }
      return 0;
  }

private:
  BasicModelBuilder& builder;
  SecT&              section;
};

template <> int
FiberSectionBuilder<2, UniaxialMaterial, FiberSection2dInt>::addHFiber(int tag, int mat, double area, const Vector& cPos) {

  UniaxialMaterial * theMaterial = builder.getTypedObject<UniaxialMaterial>(mat);
  if (theMaterial == nullptr) {
    opserr << G3_ERROR_PROMPT << "no material with tag " << mat << " for fiber " << tag << "\n";
    return -1;
  }

  section.addHFiber(*theMaterial, area, cPos(0));
  return 0;
}

template <int ndm, class MatT, class SecT> int
FiberSectionBuilder<ndm, MatT, SecT>::addHFiber(int tag, int mat, double area, const Vector& cPos) {
  opserr << G3_ERROR_PROMPT << "section does not support H fibers\n";
  return -1;
}

#endif

