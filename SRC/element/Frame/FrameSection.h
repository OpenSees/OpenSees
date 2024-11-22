//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <State.h>
#include <Field.h>
#include <material/section/SectionForceDeformation.h>

struct FrameSectionConstants {
  // n-n
  double A;
  double Ay, 
         Az;
  // m-m
  double Iy, 
         Iz, 
         Iyz;
  // w-w
  double Cw, 
         Ca;
  // n-m
  double Qy, 
         Qz, 
         Qyx, 
         Qzx;
  // n-w
  double Rw, 
         Ry, 
         Rz;
  // m-w
  double Sa, 
         Sy, 
         Sz;
};

enum FrameStress : int {
  End         =     0,
  N           =     2, //
  Vy          =     3, // 0b00000010
  Vz          =     5, // 0b00000100
  T           =     6, // 0b00000000
  My          =     4, // 0b00000000
  Mz          =     1, // 0b00000000
  R           =     7, // (Obselete, also bishear)
  Q           =     8, // (Obselete, also bimoment)
  Bimoment    =     9, // 
  Wagner      =    10, // (Obselete, this is redundant)
  Bishear     =    11,
  Max         =    12,
};

typedef       int             FrameStressLayout[10];


class FrameSection : public SectionForceDeformation {

public:
  FrameSection(int tag, int clstag, double mass=0, bool use_mass=false)
    : SectionForceDeformation(tag, clstag),
      density(mass), has_mass(use_mass)
  {}

  virtual FrameSection* getFrameCopy() =0;
  virtual FrameSection* getFrameCopy(const FrameStressLayout& layout) {
    return getFrameCopy();
  }

  virtual SectionForceDeformation* getCopy() {
    return this->getFrameCopy();
  }

  virtual int getIntegral(Field field, State state, double& value) const {
    if ((field == Field::Density) && has_mass) {
        value = density;
        return 0;
    }
    return -1;
  }

  template <int n, const FrameStressLayout& scheme>
  int setTrialState(OpenSees::VectorND<n, double> e);

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n, double> getDeformation() {

    OpenSees::VectorND<n,double> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& es = this->getSectionDeformation();
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = es(j);
    }

    return sout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n, double> getResultant() {

    OpenSees::VectorND<n,double> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& s = this->getStressResultant();
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = s(j);
    }

    return sout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getTangent(State state) {

    OpenSees::MatrixND<n,n,double> kout;

    const ID& layout = this->getType();

    int m = this->getOrder();

//  const Matrix& ks = this->getSectionTangent();
    const Matrix& ks = (state == State::Init)
                      ? this->getInitialTangent()
                      : this->getSectionTangent();

    int elem_twist    = -1,
        elem_bishear  = -1;

    int sect_bishear  = -1;

    for (int i=0; i<n; i++) {
      // Save warp location
      switch (scheme[i]) {
        case FrameStress::T:
          elem_twist = i;
          break;
        case FrameStress::Bishear:
          elem_bishear = i;
        default:
        ;
      }
      for (int j=0; j<n; j++) {
        kout(i,j) = 0.0;

        for (int k=0; k<m; k++) {
          if (layout(k) == FrameStress::Bishear)
              sect_bishear  = k;

          if (layout(k) == scheme[i]) {
            for (int l=0; l<m; l++)
              if (layout(l) == scheme[j])
                kout(i,j) = ks(k,l);
          }
        }
      }
    }

    // If element has a twisting DOF and no Bishear
    // DOF, then twist == alpha, where alpha is the
    // bishear DOF.
    if (elem_twist != -1 && sect_bishear != -1 && elem_bishear == -1)
      for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
          if (layout(j) == scheme[i])
            kout(i,elem_twist) += ks(j,sect_bishear);

    return kout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getFlexibility(State state=State::Pres);

private:
  double density;
  bool has_mass;
};

template <int n, const FrameStressLayout& scheme>
int 
FrameSection::setTrialState(OpenSees::VectorND<n, double> e) {
  double strain_data[FrameStress::Max]{};

  const int m = this->getOrder();
  Vector trial(strain_data, m);

  const ID& layout = this->getType();

  int elem_twist    = -1,
      elem_bishear  = -1,
      elem_bimoment = -1;

  for (int i=0; i<n; i++) {
    // Save warp location
    switch (scheme[i]) {
      case FrameStress::T:
        elem_twist = i;
        break;
      case FrameStress::Bishear:
        elem_bishear  = i;
        break;
      case FrameStress::Bimoment:
        elem_bimoment  = i;
        break;
      default:
      ;
    }

    for (int j=0; j<m; j++)
      if (layout(j) == scheme[i])
        trial[j] = e[i];
  }

  // If element has a twisting DOF and no Bishear
  // DOF, then twist == alpha, where alpha is the
  // bishear DOF.
  // Note that elem_twist and elem_bishear are computable,
  // at compile time, so this branch can theoretically be 
  // optimized out by the compiler, however this might be 
  // optimistic
  if (elem_twist != -1 && elem_bishear == -1) {

    for (int j=0; j<m; j++)
      switch (layout(j)) {
        case FrameStress::Bishear:
          trial[j] = e[elem_twist];
          break;
        case FrameStress::Bimoment:
          if (elem_bimoment != -1)
            trial[j] = e[elem_bimoment];
          break;
      }
  }

  return this->setTrialSectionDeformation(trial);
}

template <int n, const FrameStressLayout& scheme>
OpenSees::MatrixND<n,n, double> 
FrameSection::getFlexibility(State state)
{
  OpenSees::MatrixND<n,n,double> fout;

  const ID& layout = this->getType();

  int m = this->getOrder();

  const Matrix& ks = (state == State::Init)
                    ? this->getInitialFlexibility()
                    : this->getSectionFlexibility();

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {  
      fout(i,j) = 0.0;
      for (int k=0; k<m; k++) {
        if (layout(k) == scheme[i]) {
          for (int l=0; l<m; l++)
            if (layout(l) == scheme[j]) {
              fout(i,j) = ks(k,l);
            }
        }
      }
    }
  }

  return fout;
}
