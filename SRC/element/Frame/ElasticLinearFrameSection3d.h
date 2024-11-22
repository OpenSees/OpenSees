//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef ElasticLinearFrameSection3d_h
#define ElasticLinearFrameSection3d_h

#include <array>
#include <memory>
#include <MatrixND.h>
#include <FrameSection.h>

class Matrix;
class Vector;
class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class ElasticLinearFrameSection3d : public FrameSection
{

 public:
  // Constructor like ElasticSection3d
  ElasticLinearFrameSection3d(int tag, double E, double A, 
                              double Iz, double Iy, double G, double J,
                              double mass, bool use_mass)
    : ElasticLinearFrameSection3d(tag, 
                                    E,   G, 
                                    A, 0.0, 0.0,
                                   Iy,  Iz, 0.0,
                                  0.0, Iy+Iz-J,       // Cw, Ca
                                  0.0, 0.0, 0.0, 0.0, // Q
                                  0.0, 0.0, 0.0,      // R
                           -(Iy+Iz-J), 0.0, 0.0,      // S
                                  mass, use_mass)
  {}

  // Constructor like ElasticSection2d
  ElasticLinearFrameSection3d(int tag, double E, double A, double Iz, double mass, bool use_mass)
    : ElasticLinearFrameSection3d(tag, 
                                    E, 0.0,           // Material
                                    A, 0.0, 0.0,      // Axial and Shear
                                  0.0,  Iz, 0.0,      // Flexure
                                  0.0, 0.0,           // Warping (Cw, Ca)
                                  0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0,
                                  mass, use_mass)
  {}

  // Full constructor
  ElasticLinearFrameSection3d(int tag, 
      double E, 
      double G, 
      // n-n
      double A,
      double Ay,
      double Az,
      // m-m
      double Iy,
      double Iz,
      double Iyz,
      // w-w
      double Cw,
      double Ca,
      // n-m
      double Qy,
      double Qz,
      double Qyx,
      double Qzx,
      // n-w
      double Rw,
      double Ry,
      double Rz,
      // m-w
      double Sa,
      double Sy,
      double Sz,
      //
      double mass,
      bool   use_mass
  );

  ElasticLinearFrameSection3d();
  ~ElasticLinearFrameSection3d();

  const char *getClassType() const {return "ElasticLinearFrameSection3d";};
  
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation();
  
  virtual int getIntegral(Field field, State state, double&) const override final;

  const Vector &getStressResultant();
  const Matrix &getSectionTangent();
  const Matrix &getInitialTangent();
  const Matrix &getSectionFlexibility();
  const Matrix &getInitialFlexibility();
  
  FrameSection *getFrameCopy();
  virtual FrameSection* getFrameCopy(const FrameStressLayout& layout);
  const ID &getType();
  int getOrder() const;
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
               FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
                                              bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);

 protected:
  
 private:
  
  void getConstants(FrameSectionConstants& consts) const; 

  double E, G;

  constexpr static int nr = 8;

  Vector  v;
  Matrix  M,         // Generic matrix for returning Ks or Fs (nr x nr)
         *Ksen;      // Tangent sensitivity (nr x nr)

  OpenSees::VectorND<nr> s;
  OpenSees::VectorND<nr> e;                      // section trial deformations
  std::shared_ptr<OpenSees::MatrixND<nr,nr>> Ks;
  Matrix* Fs = nullptr;

  int parameterID;

  static ID layout;
  std::array<double, 2> centroid;
};

#endif
