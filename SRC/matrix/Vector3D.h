//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Implementation of a 3D Vector
//
// Written: cmp
//
#ifndef Vector3D_h
#define Vector3D_h

#include <VectorND.h>
#include <type_traits>

using Vector3D = OpenSees::VectorND<3, double>;

static_assert(std::is_trivially_copyable<Vector3D>::value, "Vector3D is not trivially copyable.");
static_assert(std::is_trivial<Vector3D>::value, "Vector3D is not trivial.");
static_assert(std::is_standard_layout<Vector3D>::value, "Vector3D is not standard layout.");
static_assert(std::is_aggregate<Vector3D>::value, "Vector3D is not an aggregate type.");

#endif // Vector3D_h
