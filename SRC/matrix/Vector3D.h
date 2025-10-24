//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Implementation of a 3D Vector. The purpose of this file is to ensure that
// the VectorND class is as transparent to the compiler as possible.
// In particular, we want to ensure that the Vector3D class is 
// - trivially copyable (i.e., can be copied with memcpy)
// - trivial (i.e., has no user-defined constructors, destructors, or copy)
// - standard layout (i.e., has no virtual functions, no base classes, and no 
//
// Written: cmp
//
#ifndef Vector3D_h
#define Vector3D_h

#include <VectorND.h>
#include <type_traits>

using Vector3D = OpenSees::VectorND<3, double>;

static_assert(std::is_standard_layout<Vector3D>::value, "Vector3D is not standard layout.");
static_assert(std::is_aggregate<Vector3D>::value, "Vector3D is not an aggregate type.");

static_assert(std::is_trivially_copyable<Vector3D>::value, "Vector3D is not trivially copyable.");
static_assert(std::is_trivially_move_assignable<Vector3D>::value, "Vector3D is not trivially move assignable.");
static_assert(std::is_trivial<Vector3D>::value, "Vector3D is not trivial.");

static_assert(std::is_nothrow_constructible<Vector3D>::value, "Vector3D is not nothrow constructible.");
static_assert(std::is_nothrow_move_assignable<Vector3D>::value, "Vector3D is not nothrow move assignable.");

#endif // Vector3D_h
