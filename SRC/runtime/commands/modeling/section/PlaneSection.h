//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#pragma once
#include <TaggedObject.h>


template<typename MatT>
class PlaneSection: public TaggedObject {

  public:

    PlaneSection(int tag, MatT& material, double thickness) 
    : TaggedObject(tag)
    , material(&material)
    , thickness(thickness) 
    {
    }

    double getThickness() const
    {
      return thickness;
    }

    MatT* getMaterial() const
    {
      return material;
    }

  private:
    const double thickness;
    MatT* material;
};

