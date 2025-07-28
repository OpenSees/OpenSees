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
// Written: Claudio M. Perez
//
#pragma once
#include <map>
#include <TaggedObject.h>
#include <VectorND.h>
#include <Vector3D.h>

#include <FrameTransform.h>
#include <LinearFrameTransf.h>
#include <SouzaFrameTransf.h>
#include <PDeltaFrameTransf3d.h>
#include <EuclidFrameTransf.h>
#include <Isometry/RankinIsometry.h>
#include <Isometry/CrisfieldIsometry.h>
#include <Isometry/BattiniIsometry.h>

namespace OpenSees {

class FrameTransformBuilder : public TaggedObject {
public:
    FrameTransformBuilder(int ndm, int t, const char *n) 
    : ndm(ndm), 
      TaggedObject(t), 
      vz{{0, 0, 0}}, offsets{}, offset_flags(0) {
      snprintf(name, sizeof(name), "%s", n);
    }

    virtual ~FrameTransformBuilder() {}
  
    template<int nn, int ndf>
    FrameTransform<nn, ndf> *
    create()
    {
      std::array<Vector3D, nn> *offset_array = nullptr;
      std::array<Vector3D, nn>  offset_data{};

      if (offsets.size() > 0) {
        offset_array = &offset_data;
        for (int i = 0; i < nn; ++i) {
          auto it = offsets.find(i+1);
          if (it != offsets.end())
            (*offset_array)[i] = it->second;
        }
      }

      int tag = this->getTag();
      if (strstr(name, "Linear") != nullptr)
        return new LinearFrameTransf<nn, ndf> (tag, vz, offset_array, offset_flags);

      else if (strcmp(name, "Corotational") == 0) {
        if constexpr (ndf == 6)
          return new SouzaFrameTransf<nn, ndf> (tag, vz, offset_array, offset_flags);
        else 
          return nullptr;
      }

      else if (strstr(name, "PDelta") != nullptr)
        return new PDeltaFrameTransf<nn, ndf> (tag, vz, offset_array, offset_flags);

      else if (strcmp(name, "Corotational02") == 0 || strcmp(name, "Isometric") == 0 || strstr(name, "Rigid") != nullptr)
      {
        if (getenv("Battini"))
          return new EuclidFrameTransf<nn, ndf, BattiniIsometry<nn>> (tag, vz, offset_array, offset_flags);
        else if (getenv("Crisfield"))
          return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,true>> (tag, vz, offset_array, offset_flags);
        else if (getenv("Crisfield02"))
          return new EuclidFrameTransf<nn, ndf, CrisfieldIsometry<nn,false>> (tag, vz, offset_array, offset_flags);
        else
          return new EuclidFrameTransf<nn, ndf, RankinIsometry<nn>> (tag, vz, offset_array, offset_flags);
      }

      return nullptr;
    }

    virtual void 
    Print(OPS_Stream&s, int flag) {
      if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"" << name << "\"";
        s << ", ";
        s << "\"vz\": [" << vz[0] << ", " << vz[1] << ", " << vz[2] << "]";
        if (offsets.size() > 0) {
          s << ", \"offsets\": {";
          for (auto it = offsets.begin(); it != offsets.end(); ++it) {
            s << it->first << ": [" << it->second[0] << ", " << it->second[1] << ", " << it->second[2] << "]";
            if (std::next(it) != offsets.end())
                s << ", ";
          }
        }
        s << "}";
      }
    }

    int ndm;
    int tag;
    char name[128];
    Vector3D vz;
    std::map<int, Vector3D> offsets;
    int offset_flags;
};

} // namespace OpenSees