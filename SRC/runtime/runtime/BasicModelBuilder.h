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
// Description: This file contains the class definition for
// BasicModelBuilder. A BasicModelBuilder aims to be a threadsafe
// alternative to the TclBasicBuilder class. This class adds the commands to
// create the model for the standard models that can be generated using the
// elements released with the g3 framework.
//
// Written: cmp
// Created: 10/21
//
#ifndef TCLSAFEBUILDER_H
#define TCLSAFEBUILDER_H

#include <typeinfo>
#include <string>
#include <unordered_map>

#include <TaggedObject.h>

class LoadPattern;
class MultiSupportPattern;
class OPS_Stream;
class ID;
class Domain;
struct Tcl_Interp;

class BasicModelBuilder {
public:

// Constructors / Destructors
  BasicModelBuilder(Domain &domain, Tcl_Interp *interp, int ndm, int ndf);
  ~BasicModelBuilder();



// Options
  void letClobber(bool option);
  bool canClobber();

// Methods

  int getNDM() const;
  int getNDF() const;
  Domain *getDomain() const;

  int  getCurrentSectionBuilder(int&);
  void setCurrentSectionBuilder(int);

  LoadPattern* getEnclosingPattern();
  int setEnclosingPattern(LoadPattern*);

  int incrNodalLoadTag();
  int decrNodalLoadTag();
  int getNodalLoadTag();

  //
  // Managing tagged objects
  //
  template<class T> int addTypedObject(int tag, T* obj) {
    return addRegistryObject(typeid(T).name(), nullptr, tag, obj);
  }

  template<class T, const char* specialize=nullptr> int addTaggedObject(T& obj) {
    int tag = obj.getTag();
    return addRegistryObject(typeid(T).name(), specialize, tag, &obj);
  }

  constexpr static int SilentLookup = 1;
  template <class T>
  int printRegistry(OPS_Stream& stream, int flag) const {
    return printRegistry(typeid(T).name(), stream, flag);
  }

  template<class T, const char* specialize=nullptr> T* 
  getTypedObject(int tag, int flags=0) const {
    return (T*)getRegistryObject(typeid(T).name(), specialize, tag, flags);
  }

  template<class T> int 
  removeObject(int tag, int flags=0) {
    return removeRegistryObject(typeid(T).name(), tag, flags);
  }

  template <class T> int findFreeTag(int &tag) const {
    return findFreeTag(typeid(T).name(), tag);
  }

  int addSP_Constraint(int axisDirn, 
         double axisValue, 
         const ID &fixityCodes, 
         double tol=1e-10);

  int buildFE_Model();

//
private:
  int   addRegistryObject(const char*, const char*, int tag, void* obj); 
  void* getRegistryObject(const char*, const char*, int tag, int flags) const;
  int   removeRegistryObject(const char*, int tag, int flags);
  int   findFreeTag(const char*, int& tag) const;
  int   printRegistry(const char *, OPS_Stream& stream, int flag) const ;


  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  Tcl_Interp *theInterp;
  Domain *theDomain     = nullptr;

  int next_node_load          = 0;

  // Options
  bool no_clobber = true;

  // previously extern variables
  LoadPattern *tclEnclosingPattern = nullptr;
//MultiSupportPattern *theTclMultiSupportPattern = nullptr;

  bool  section_builder_is_set   = false;
  int   current_section_builder  = 0;

// OBJECT CONTAINERS
  std::unordered_map<std::string, std::unordered_map<int, TaggedObject*>> m_registry;

};

#endif

