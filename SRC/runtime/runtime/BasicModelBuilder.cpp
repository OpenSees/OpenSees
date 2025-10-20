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
// A BasicModelBuilder stores intermediate "reference" objects like
// materials and sections that are used
// to construct other objects like elements.
//
// Written: cmp
//
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <initializer_list>
#include <string>
#include <unordered_map>

#include <modeling/commands.h>

#include <runtimeAPI.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>

#include <CrdTransf.h>
#include <BasicModelBuilder.h>

#include <tcl.h> // For TCL_OK/ERROR


BasicModelBuilder::BasicModelBuilder(Domain &domain,
                                     Tcl_Interp *interp, 
                                     int NDM, int NDF)
  : ndm(NDM), ndf(NDF), theInterp(interp),
    section_builder_is_set(false),
    theDomain(&domain),
    tclEnclosingPattern(nullptr),
    next_node_load(0)
    // , next_elem_load(0)
{
  using namespace OpenSees;

  static int ncmd = sizeof(ModelBuilderCommands)/sizeof(decltype(ModelBuilderCommands[0])); // CommandTableEntry);

  Tcl_CreateCommand(interp, "wipe", TclCommand_wipeModel, (ClientData)this, nullptr);

  for (int i = 0; i < ncmd; i++)
    Tcl_CreateCommand(interp, 
        ModelBuilderCommands[i].name, 
        ModelBuilderCommands[i].func, 
        (ClientData) this, nullptr);

  tclEnclosingPattern = nullptr;

  Tcl_SetAssocData(interp, "OPS::theTclBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theBasicModelBuilder", NULL, (ClientData)this);
  Tcl_SetAssocData(interp, "OPS::theTclDomain", NULL, (ClientData)&theDomain);

}

BasicModelBuilder::~BasicModelBuilder()
{

  for (auto& [part, val] : m_registry ) {
    for (auto& [tag, obj] : val)
      delete obj;
  }

  // set the pointers to 0
  theDomain = nullptr;
  tclEnclosingPattern = nullptr;

  using namespace OpenSees;

  static int ncmd = sizeof(ModelBuilderCommands)/sizeof(decltype(ModelBuilderCommands[0]));
  for (int i = 0; i < ncmd; i++)
    Tcl_DeleteCommand(theInterp, ModelBuilderCommands[i].name);
}


int
BasicModelBuilder::buildFE_Model()
{
  return 0;
}


int
BasicModelBuilder::getNDM() const
{
  return ndm;
}

int
BasicModelBuilder::getNDF() const
{
  return ndf;
}


void
BasicModelBuilder::letClobber(bool let_clobber)
{
  no_clobber = !let_clobber;
}

bool
BasicModelBuilder::canClobber()
{
  return !no_clobber;
}

int
BasicModelBuilder::incrNodalLoadTag()
{
  return ++next_node_load;
}

int
BasicModelBuilder::decrNodalLoadTag()
{
  return --next_node_load;
}

int
BasicModelBuilder::getNodalLoadTag() 
{
  return   next_node_load;
}


int
BasicModelBuilder::addSP_Constraint(int axisDirn, double axisValue, const ID &fixityCodes, double tol)
{
  return theDomain->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
}

LoadPattern *
BasicModelBuilder::getEnclosingPattern()
{
  return tclEnclosingPattern;
}

int
BasicModelBuilder::setEnclosingPattern(LoadPattern* pat)
{
  tclEnclosingPattern = pat;
  return 1;
}

int
BasicModelBuilder::getCurrentSectionBuilder(int& tag)
{
  if (section_builder_is_set) {
    tag = current_section_builder;
    return  0;
  } else
    return -1;
}

void 
BasicModelBuilder::setCurrentSectionBuilder(int tag)
{
  section_builder_is_set   = true;
  current_section_builder  = tag;
}

Domain *
BasicModelBuilder::getDomain() const 
{
  return theDomain;
}


int 
BasicModelBuilder::printRegistry(const char *partition, OPS_Stream& stream, int flag) const 
{
  int count = 0;
  auto iter = m_registry.find(partition);
  if (iter == m_registry.end()) {
    return count;
  }

  for (auto const& [key, val] : iter->second) {
    if (count != 0)
      stream << ",\n";

    val->Print(stream, flag);
    count++;
  }

  return count;
}

void* 
BasicModelBuilder::getRegistryObject(const char* type, const char* specialize, int tag, int flags) const
{
  std::string partition = std::string{type};
  if (specialize)
    partition += std::string{specialize};

  auto iter = m_registry.find(partition);
  if (iter == m_registry.end()) {
    if (flags == 0)
      opserr << "No objects of type \"" << partition
             << "\" have been created.\n";
    return nullptr;
  }

  auto iter_objs = iter->second.find(tag) ;
  if (iter_objs == iter->second.end()) {
    if (flags == 0)
      opserr << "No object with tag \"" << tag << "\" in partition \"" 
             << partition << "\"\n";
    return nullptr;
  }

  return (void*)iter_objs->second;
}

int
BasicModelBuilder::addRegistryObject(const char* type, const char* specialize, int tag, void *obj)
{
  std::string partition = std::string{type};
  if (specialize)
    partition += std::string{specialize};

  m_registry[partition][tag] = (TaggedObject*)obj;
  return TCL_OK;
}

int
BasicModelBuilder::findFreeTag(const char* partition, int& tag) const
{
  tag = 0;
  const auto iter = m_registry.find(std::string{partition});
  // If we dont have a table with partition name, no objects
  // have been created and tag = 0 works; return success.
  if (iter == m_registry.end())
    return 0;


  // Otherwise, find something larger than all existing tags
  const std::unordered_map<int, TaggedObject*>& table = iter->second;
  for (auto const& [key, val] : table)
    if (key > tag)
      tag = key + 1;

  return 0;
}

int
BasicModelBuilder::removeRegistryObject(const char* partition, int tag, int flags) 
{
  const auto iter = m_registry.find(std::string{partition});
  if (iter == m_registry.end()) {
    if (flags == 0)
      opserr << "No objects of type \"" << partition
             << "\" have been created.\n";
    return -1;
  }
  std::unordered_map<int, TaggedObject*>& table = iter->second;
  if (table.find(tag) != table.end()) {
      table.erase(tag);
      return 0;
  }
  return -1;
}


