#ifndef DomainWrapper_H
#define DomainWrapper_H
#pragma once
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <LoadPatternIter.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>
#include <RigidRod.h>

#include "BaseDomainWrapper.h"
#include "../../elements/ElementWrapper.h"
#include "../../recorders/RecorderWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../constraints/ConstraintWrapper.h"
#include "../patterns/LoadPatternWrapper.h"
#include "../nodes/NodeWrapper.h"
#include "../../taggeds/TaggedObjectWrapper.h"

using namespace System;
using namespace OpenSees::Elements;
using namespace OpenSees::Components;
using namespace OpenSees::Components::LoadPatterns;
using namespace OpenSees::Components::Loads;
using namespace OpenSees::Components::Constraints;
using namespace OpenSees::Recorders;
namespace OpenSees {
	namespace Components {

		public ref class DomainWrapper : BaseDomainWrapper
		{
		public:
			DomainWrapper();
			~DomainWrapper();
			bool^ AddNode(array<NodeWrapper^> ^nodes);
			bool^ AddNode(NodeWrapper^ node);
			bool^ AddElement(array<ElementWrapper^> ^elements);
			bool^ AddElement(ElementWrapper^ element);
			bool^ AddSP_Constraint(array<SP_ConstraintWrapper^> ^sps);
			bool^ AddSP_Constraint(SP_ConstraintWrapper^ sP_Constraint);
			bool^ AddMP_Constraint(array<MP_ConstraintWrapper^> ^mps);
			bool^ AddMP_Constraint(MP_ConstraintWrapper^ sP_Constraint);
			bool^ AddLoadPattern(LoadPatternWrapper^ loadPattern);
			bool^ AddNodalLoad(NodalLoadWrapper^ nodalLoad, int loadPatternTag);
			bool^ AddNodalLoad(array<NodalLoadWrapper^>^ nodalLoads, int loadPatternTag);
			bool^ AddElementLoad(ElementalLoadWrapper^ eleLoad, int loadPatternTag) {
				return _Domain->addElementalLoad(eleLoad->_ElementalLoad, loadPatternTag);
			}
			bool^ AddRecorder(RecorderWrapper^ recorder) {
				if (_Domain->addRecorder(*recorder->_Recorder) == 0)
					return true;
				else
					return false;
			}

			bool^ RemoveRecorders() {
				if (_Domain->removeRecorders() == 0)
					return true;
				else
					return false;
			}
			void ClearAll() {
				_Domain->clearAll();
			}
			ElementWrapper^ RemoveElement(int tag) {
				Element* ret = _Domain->removeElement(tag);
				if (ret == 0) return nullptr;
				ElementWrapper^ wret = gcnew ElementWrapper();
				wret->_Element = ret;
				wret->_TaggedObject = ret;
				return wret;
			}
			NodeWrapper^ RemoveNode(int tag) {
				Node* ret = _Domain->removeNode(tag);
				if (ret == 0) return nullptr;
				NodeWrapper^ wret = gcnew NodeWrapper();
				wret->_Node = ret;
				wret->_TaggedObject = ret;
				return wret;
			}

			LoadPatternWrapper^ RemoveLoadPattern(int tag) {
				LoadPattern* ret = _Domain->removeLoadPattern(tag);
				if (ret == 0) return nullptr;
				LoadPatternWrapper^ wret = gcnew LoadPatternWrapper();
				wret->_LoadPattern = ret;
				wret->_TaggedObject = ret;
				return wret;
			}

			NodeWrapper^ GetNode(int tag) {
				Node* node = _Domain->getNode(tag);
				if (node != 0)
				{
					NodeWrapper^ _node = gcnew NodeWrapper();
					_node->_Node = node;
					_node->_TaggedObject = node;
					return _node;
				}
				else
					return nullptr;
			}

			ElementWrapper^ GetElement(int tag) {
				Element* ele = _Domain->getElement(tag);
				if (ele != 0)
				{
					ElementWrapper^ _ele = gcnew ElementWrapper();
					_ele->_Element = ele;
					_ele->_TaggedObject= ele;
					return _ele;
				}
				else
					return nullptr;
			}

			LoadPatternWrapper^ GetLoadPattern(int tag) {
				LoadPattern* pattern = _Domain->getLoadPattern(tag);
				if (pattern != 0)
				{
					LoadPatternWrapper^ _pattern = gcnew LoadPatternWrapper();
					_pattern->_LoadPattern = pattern;
					_pattern->_TaggedObject = pattern;
					return _pattern;
				}
				else
					return nullptr;
			}

			int GetNumNodes() { return _Domain->getNumNodes(); }
			
			array<NodeWrapper^>^ GetNodes() {
				int num = _Domain->getNumNodes();
				array<NodeWrapper^>^ _nodes = gcnew array<NodeWrapper^>(num);
				NodeIter& iterator = _Domain->getNodes();
				Node* obj = 0;
				int i = 0;
				while ((obj = iterator()) != 0)
				{
					_nodes[i] = gcnew NodeWrapper();
					_nodes[i]->_Node = obj;
					_nodes[i]->_TaggedObject = obj;
					i++;
				}
				return _nodes;
			}

			int GetNumElements() { return _Domain->getNumElements(); }
			
			array<ElementWrapper^>^ GetElements() {
				int num = _Domain->getNumElements();
				array<ElementWrapper^>^ _elements = gcnew array<ElementWrapper^>(num);
				ElementIter& iterator = _Domain->getElements();
				Element* obj = 0;
				int i = 0;
				while ((obj = iterator()) != 0)
				{
					_elements[i] = gcnew ElementWrapper();
					_elements[i]->_Element = obj;
					_elements[i]->_TaggedObject = obj;
					i++;
				}
				return _elements;
			}

			int GetNumLoadPatterns() { return _Domain->getNumLoadPatterns(); }
			
			array<LoadPatternWrapper^>^ GetLoadPatterns() {
				int num = _Domain->getNumLoadPatterns();
				array<LoadPatternWrapper^>^ _patterns = gcnew array<LoadPatternWrapper^>(num);
				LoadPatternIter& iterator = _Domain->getLoadPatterns();
				LoadPattern* obj = 0;
				int i = 0;
				while ((obj = iterator()) != 0)
				{
					_patterns[i] = gcnew LoadPatternWrapper();
					_patterns[i]->_LoadPattern = obj;
					_patterns[i]->_TaggedObject = obj;
					i++;
				}
				return _patterns;
			}


			void RevertToStart() {
				_Domain->revertToStart();
			}

			void RevertToLastCommit() {
				_Domain->revertToLastCommit();
			}

			array<SP_ConstraintWrapper^>^ GetSPs() {
				int numNodes = _Domain->getNumSPs();
				array<SP_ConstraintWrapper^>^ _sps = gcnew array<SP_ConstraintWrapper^>(numNodes);
				SP_ConstraintIter& spIter = _Domain->getSPs();

				SP_Constraint* obj = 0;
				int i = 0;
				while ((obj = spIter()) != 0)
				{
					_sps[i] = gcnew SP_ConstraintWrapper();
					_sps[i]->_SP_Constraint = obj;
					_sps[i]->_TaggedObject = obj;
					i++;
				}
				return _sps;
			}

			array<MP_ConstraintWrapper^>^ GetMPs() {
				int numNodes = _Domain->getNumSPs();
				array<MP_ConstraintWrapper^>^ _mps = gcnew array<MP_ConstraintWrapper^>(numNodes);
				MP_ConstraintIter& mpIter = _Domain->getMPs();

				MP_Constraint* obj = 0;
				int i = 0;
				while ((obj = mpIter()) != 0)
				{
					_mps[i] = gcnew MP_ConstraintWrapper();
					_mps[i]->_MP_Constraint = obj;
					_mps[i]->_TaggedObject = obj;
				}
				return _mps;
			}

			void CreateRigidDiaphragm(int nodeR, IDWrapper^ nodeC, int perpDirnToPlaneConstrained) {
				RigidDiaphragm* rb = new RigidDiaphragm(*this->_Domain, nodeR, *nodeC->_ID, perpDirnToPlaneConstrained);
			};
			void CreateRigidBeam(int nR, int nC) {
				RigidBeam* rb = new RigidBeam(*this->_Domain, nR, nC);
			};
			void CreateRigidRod(int nodeR, int nodeC) {
				RigidRod* rb = new RigidRod(*this->_Domain, nodeR, nodeC);
			};

			void Print(int flag)
			{
				_Domain->Print(opserr, flag);
			};
			void CalculateNodalReactions(int flag)
			{
				_Domain->calculateNodalReactions(flag);
			};
			void SetLoadConst() {
				_Domain->setLoadConstant();
			}
			void SetCurrentTime(double t)
			{
				_Domain->setCurrentTime(t);
			}
			void SetCommittedTime(double t)
			{
				_Domain->setCommittedTime(t);
			}
			void SetMass(int nodeId, MatrixWrapper^ matrix) {
				_Domain->setMass(*matrix->_Matrix, nodeId);
			}
			void SetMass(array<int>^ nodeIds, MatrixWrapper^ matrix) {
				for (int i = 0; i < nodeIds->Length; i++)
					_Domain->setMass(*matrix->_Matrix, nodeIds[i]);
			}
			void SetRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc) {
				_Domain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
			}
			double GetTime()
			{
				return _Domain->getCurrentTime();
			}
			array<double>^ GetEigenValues()
			{
				Vector vec = _Domain->getEigenvalues();
				int size = vec.Size();
				array<double>^ _vec = gcnew array<double>(size);
				for (int i = 0; i < size; i++)
					_vec[i] = vec[i];
				return _vec;
			}
			
			delegate void DomainChanged(int^);
			event EventHandler^ AddNodeEventHandler;


		internal:
			void SetDomain(Domain* theDomain) {
				_Domain = theDomain;
			}
		private:
			
		};
	}
}
#endif