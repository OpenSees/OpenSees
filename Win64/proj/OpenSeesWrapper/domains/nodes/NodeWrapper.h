#pragma once
#include <Node.h>
#include <Vector.h>
#include <DOF_Group.h>
#include "../components/DomainComponentWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../matrix/IDWrapper.h"

using namespace OpenSees;
using namespace System;
namespace OpenSees {
	namespace Components {
		public ref class NodeWrapper : DomainComponentWrapper
		{
		internal:
			Node * _Node;
			NodeWrapper(Node * _Node) {
				this->_Node = _Node;
				this->_TaggedObject = _Node;
			};
		public:
			NodeWrapper();
			NodeWrapper(int tag, int ndf, double crd1);
			NodeWrapper(int tag, int ndf, double crd1, double crd2);
			NodeWrapper(int tag, int ndf, double crd1, double crd2, double crd3);
			~NodeWrapper();
			
			int GetNumberDOF(void) {
				return this->_Node->getNumberDOF();
			};

			VectorWrapper^ GetCrdsVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getCrds());
			}

			array<double>^ GetCrds()
			{
				return VectorWrapper::GetArray(_Node->getCrds());
			}

			// commit data
			VectorWrapper^ GetCommitDispVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getDisp());
			}

			array<double>^ GetCommitDisp()
			{
				return VectorWrapper::GetArray(_Node->getDisp());
			}

			VectorWrapper^ GetCommitVelVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getVel());
			}

			array<double>^ GetCommitVel()
			{
				return VectorWrapper::GetArray(_Node->getVel());
			}

			VectorWrapper^ GetCommitAccelVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getAccel());
			}

			array<double>^ GetCommitAccel()
			{
				return VectorWrapper::GetArray(_Node->getAccel());
			}

			// trial data
			VectorWrapper^ GetTrialDispVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getTrialDisp());
			}

			array<double>^ GetTrialDisp()
			{
				return VectorWrapper::GetArray(_Node->getTrialDisp());
			}

			VectorWrapper^ GetTrialVelVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getTrialVel());
			}

			array<double>^ GetTrialVel()
			{
				return VectorWrapper::GetArray(_Node->getTrialVel());
			}

			VectorWrapper^ GetTrialAccelVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getTrialAccel());
			}

			array<double>^ GetTrialAccel()
			{
				return VectorWrapper::GetArray(_Node->getTrialAccel());
			}

			// incr data
			VectorWrapper^ GetIncrDispVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getIncrDisp());
			}

			array<double>^ GetIncrDisp()
			{
				return VectorWrapper::GetArray(_Node->getIncrDisp());
			}

			VectorWrapper^ GetIncrDeltaDispVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getIncrDeltaDisp());
			}

			array<double>^ GetIncrDeltaDisp()
			{
				return VectorWrapper::GetArray(_Node->getIncrDeltaDisp());
			}

			// RV
			VectorWrapper^ GetRVVector(VectorWrapper^ accel)
			{
				return VectorWrapper::GetVectorWrapper(_Node->getRV(*accel->_Vector));
			}

			array<double>^ GetRV(VectorWrapper^ accel)
			{
				return VectorWrapper::GetArray(_Node->getRV(*accel->_Vector));
			}

			// ID
			IDWrapper^ GetIDWrapper() {			
				DOF_Group* dofgrp = _Node->getDOF_GroupPtr();
				return IDWrapper::GetIDWrapper(dofgrp->getID());
			};

			array<int>^ GetIds()
			{
				DOF_Group* dofgrp = _Node->getDOF_GroupPtr();
				return IDWrapper::GetArray(dofgrp->getID());
			}

			// reactions
			VectorWrapper^ GetReactionsVector()
			{
				return VectorWrapper::GetVectorWrapper(_Node->getReaction());
			}

			array<double>^ GetReactions()
			{
				return VectorWrapper::GetArray(_Node->getReaction());
			}

			// eigen vectors
			MatrixWrapper^ GetEigenvectorsMatrix()
			{
				return MatrixWrapper::GetMatrixWrapper(_Node->getEigenvectors());
			}

			array<double,2>^ GetEigenvectors()
			{
				return MatrixWrapper::GetArray(_Node->getEigenvectors());
			}

			// commit revert 
			int CommitState() {
				return _Node->commitState();
			}

			int RevertToLastCommit() {
				return _Node->revertToLastCommit();
			}

			int RevertToStart() {
				return _Node->revertToStart();
			}

			/*VectorWrapper^  GetDisp() {
				Vector mat = _Node->getDisp();
				VectorWrapper^ wmat = gcnew VectorWrapper();
				wmat->_Vector = &mat;
				return wmat;
			};

			VectorWrapper^  GetVel() {
				Vector mat = _Node->getVel();
				VectorWrapper^ wmat = gcnew VectorWrapper();
				wmat->_Vector = &mat;
				return wmat;
			};

			VectorWrapper^  GetAccel() {
				Vector mat = _Node->getAccel();
				VectorWrapper^ wmat = gcnew VectorWrapper();
				wmat->_Vector = &mat;
				return wmat;
			};

			VectorWrapper^  GetReactions() {
				Vector mat = _Node->getReaction();
				VectorWrapper^ wmat = gcnew VectorWrapper();
				wmat->_Vector = &mat;
				return wmat;
			};*/

			



			[ObsoleteAttribute("Use GetCommitDisp")]
			array<double>^ GetDisp()
			{
				Vector vec = _Node->getDisp();
				int size = vec.Size();
				array<double>^ _vec = gcnew array<double>(size);
				for (int i = 0; i < size; i++)
					_vec[i] = vec[i];
				return _vec;
			}

			[ObsoleteAttribute("Use GetCommitVelc")]
			array<double>^ GetVel()
			{
				Vector vec = _Node->getVel();
				int size = vec.Size();
				array<double>^ _vec = gcnew array<double>(size);
				for (int i = 0; i < size; i++)
					_vec[i] = vec[i];

				return _vec;
			}

			[ObsoleteAttribute("Use GetCommitAcc")]
			array<double>^ GetAcc()
			{
				Vector vec = _Node->getAccel();
				int size = vec.Size();
				array<double>^ _vec = gcnew array<double>(size);
				for (int i = 0; i < size; i++)
					_vec[i] = vec[i];
				return _vec;
			}
		};
	}
}



