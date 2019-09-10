#pragma once
#include <SP_Constraint.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MP_Constraint.h>
#include <Pressure_Constraint.h>
#include <joint/MP_Joint2D.h>
#include <joint/MP_Joint3D.h>

#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../components/DomainComponentWrapper.h"

namespace OpenSees {
	namespace Components {
		namespace Constraints {
			public ref class SP_ConstraintWrapper : DomainComponentWrapper
			{
			public:
				SP_ConstraintWrapper();
				SP_ConstraintWrapper(int nodeTag, int ndof, double value, bool isConstant);
				int GetNodeTag() {
					return _SP_Constraint->getNodeTag();
				}

				int GetDOF_Number() {
					return _SP_Constraint->getDOF_Number();
				}

				double GetValue() {
					return _SP_Constraint->getValue();
				}

			

				~SP_ConstraintWrapper() {
					if (_SP_Constraint != 0)
						delete _SP_Constraint;
				};
			internal:
				SP_Constraint * _SP_Constraint;
			};

			public ref class ImposedMotionSPWrapper : SP_ConstraintWrapper
			{
			public:
				ImposedMotionSPWrapper(int nodeTag, int ndof, int patternTag, int theGroundMotionTag);
				~ImposedMotionSPWrapper() {
					if (_SP_Constraint != 0)
						delete _SP_Constraint;
				};
			internal:

			};

			public ref class ImposedMotionSP1Wrapper : SP_ConstraintWrapper
			{
			public:
				ImposedMotionSP1Wrapper(int nodeTag, int ndof, int patternTag, int groundMotionTag);
				~ImposedMotionSP1Wrapper() {
					if (_SP_Constraint != 0)
						delete _SP_Constraint;
				};
			internal:

			};

			public ref class MP_ConstraintWrapper : DomainComponentWrapper
			{
			public:
				MP_ConstraintWrapper() {};
				MP_ConstraintWrapper(int nodeRetain,
					int nodeConstr,
					IDWrapper^ constrainedDOF,
					IDWrapper^ retainedDOF,
					int classTag);
				MP_ConstraintWrapper(int nodeRetain,
					int nodeConstr,
					MatrixWrapper^ constrnt,
					IDWrapper^ constrainedDOF,
					IDWrapper^ retainedDOF);

				~MP_ConstraintWrapper() {
					if (_MP_Constraint != 0)
						delete _MP_Constraint;
				};
			internal:
				MP_Constraint * _MP_Constraint;
			};

		}
	}
}

//namespace OpenSees {
//	namespace MultiPoint {
//		public ref class MP_ConstraintWrapper : DomainComponentWrapper
//		{
//		public:
//			MP_ConstraintWrapper() {};
//			MP_ConstraintWrapper(int nodeRetain,
//				int nodeConstr,
//				IDWrapper^ constrainedDOF,
//				IDWrapper^ retainedDOF,
//				int classTag);
//			MP_ConstraintWrapper(int nodeRetain,
//				int nodeConstr,
//				MatrixWrapper^ constrnt,
//				IDWrapper^ constrainedDOF,
//				IDWrapper^ retainedDOF);
//			void static CreateRigidBeam(DomainWrapper^ domain, int nR, int nC) {
//				RigidBeam* rb = new RigidBeam(*domain->_Domain, nR, nC);
//			}
//
//			static void CreateRigidDiaphragm(DomainWrapper^ domain, int nodeR, IDWrapper^ nodeC,
//				int perpDirnToPlaneConstrained) {
//				RigidDiaphragm* rb = new RigidDiaphragm(*domain->_Domain, nodeR, *nodeC->_ID, perpDirnToPlaneConstrained);
//			}
//
//			static void CreateRigidRod(DomainWrapper^ domain, int nodeR, int nodeC) {
//				RigidRod* rb = new RigidRod(*domain->_Domain, nodeR, nodeC);
//			}
//
//			~MP_ConstraintWrapper();
//		internal:
//			MP_Constraint * _MP_Constraint;
//		};
//	}
//}
//
//namespace OpenSees {
//	namespace Pressure {
//		public ref class Pressure_ConstraintWrapper : DomainComponentWrapper
//		{
//		public:
//			Pressure_ConstraintWrapper() {};
//			Pressure_ConstraintWrapper(int classTag, int nodeId, int ptag, int ndf) {};
//			Pressure_ConstraintWrapper(int nodeId, int ptag, int ndf) {};
//			Pressure_ConstraintWrapper(int nodeId, int ndf) {};
//			~Pressure_ConstraintWrapper() {};
//		internal:
//			Pressure_Constraint * _Pressure_Constraint;
//		};
//	}
//}
