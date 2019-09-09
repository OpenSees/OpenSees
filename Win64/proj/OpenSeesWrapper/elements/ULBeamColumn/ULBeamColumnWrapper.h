#pragma once

#include <CyclicModel.h>
#include <BilinearCyclic.h>
#include <LinearCyclic.h>
#include <QuadraticCyclic.h>
#include <Elastic2dGNL.h>
#include <UpdatedLagrangianBeam2D.h>
#include <Elastic2dGNL.h>
#include <Inelastic2DYS01.h>
#include <Inelastic2DYS02.h>
#include <Inelastic2DYS03.h>
#include <InelasticYS2DGNL.h>

#include <FourNodeQuadWithSensitivity.h>
#include <NineNodeMixedQuad.h>

#include "../ElementWrapper.h"
#include "../../materials/yieldSurface/YieldSurfaceWrapper.h";
#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::YieldSurfaces;
namespace OpenSees {
	namespace Elements {
#pragma region CyclicModel
		public ref class CyclicModelWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			CyclicModelWrapper() {};
			~CyclicModelWrapper() {
				if (_CyclicModel != 0)
					delete _CyclicModel;
			};
		internal:
			CyclicModel * _CyclicModel;
		};

		public ref class BilinearCyclicWrapper abstract : CyclicModelWrapper
		{
		public:
			BilinearCyclicWrapper(int tag, double weight) {};
			~BilinearCyclicWrapper() {
				if (_CyclicModel != 0)
					delete _CyclicModel;
			};
		internal:
			
		};

		public ref class LinearCyclicWrapper abstract : CyclicModelWrapper
		{
		public:
			LinearCyclicWrapper(int tag) {};
			~LinearCyclicWrapper() {
				if (_CyclicModel != 0)
					delete _CyclicModel;
			};
		internal:
			
		};

		public ref class QuadraticCyclicWrapper abstract : CyclicModelWrapper
		{
		public:
			QuadraticCyclicWrapper(int tag, double wt, double qy) {};
			~QuadraticCyclicWrapper() {
				if (_CyclicModel != 0)
					delete _CyclicModel;
			};
		internal:
			
		};

#pragma endregion
#pragma region beamcolumns


		public ref class UpdatedLagrangianBeam2DWrapper abstract : ElementWrapper
		{
		public:
			UpdatedLagrangianBeam2DWrapper() {};
			~UpdatedLagrangianBeam2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		public ref class InelasticYS2DGNLWrapper abstract : UpdatedLagrangianBeam2DWrapper
		{
		public:
			InelasticYS2DGNLWrapper() {};
			InelasticYS2DGNLWrapper(int tag,
				int Nd1, int Nd2,
				YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
				int rf_algo, bool islinear, double rho) {};
			~InelasticYS2DGNLWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		public ref class Elastic2dGNLWrapper : UpdatedLagrangianBeam2DWrapper
		{
		public:
			Elastic2dGNLWrapper(int tag, double A, double E, double I, int Nd1, int Nd2,
				bool islinear, double rho);
			~Elastic2dGNLWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		public ref class Inelastic2DYS01Wrapper : InelasticYS2DGNLWrapper
		{
		public:
			Inelastic2DYS01Wrapper(int tag, double A, double E, double Iz,
				int Nd1, int Nd2,
				YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
				int rf_algo , bool islinear , double rho );
			~Inelastic2DYS01Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		public ref class Inelastic2DYS02Wrapper : InelasticYS2DGNLWrapper
		{
		public:
			Inelastic2DYS02Wrapper(int tag, double A, double E, double Iz,
				int Nd1, int Nd2,
				YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
				CyclicModelWrapper^ cycModel,
				double del_p_max,
				double Alpha, double Beta, int rf_algo,
				bool islinear, double rho);
			~Inelastic2DYS02Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		public ref class Inelastic2DYS03Wrapper : InelasticYS2DGNLWrapper
		{
		public:
			Inelastic2DYS03Wrapper(int tag, double a_ten, double a_com, double e,
				double iz_pos, double iz_neg, int Nd1, int Nd2,
				YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
				int rf_algo, bool islinear, double rho);
			~Inelastic2DYS03Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		internal:
		};

		

#pragma endregion

	}
}
