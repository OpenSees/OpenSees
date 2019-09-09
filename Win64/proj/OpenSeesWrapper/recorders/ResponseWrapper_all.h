#pragma once

#include "ResponseWrapper.h"
#include <elementresponse.h>
#include <fiberresponse.h>
#include <materialresponse.h>

#include "../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../materials/ndmaterial/NDMaterialWrapper.h"
#include "../materials/section/fiber/FiberWrapper.h"
#include "../elements/ElementWrapper.h"
#include "../matrix/IDWrapper.h"
#include "../matrix/MatrixWrapper.h"
#include "../matrix/VectorWrapper.h"


using namespace System;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;
using namespace OpenSees::Materials::Sections::Repres;

namespace OpenSees {
	namespace Recorders {
		
		public ref class ElementResponseWrapper : ResponseWrapper
		{
		public:
			ElementResponseWrapper(ElementWrapper^ ele, int id) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, int val) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, double val) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, IDWrapper^ val) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, VectorWrapper^ val) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, MatrixWrapper^ val) ;
			ElementResponseWrapper(ElementWrapper^ ele, int id, VectorWrapper^ val1, IDWrapper^ val2) ;
			~ElementResponseWrapper() {
				if (_Response != 0)
					delete _Response;
			};
		internal:
			
		};

		public ref class FiberResponseWrapper : ResponseWrapper
		{
		public:
			FiberResponseWrapper(FiberWrapper^ fib, int id) ;
			FiberResponseWrapper(FiberWrapper^ fib, int id, int val) ;
			FiberResponseWrapper(FiberWrapper^ fib, int id, double val) ;
			FiberResponseWrapper(FiberWrapper^ fib, int id, IDWrapper^ val) ;
			FiberResponseWrapper(FiberWrapper^ fib, int id, VectorWrapper^ val) ;
			FiberResponseWrapper(FiberWrapper^ fib, int id, MatrixWrapper^ val) ;
			~FiberResponseWrapper() {
				if (_Response != 0)
					delete _Response;
			};
		internal:

		};

		public ref class NDMaterialResponseWrapper : ResponseWrapper
		{
		public:
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id) ;
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id, int val) ;
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id, double val) ;
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id, IDWrapper^ val) ;
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id, VectorWrapper^ val) ;
			NDMaterialResponseWrapper(NDMaterialWrapper^ mat, int id, MatrixWrapper^ val) ;
			~NDMaterialResponseWrapper() {
				if (_Response != 0)
					delete _Response;
			};
		internal:

		};

		public ref class UniaxialMaterialResponseWrapper : ResponseWrapper
		{
		public:
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id);
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id, int val);
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id, double val);
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id, IDWrapper^ val);
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id, VectorWrapper^ val);
			UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ mat, int id, MatrixWrapper^ val);
			~UniaxialMaterialResponseWrapper() {
				if (_Response != 0)
					delete _Response;
			};
		internal:

		};
	}
}