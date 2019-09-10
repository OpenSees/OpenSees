#pragma once
#include <information.h>

#include "../OPS.h"
#include "../matrix/IDWrapper.h"
#include "../matrix/MatrixWrapper.h"
#include "../matrix/VectorWrapper.h"


using namespace System;
using namespace OpenSees;


namespace OpenSees {
	namespace Recorders {
		public enum class InfoTypeWrapper {
			UnknownType = 0, IntType = 1, DoubleType = 2,
			IdType = 3, VectorType = 4, MatrixType = 5
		};
		public ref class InformationWrapper 
		{
		public:
			InformationWrapper();
			InformationWrapper(int val);
			InformationWrapper(double val);
			InformationWrapper(IDWrapper^ val);
			InformationWrapper(VectorWrapper^ val);
			InformationWrapper(MatrixWrapper^ val);
			InformationWrapper(IDWrapper^ val1, VectorWrapper^ val2);
			int SetInt(int val) {
				return _Information->setInt(val);
			}

			int SetDouble(double val) {
				return _Information->setDouble(val);
			}

			int SetID(IDWrapper^ id) {
				return _Information->setID(*id->_ID);
			}

			int SetVector(VectorWrapper^ vec) {
				return _Information->setVector(*vec->_Vector);
			}

			int SetMatrix(MatrixWrapper^ mat) {
				return _Information->setMatrix(*mat->_Matrix);
			}

			int SetString(String^ str) {
				return _Information->setString(OPS::StringToChar(str));
			}

			void Print() {
				_Information->Print(opserr);
			}

			~InformationWrapper() {
				if (_Information != 0)
					delete _Information;
			};

			property InfoTypeWrapper^ Type {
				InfoTypeWrapper^ get() { 
					int type = _Information->theType;
					return static_cast<InfoTypeWrapper>(type);
				}
			}

			property int Int {
				int get() { return _Information->theInt; }
			}

			property double Double {
				double get() { return _Information->theDouble; }
			}

			property IDWrapper^ Id {
				IDWrapper^ get() { 
					IDWrapper^ _id = gcnew IDWrapper(_Information->theID);
					return _id;
				}
			}

			property VectorWrapper^ Vector {
				VectorWrapper^ get() {
					VectorWrapper^ _id = gcnew VectorWrapper(_Information->theVector);
					return _id;
				}
			}

			property MatrixWrapper^ Matrix {
				MatrixWrapper^ get() {
					MatrixWrapper^ _id = gcnew MatrixWrapper(_Information->theMatrix);
					return _id;
				}
			}

			property String^ String {
				System::String^ get() {
					System::String^ _id = gcnew System::String(_Information->theString);
					return _id;
				}
			}

		internal:
			Information *_Information;
		};
	}
}