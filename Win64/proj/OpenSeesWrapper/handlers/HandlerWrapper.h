#pragma once
#include <OPS_Stream.h>
#include <DataFileStream.h>
#include <BinaryFileStream.h>
#include <XmlFileStream.h>
#include <TCP_Stream.h>

#include "../actors/IMovableObjectWrapper.h"
#include "../matrix/VectorWrapper.h"
#include "../OPS.h"

using namespace System;
using namespace OpenSees;

namespace OpenSees {
	namespace Handlers {

		public enum class  OpenModeWrapper : int { OVERWRITE, APPEND };
		public enum class  FloatFieldWrapper { FIXEDD, SCIENTIFIC };

		public ref class OPS_StreamWrapper : IMovableObjectWrapper
		{
		public:
			OPS_StreamWrapper() {};

			// output format
			int SetFile(String^ filename, int mode, bool echo) {
				char* _filename = OPS::StringToChar(filename);
				return _OPS_StreamPtr->setFile(_filename, (openMode)mode, echo);
			};

			int SetPrecision(int precision) {
				return _OPS_StreamPtr->setPrecision(precision);
			};

			int SetFloatField(int _floatField) {
				return _OPS_StreamPtr->setFloatField((floatField)_floatField);
			};

			int Precision(int precision) {
				return _OPS_StreamPtr->precision(precision);
			};

			int Width(int width) {
				return _OPS_StreamPtr->width(width);
			};

			int Tag(String^ tagName) {
				return _OPS_StreamPtr->tag(OPS::StringToChar(tagName));
			};

			int Tag(String^ tagName, String^ value) {
				return _OPS_StreamPtr->tag(OPS::StringToChar(tagName), OPS::StringToChar(value));
			};

			int EndTag() {
				return _OPS_StreamPtr->endTag();
			};

			int Attr(String^ name, int n) {
				return _OPS_StreamPtr->attr(OPS::StringToChar(name), n);
			};

			int Attr(String^ name, double val) {
				return _OPS_StreamPtr->attr(OPS::StringToChar(name), val);
			};

			int Attr(String^ name, String^ value) {
				return _OPS_StreamPtr->attr(OPS::StringToChar(name), OPS::StringToChar(value));
			};

			int Write(VectorWrapper^ vec) {
				return _OPS_StreamPtr->write(*vec->_Vector);
			};

			~OPS_StreamWrapper() {
				if (_OPS_StreamPtr != 0)
					delete _OPS_StreamPtr;
			};
		internal:
			OPS_Stream * _OPS_StreamPtr;
		private:

		};

		public ref class DataFileStreamWrapper : OPS_StreamWrapper
		{
		public:
			DataFileStreamWrapper(String^ filename);
			DataFileStreamWrapper(String^ filename, int mode, int indent, int doCSV,
				bool closeOnWrite, int precision, bool doScientific);
			~DataFileStreamWrapper() {
				if (_OPS_StreamPtr != 0)
					delete _OPS_StreamPtr;
			};


		};

		public ref class BinaryFileStreamWrapper : OPS_StreamWrapper
		{
		public:
			BinaryFileStreamWrapper(String^ filename);
			BinaryFileStreamWrapper(String^ filename, int mode);
			~BinaryFileStreamWrapper() {
				if (_OPS_StreamPtr != 0)
					delete _OPS_StreamPtr;
			};


		};

		public ref class XmlFileStreamWrapper : OPS_StreamWrapper
		{
		public:
			XmlFileStreamWrapper(String^ filename, int mode);
			~XmlFileStreamWrapper() {
				if (_OPS_StreamPtr != 0)
					delete _OPS_StreamPtr;
			};
		};

		public ref class TCP_StreamWrapper : OPS_StreamWrapper
		{
		public:
			TCP_StreamWrapper(unsigned int other_Port,
				String^ other_InetAddr,
				bool checkEndianness);
			~TCP_StreamWrapper() {
				if (_OPS_StreamPtr != 0)
					delete _OPS_StreamPtr;
			};
		};
	}
}