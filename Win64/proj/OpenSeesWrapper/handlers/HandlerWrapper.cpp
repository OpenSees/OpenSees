#include "stdafx.h"
#include "HandlerWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees::Handlers;

DataFileStreamWrapper::DataFileStreamWrapper(String^ filename) {
	_OPS_StreamPtr = new DataFileStream((char*)(void*)Marshal::StringToHGlobalAnsi(filename));
}

DataFileStreamWrapper::DataFileStreamWrapper(String^ filename, int mode, int indent, int doCSV,
	bool closeOnWrite, int precision, bool doScientific) {
	_OPS_StreamPtr = new DataFileStream((char*)(void*)Marshal::StringToHGlobalAnsi(filename), (openMode)mode, indent, doCSV, closeOnWrite, precision, doScientific);
}

BinaryFileStreamWrapper::BinaryFileStreamWrapper(String^ filename) {
	_OPS_StreamPtr = new BinaryFileStream((char*)(void*)Marshal::StringToHGlobalAnsi(filename));
}

BinaryFileStreamWrapper::BinaryFileStreamWrapper(String^ filename, int mode) {
	_OPS_StreamPtr = new BinaryFileStream((char*)(void*)Marshal::StringToHGlobalAnsi(filename), (openMode)mode);
}

XmlFileStreamWrapper::XmlFileStreamWrapper(String^ filename, int mode) {
	_OPS_StreamPtr = new XmlFileStream((char*)(void*)Marshal::StringToHGlobalAnsi(filename), (openMode)mode);
}

TCP_StreamWrapper::TCP_StreamWrapper(unsigned int other_Port,
	String^ other_InetAddr,
	bool checkEndianness) {
	_OPS_StreamPtr = new TCP_Stream(other_Port, (char*)(void*)Marshal::StringToHGlobalAnsi(other_InetAddr), checkEndianness);
}



