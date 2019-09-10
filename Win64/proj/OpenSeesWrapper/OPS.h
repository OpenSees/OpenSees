#pragma once
// OPS 1
using namespace System;
using namespace System::Runtime::InteropServices;
namespace OpenSees {
	public ref class OPS
	{
	public:
		static property System::Object^ Null {
			System::Object^ get() { return nullptr; }
		}
	internal:
		static char* StringToChar(String^ str) {
			char* _str = (char*)(void*)Marshal::StringToHGlobalAnsi(str);
			return _str;
		}
	};
}

