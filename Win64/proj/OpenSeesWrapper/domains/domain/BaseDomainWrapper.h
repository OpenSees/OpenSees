#ifndef BaseDomainWrapper_H
#define BaseDomainWrapper_H
#pragma once
#include <Domain.h>

using namespace System;

namespace OpenSees {
	namespace Components {
		
		public ref class BaseDomainWrapper 
		{
		public:
			BaseDomainWrapper() {};
			~BaseDomainWrapper() {};

		internal:
			Domain * _Domain;
		};
	}
}
#endif