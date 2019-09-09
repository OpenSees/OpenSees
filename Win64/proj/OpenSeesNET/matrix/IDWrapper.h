#pragma once
#include <ID.h>

using namespace System;
namespace OpenSees {

	public ref class IDWrapper
	{
	public:
		IDWrapper(array<int>^ data);
		IDWrapper(int size);
		IDWrapper();
		array<int>^ ToArray() {
			int size = _ID->Size();
			array<int>^ data = gcnew array<int>(size);
			for (int i = 0; i < size; i++)
				data[i] = _ID->operator[](i);
			return data;
		}

		int Size()
		{
			return _ID->Size();
		}

		void Zero()
		{
			_ID->Zero();
		}

		int Get(int index) {
			return _ID->operator[](index);
		};

		void Set(int index, int value)
		{
			_ID->operator()(index) = value;
		};

		property int default[int]
		{
			int get(int index) {
				return _ID->operator[](index);
			}
			void set(int index, int value)
			{
				_ID->operator()(index) = value;
			}
		}

		int SetData(array<int>^ data, int selfFact, int otherFact) {
			if (this->Size() != data->Length) return -1;
			for (int i = 0; i < data->Length; i++)
				if (selfFact == 0)
				{
					if (otherFact == 1)
						this->_ID->operator()(i) = data[i];
					else
						this->_ID->operator()(i) = data[i]* otherFact;
				}
				else
				{
					if (otherFact == 1)
						this->_ID->operator()(i) = this->_ID->operator()(i)*selfFact + data[i];
					else
						this->_ID->operator()(i) = this->_ID->operator()(i)*selfFact + data[i] * otherFact;
				}
			return 0;
		}


		static array<int>^ GetArray(IDWrapper^ id) {
			int sz = id->Size();
			array<int>^ v = gcnew array<int>(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = id[i];
			}
			return v;
		}

		~IDWrapper();
		void Print() {
			this->_ID->Print();
		}
		
	internal:
		ID * _ID;
		IDWrapper(ID* id) {
			_ID = id;
		};
		IDWrapper(int* data, int size) {
			_ID = new ID(data, size);
		};

		//static const ID& GetID(IDWrapper^ id) {
		//	// TODO :: not tested
		//	int sz = id->Size();
		//	ID _id = ID(sz);
		//	for (int i = 0; i < sz; i++)
		//	{
		//		_id(i) = id[0];
		//	}
		//	return _id;
		//};

		static IDWrapper^ GetIDWrapper(ID* id) {
			if (id == 0) return gcnew IDWrapper();
			int sz = id->Size();
			IDWrapper^ v = gcnew IDWrapper(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = id->operator()(i);
			}
			return v;
		};

		static IDWrapper^ GetIDWrapper(const ID &id) {
			int sz = id.Size();
			IDWrapper^ v = gcnew IDWrapper(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = id(i);
			}
			return v;
		};



		static array<int>^ GetArray(const ID &id) {
			int sz = id.Size();
			array<int>^ v = gcnew array<int>(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = id(i);
			}
			return v;
		}

		static array<int>^ GetArray(ID *id) {
			int sz = id->Size();
			array<int>^ v = gcnew array<int>(sz);
			for (int i = 0; i < sz; i++)
			{
				v[i] = id->operator()(i);
			}
			return v;
		}


	private:

	};
}