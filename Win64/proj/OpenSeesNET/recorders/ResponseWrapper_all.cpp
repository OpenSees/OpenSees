#include "stdafx.h"

#include "ResponseWrapper_all.h"

using namespace OpenSees::Recorders;

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id) {
	_Response = new ElementResponse(ele->_Element, id);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, int val) {
	_Response = new ElementResponse(ele->_Element, id, val);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, double val) {
	_Response = new ElementResponse(ele->_Element, id, val);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, IDWrapper^ val) {
	_Response = new ElementResponse(ele->_Element, id, *val->_ID);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, VectorWrapper^ val) {
	_Response = new ElementResponse(ele->_Element, id, *val->_Vector);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, MatrixWrapper^ val) {
	_Response = new ElementResponse(ele->_Element, id, *val->_Matrix);
}

ElementResponseWrapper::ElementResponseWrapper(ElementWrapper^ ele, int id, VectorWrapper^ val1, IDWrapper^ val2) {
	_Response = new ElementResponse(ele->_Element, id, *val1->_Vector, *val2->_ID);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id) {
	_Response = new FiberResponse(ele->_Fiber, id);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id, int val) {
	_Response = new FiberResponse(ele->_Fiber, id, val);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id, double val) {
	_Response = new FiberResponse(ele->_Fiber, id, val);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id, IDWrapper^ val) {
	_Response = new FiberResponse(ele->_Fiber, id, *val->_ID);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id, VectorWrapper^ val) {
	_Response = new FiberResponse(ele->_Fiber, id, *val->_Vector);
}

FiberResponseWrapper::FiberResponseWrapper(FiberWrapper^ ele, int id, MatrixWrapper^ val) {
	_Response = new FiberResponse(ele->_Fiber, id, *val->_Matrix);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id, int val) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id, val);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id, double val) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id, val);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id, IDWrapper^ val) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id, *val->_ID);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id, VectorWrapper^ val) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id, *val->_Vector);
}

UniaxialMaterialResponseWrapper::UniaxialMaterialResponseWrapper(UniaxialMaterialWrapper^ ele, int id, MatrixWrapper^ val) {
	_Response = new MaterialResponse(ele->_UniaxialMaterial, id, *val->_Matrix);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id) {
	_Response = new MaterialResponse(ele->_NDMaterial, id);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id, int val) {
	_Response = new MaterialResponse(ele->_NDMaterial, id, val);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id, double val) {
	_Response = new MaterialResponse(ele->_NDMaterial, id, val);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id, IDWrapper^ val) {
	_Response = new MaterialResponse(ele->_NDMaterial, id, *val->_ID);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id, VectorWrapper^ val) {
	_Response = new MaterialResponse(ele->_NDMaterial, id, *val->_Vector);
}

NDMaterialResponseWrapper::NDMaterialResponseWrapper(NDMaterialWrapper^ ele, int id, MatrixWrapper^ val) {
	_Response = new MaterialResponse(ele->_NDMaterial, id, *val->_Matrix);
}



