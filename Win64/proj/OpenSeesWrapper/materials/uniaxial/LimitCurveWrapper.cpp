#include "stdafx.h"
#include "LimitCurveWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Components;

AxialCurveWrapper::AxialCurveWrapper(int tag, int eleTag, BaseDomainWrapper^ theDomain,
	double Fsw, double Kdeg, double Fres,
	int defType, int forType,
	int ndI, int ndJ, int dof, int perpDirn,
	double delta, int eleRemove) {
	_LimitCurve = new AxialCurve(0, tag, eleTag, theDomain->_Domain, Fsw, Kdeg, Fres, defType, forType, ndI, ndJ, dof, perpDirn, delta, eleRemove);
}
//
RotationShearCurveWrapper::RotationShearCurveWrapper(int crvTag, int eleTag,
	int ndI, int ndj, int rotAxis, double Vn, double Vr, double Kdeg, double rotLim, int defTyp,
	double b, double d, double h, double L, double st, double As, double Acc, double ld, double db, double rhot, double fc,
	double fy, double fyt, double delta, BaseDomainWrapper^ theDom, BaseElementWrapper^ theEle, NodeWrapper^ theNdI, NodeWrapper^ theNdJ) {
	_LimitCurve = new RotationShearCurve(crvTag, eleTag, ndI, ndj, rotAxis, Vn, Vr, Kdeg, rotLim, defTyp, b, d, h, L, st, As, Acc, ld, db, rhot, fc, fy, fyt, delta, theDom->_Domain, theEle->_Element, theNdI->_Node, theNdJ->_Node);
}

ThreePointCurveWrapper::ThreePointCurveWrapper(int tag, int eleTag, BaseDomainWrapper^ theDomain,
	double x1, double y1,
	double x2, double y2,
	double x3, double y3,
	double Kdeg, double Fres,
	int defType, int forType,
	int ndI, int ndJ, int dof, int perpDirn) {
	_LimitCurve = new ThreePointCurve(tag, eleTag, theDomain->_Domain, x1, y1, x2, y2, x3, y3, Kdeg, Fres, defType, forType, ndI, ndJ, dof, perpDirn);
}










