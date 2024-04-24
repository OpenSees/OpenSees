/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Massimo Petracca - ASDEA Software, Italy (2024)
//
// Parallel3DMaterial is an aggregation
// of NDMaterials (3D only) objects all considered acting in Parallel.
//

#ifndef Parallel3DMaterial_h
#define Parallel3DMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>
#include <map>
#include <memory>

namespace Parallel3DUtils {
	class ResponseWrapper;
}

class Parallel3DMaterial : public NDMaterial
{
public:
	// life-cycle
	Parallel3DMaterial(
		int tag,
		const std::vector<NDMaterial*>& theMaterials,
		const std::vector<double>& theWeights);
	Parallel3DMaterial();
	~Parallel3DMaterial();

	// info
	const char* getClassType(void) const { return "Parallel3DMaterial"; };

	// density
	double getRho(void);

	// set state
	int setTrialStrain(const Vector& v);
	int setTrialStrain(const Vector& v, const Vector& r);

	// get state
	const Vector& getStrain(void);
	const Vector& getStress(void);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* code);
	const char* getType(void) const;
	int getOrder(void) const;
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

private:
	void computeInitialTangent();
	void computeTangent();
	void computeStress();

private:
	// the list of NDMaterials (size = number of materials)
	std::vector<NDMaterial*> m_materials;
	// the list of weights (size = number of materials)
	std::vector<double> m_weights;
	// data
	Vector m_strain = Vector(6);
	Vector m_strain_commit = Vector(6);
	Vector m_stress = Vector(6);
	Matrix m_tangent = Matrix(6, 6);
	Matrix m_initial_tangent = Matrix(6, 6);
	// responses for homogenized outputs
	std::map<int, std::shared_ptr<Parallel3DUtils::ResponseWrapper>> m_response_map;
};
#endif