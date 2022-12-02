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

// Massimo Petracca - ASDEA Software, Italy (2022)
//
// Series3DMaterial is an aggregation
// of NDMaterials (3D only) objects all considered acting in Series.
//

#ifndef Series3DMaterial_h
#define Series3DMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>
#include <map>
#include <memory>

namespace Series3DUtils {
	class SolverWrapper;
	class ResponseWrapper;
}

class Series3DMaterial : public NDMaterial
{
public:
	enum IterativeTangentType {
		IT_Tangent,
		IT_Initial,
		IT_StabilizedTangent
	};

public:
	// life-cycle
	Series3DMaterial(
		int tag,
		const std::vector<NDMaterial*>& theMaterials,
		const std::vector<double>& theWeights,
		int theMaxNumberOfIterations = 10,
		double theRelativeTolerance = 1.0e-4,
		double theAbsoluteTolerance = 1.0e-8,
		bool verbose = false);
	Series3DMaterial();
	~Series3DMaterial();

	// info
	const char* getClassType(void) const { return "Series3DMaterial"; };

	// density
	double getRho(void);

	// set state
	int setTrialStrain(const Vector& strain);

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
	virtual int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

private:
	bool imposeIsoStressCondition(IterativeTangentType ittype);
	double computeResidualNorm();
	const Matrix& computeDenominator(IterativeTangentType ittype);
	const Vector& computeWeightedStrainResidual();
	bool solveForLagrangeMultipliers(
		const Vector& ewg,
		IterativeTangentType ittype,
		Series3DUtils::SolverWrapper& solver);
	bool solveForStrainVectors(
		const Vector& ewg,
		IterativeTangentType ittype,
		Series3DUtils::SolverWrapper& solver,
		const Matrix &D);
	const Matrix& getMaterialTangent(NDMaterial* mat, IterativeTangentType ittype) const;
	void computeHomogenizedTangent(IterativeTangentType ittype);
	void computeHomogenizedStress();

private:
	// the list of NDMaterials (size = number of materials)
	std::vector<NDMaterial*> m_materials;
	// the list of weights (size = number of materials)
	std::vector<double> m_weights;
	// the maximum number of iterations
	int m_max_iter = 10;
	// the relative tolerance (iso-stress constraint)
	double m_rel_tol = 1.0e-4;
	// the absolute tolerance (iso-stress constraint)
	double m_abs_tol = 1.0e-8;
	// the verbose flag
	bool m_verbose = false;
	// the unknown lagrange multipliers (size = strain-size, 6 for 3D problems)
	Vector m_lambda = Vector(6);
	Vector m_lambda_commit = Vector(6);
	// the homogenized strain, stress and tangent
	Vector m_strain = Vector(6);
	Vector m_strain_commit = Vector(6);
	Vector m_stress = Vector(6);
	Vector m_stress_commit = Vector(6);
	Matrix m_tangent = Matrix(6, 6);
	Matrix m_initial_tangent = Matrix(6, 6);
	// stabilization term
	double m_stab = 0.01;
	// needed to restore the last committed stage on each material
	// before solving for the iso-stress condition.
	std::vector<Vector> m_mat_strain_commit;
	// responses for homogenized outputs
	std::map<int, std::shared_ptr<Series3DUtils::ResponseWrapper>> m_response_map;
};
#endif