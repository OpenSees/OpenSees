#ifndef SparsePythonCOOEigenSOE_h
#define SparsePythonCOOEigenSOE_h

#include "SparsePythonEigenCommon.h"

#include <EigenSOE.h>

#include <vector>

class Channel;
class FEM_ObjectBroker;
class Graph;
class ID;
class Matrix;
class SparsePythonCOOEigenSolver;

class SparsePythonCOOEigenSOE : public EigenSOE
{
  public:
    using MatrixStatus = SparsePythonEigenMatrixStatus;

    explicit SparsePythonCOOEigenSOE(SparsePythonCOOEigenSolver &theSolver);
    SparsePythonCOOEigenSOE();
    ~SparsePythonCOOEigenSOE() override;

    SparsePythonEigenStorageScheme getStorageScheme() const { return SparsePythonEigenStorageScheme::COO; }

    int getNumEqn(void) const;
    int setSize(Graph &theGraph) override;

    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addM(const Matrix &, const ID &, double fact = 1.0) override;

    void zeroA(void) override;
    void zeroM(void) override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

    int setPythonSolver(SparsePythonCOOEigenSolver &newSolver);

    MatrixStatus getMatrixStatus(void) const { return matrixStatus; }

    const std::vector<int> &getRowIndices(void) const { return rowIndices; }
    const std::vector<int> &getColIndices(void) const { return colIndices; }

    std::vector<double> &getKValues(void) { return kValues; }
    const std::vector<double> &getKValues(void) const { return kValues; }

    std::vector<double> &getMValues(void) { return mValues; }
    const std::vector<double> &getMValues(void) const { return mValues; }

    std::size_t getNNZ(void) const { return kValues.size(); }

  private:
    void buildCOOPattern(Graph &theGraph);
    void updateMatrixStatus();

    std::vector<int> rowIndices;
    std::vector<int> colIndices;
    std::vector<int> rowOffsets;

    std::vector<double> kValues;
    std::vector<double> mValues;

    MatrixStatus matrixStatus;

    friend class SparsePythonCOOEigenSolver;
};

#endif /* SparsePythonCOOEigenSOE_h */


