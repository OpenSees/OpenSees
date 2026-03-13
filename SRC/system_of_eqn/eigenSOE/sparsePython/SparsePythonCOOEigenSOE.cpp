#include "SparsePythonCOOEigenSOE.h"

#include "SparsePythonCOOEigenSolver.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <ID.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <Vertex.h>
#include <VertexIter.h>

#include <algorithm>

SparsePythonCOOEigenSOE::SparsePythonCOOEigenSOE(SparsePythonCOOEigenSolver &theSolver)
    : EigenSOE(theSolver, EigenSOE_TAGS_SparsePythonCOOEigenSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED)
{
    theSolver.setEigenSOE(*this);
}

SparsePythonCOOEigenSOE::SparsePythonCOOEigenSOE()
    : EigenSOE(EigenSOE_TAGS_SparsePythonCOOEigenSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED)
{
}

SparsePythonCOOEigenSOE::~SparsePythonCOOEigenSOE() = default;

int
SparsePythonCOOEigenSOE::getNumEqn(void) const
{
    return static_cast<int>(rowOffsets.size() > 0 ? rowOffsets.size() - 1 : 0);
}

void
SparsePythonCOOEigenSOE::buildCOOPattern(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    if (size <= 0) {
        rowIndices.clear();
        colIndices.clear();
        kValues.clear();
        mValues.clear();
        rowOffsets.assign(1u, 0);
        return;
    }

    const std::size_t nnz = static_cast<std::size_t>(size) + 2 * theGraph.getNumEdge();

    rowIndices.resize(nnz);
    colIndices.resize(nnz);
    kValues.assign(nnz, 0.0);
    mValues.assign(nnz, 0.0);
    rowOffsets.assign(static_cast<std::size_t>(size) + 1u, 0);

    int last = 0;
    rowOffsets[0] = 0;
    for (int row = 0; row < size; ++row) {
        Vertex *theVertex = theGraph.getVertexPtr(row);
        if (theVertex == nullptr) {
            opserr << "WARNING: SparsePythonCOOEigenSOE::buildCOOPattern - "
                   << "vertex " << row << " not found in graph, size set to 0" << endln;
            rowIndices.assign(0u, 0);
            colIndices.assign(0u, 0);
            kValues.clear();
            mValues.clear();
            rowOffsets.assign(1u, 0);
            return;
        }

        const ID &theAdjacency = theVertex->getAdjacency();
        ID columns(0, theAdjacency.Size() + 1);

        columns.insert(theVertex->getTag());
        for (int i = 0; i < theAdjacency.Size(); ++i) {
            columns.insert(theAdjacency(i));
        }

        for (int i = 0; i < columns.Size(); ++i) {
            rowIndices[static_cast<std::size_t>(last)] = row;
            colIndices[static_cast<std::size_t>(last)] = columns(i);
            ++last;
        }

        rowOffsets[static_cast<std::size_t>(row + 1)] = last;
    }

    if (static_cast<std::size_t>(last) < nnz) {
        rowIndices.resize(static_cast<std::size_t>(last));
        colIndices.resize(static_cast<std::size_t>(last));
        kValues.resize(static_cast<std::size_t>(last), 0.0);
        mValues.resize(static_cast<std::size_t>(last), 0.0);
    }
}

int
SparsePythonCOOEigenSOE::setSize(Graph &theGraph)
{
    buildCOOPattern(theGraph);

    EigenSolver *theSolver = this->getSolver();
    if (theSolver != nullptr) {
        const int solverOK = theSolver->setSize();
        if (solverOK < 0) {
            opserr << "WARNING: SparsePythonCOOEigenSOE::setSize - solver failed setSize()" << endln;
            return solverOK;
        }
    }

    matrixStatus = MatrixStatus::STRUCTURE_CHANGED;
    return 0;
}

void
SparsePythonCOOEigenSOE::updateMatrixStatus()
{
    if (matrixStatus != MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }
}

int
SparsePythonCOOEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "SparsePythonCOOEigenSOE::addA - Matrix and ID not of similar sizes\n";
        return -1;
    }

    const int size = getNumEqn();
    if (size == 0 || rowOffsets.size() != static_cast<std::size_t>(size) + 1u) {
        return 0;
    }

    for (int i = 0; i < idSize; ++i) {
        const int row = id(i);
        if (row < 0 || row >= size) {
            continue;
        }

        const int start = rowOffsets[static_cast<std::size_t>(row)];
        const int end = rowOffsets[static_cast<std::size_t>(row + 1)];

        for (int j = 0; j < idSize; ++j) {
            const int col = id(j);
            if (col < 0 || col >= size) {
                continue;
            }

            for (int k = start; k < end; ++k) {
                if (colIndices[static_cast<std::size_t>(k)] == col) {
                    kValues[static_cast<std::size_t>(k)] += fact * m(i, j);
                    break;
                }
            }
        }
    }

    updateMatrixStatus();
    return 0;
}

int
SparsePythonCOOEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "SparsePythonCOOEigenSOE::addM - Matrix and ID not of similar sizes\n";
        return -1;
    }

    const int size = getNumEqn();
    if (size == 0 || rowOffsets.size() != static_cast<std::size_t>(size) + 1u) {
        return 0;
    }

    for (int i = 0; i < idSize; ++i) {
        const int row = id(i);
        if (row < 0 || row >= size) {
            continue;
        }

        const int start = rowOffsets[static_cast<std::size_t>(row)];
        const int end = rowOffsets[static_cast<std::size_t>(row + 1)];

        for (int j = 0; j < idSize; ++j) {
            const int col = id(j);
            if (col < 0 || col >= size) {
                continue;
            }

            for (int k = start; k < end; ++k) {
                if (colIndices[static_cast<std::size_t>(k)] == col) {
                    mValues[static_cast<std::size_t>(k)] += fact * m(i, j);
                    break;
                }
            }
        }
    }

    updateMatrixStatus();
    return 0;
}

void
SparsePythonCOOEigenSOE::zeroA(void)
{
    std::fill(kValues.begin(), kValues.end(), 0.0);
    if (matrixStatus == MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }
}

void
SparsePythonCOOEigenSOE::zeroM(void)
{
    std::fill(mValues.begin(), mValues.end(), 0.0);
    if (matrixStatus == MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }
}

int
SparsePythonCOOEigenSOE::sendSelf(int, Channel &)
{
    return 0;
}

int
SparsePythonCOOEigenSOE::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
SparsePythonCOOEigenSOE::setPythonSolver(SparsePythonCOOEigenSolver &newSolver)
{
    this->theSolver = &newSolver;
    return 0;
}


