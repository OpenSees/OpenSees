#include "SparsePythonCOOLinSOE.h"

#include "SparsePythonCOOLinSolver.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <ID.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <Vertex.h>
#include <VertexIter.h>

#include <algorithm>
#include <cmath>

SparsePythonCOOLinSOE::SparsePythonCOOLinSOE(SparsePythonCOOLinSolver &theSolver)
    : LinearSOE(theSolver, LinSOE_TAGS_SparsePythonCOOLinSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED)
{
    theSolver.setLinearSOE(*this);
    updateVectorViews();
}

SparsePythonCOOLinSOE::SparsePythonCOOLinSOE()
    : LinearSOE(LinSOE_TAGS_SparsePythonCOOLinSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED)
{
    updateVectorViews();
}

SparsePythonCOOLinSOE::~SparsePythonCOOLinSOE() = default;

int
SparsePythonCOOLinSOE::getNumEqn(void) const
{
    return static_cast<int>(sol.size());
}

void
SparsePythonCOOLinSOE::buildCOOPattern(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    if (size < 0) {
        opserr << "WARNING: SparsePythonCOOLinSOE::buildCOOPattern - invalid graph size" << endln;
        rowIndices.clear();
        colIndices.clear();
        values.clear();
        rhs.clear();
        sol.clear();
        rowOffsets.clear();
        updateVectorViews();
        return;
    }

    rhs.assign(static_cast<std::size_t>(size), 0.0);
    sol.assign(static_cast<std::size_t>(size), 0.0);

    const std::size_t nnz =
        static_cast<std::size_t>(size) + 2u * static_cast<std::size_t>(theGraph.getNumEdge());

    rowOffsets.assign(static_cast<std::size_t>(size) + 1u, 0);
    rowIndices.resize(nnz);
    colIndices.resize(nnz);
    values.assign(nnz, 0.0);

    std::size_t last = 0;
    rowOffsets[0] = 0;

    Vertex *theVertex = nullptr;
    for (int row = 0; row < size; ++row) {
        theVertex = theGraph.getVertexPtr(row);
        if (theVertex == nullptr) {
            opserr << "WARNING: SparsePythonCOOLinSOE::buildCOOPattern - vertex "
                   << row << " not found; resetting system" << endln;
            rowOffsets.assign(1u, 0);
            rowIndices.clear();
            colIndices.clear();
            values.clear();
            rhs.clear();
            sol.clear();
            updateVectorViews();
            return;
        }

        const ID &adjacency = theVertex->getAdjacency();
        ID columns(0, adjacency.Size() + 1);

        columns.insert(theVertex->getTag());
        for (int i = 0; i < adjacency.Size(); ++i) {
            columns.insert(adjacency(i));
        }

        for (int i = 0; i < columns.Size(); ++i) {
            rowIndices[last] = row;
            colIndices[last] = columns(i);
            ++last;
        }

        rowOffsets[static_cast<std::size_t>(row + 1)] = static_cast<int>(last);
    }

    if (last < nnz) {
        rowIndices.resize(last);
        colIndices.resize(last);
        values.resize(last);
    }
}

int
SparsePythonCOOLinSOE::setSize(Graph &theGraph)
{
    buildCOOPattern(theGraph);
    updateVectorViews();

    LinearSOESolver *theSolver = this->getSolver();
    if (theSolver != nullptr) {
        const int solverOK = theSolver->setSize();
        if (solverOK < 0) {
            opserr << "WARNING: SparsePythonCOOLinSOE::setSize - solver failed setSize()" << endln;
            return solverOK;
        }
    }

    matrixStatus = MatrixStatus::STRUCTURE_CHANGED;
    return 0;
}

int
SparsePythonCOOLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "SparsePythonCOOLinSOE::addA - Matrix and ID not of similar sizes" << endln;
        return -1;
    }

    const int size = getNumEqn();
    if (size == 0 || rowOffsets.size() != static_cast<std::size_t>(size) + 1u) {
        return 0;
    }

    if (fact == 1.0) {
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
                        values[static_cast<std::size_t>(k)] += m(i, j);
                        break;
                    }
                }
            }
        }
    } else {
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
                        values[static_cast<std::size_t>(k)] += fact * m(i, j);
                        break;
                    }
                }
            }
        }
    }

    if (matrixStatus != MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }

    return 0;
}

int
SparsePythonCOOLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != v.Size()) {
        opserr << "SparsePythonCOOLinSOE::addB - Vector and ID not of similar sizes" << endln;
        return -1;
    }

    for (int i = 0; i < idSize; ++i) {
        const int pos = id(i);
        if (pos < 0 || pos >= getNumEqn()) {
            continue;
        }
        rhs[static_cast<std::size_t>(pos)] += fact * v(i);
    }
    return 0;
}

int
SparsePythonCOOLinSOE::setB(const Vector &v, double fact)
{
    if (fact == 0.0) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        bVector.Zero();
        return 0;
    }

    if (v.Size() != getNumEqn()) {
        opserr << "SparsePythonCOOLinSOE::setB - incompatible sizes "
               << getNumEqn() << " and " << v.Size() << endln;
        return -1;
    }

    for (int i = 0; i < getNumEqn(); ++i) {
        rhs[static_cast<std::size_t>(i)] = fact * v(i);
    }
    return 0;
}

void
SparsePythonCOOLinSOE::zeroA(void)
{
    std::fill(values.begin(), values.end(), 0.0);
    if (matrixStatus != MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }
}

void
SparsePythonCOOLinSOE::zeroB(void)
{
    std::fill(rhs.begin(), rhs.end(), 0.0);
}

const Vector &
SparsePythonCOOLinSOE::getX(void)
{
    return xVector;
}

const Vector &
SparsePythonCOOLinSOE::getB(void)
{
    return bVector;
}

double
SparsePythonCOOLinSOE::normRHS(void)
{
    double norm = 0.0;
    for (double value : rhs) {
        norm += value * value;
    }
    return std::sqrt(norm);
}

void
SparsePythonCOOLinSOE::setX(int loc, double value)
{
    if (loc >= 0 && loc < getNumEqn()) {
        sol[static_cast<std::size_t>(loc)] = value;
    }
}

void
SparsePythonCOOLinSOE::setX(const Vector &X)
{
    if (X.Size() != getNumEqn()) {
        return;
    }
    for (int i = 0; i < getNumEqn(); ++i) {
        sol[static_cast<std::size_t>(i)] = X(i);
    }
}

int
SparsePythonCOOLinSOE::solve(void)
{
    int status = this->LinearSOE::solve();
    if (status < 0) {
        return status;
    }
    matrixStatus = MatrixStatus::UNCHANGED;
    return status;
}

int
SparsePythonCOOLinSOE::sendSelf(int, Channel &)
{
    return 0;
}

int
SparsePythonCOOLinSOE::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

int
SparsePythonCOOLinSOE::setPythonSolver(SparsePythonCOOLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    if (getNumEqn() != 0) {
        const int solverOK = newSolver.setSize();
        if (solverOK < 0) {
            opserr << "WARNING: SparsePythonCOOLinSOE::setPythonSolver - solver failed setSize()" << endln;
            return -1;
        }
    }
    return this->LinearSOE::setSolver(newSolver);
}

void
SparsePythonCOOLinSOE::updateVectorViews()
{
    if (sol.empty()) {
        xVector.resize(0);
    } else {
        xVector.setData(sol.data(), static_cast<int>(sol.size()));
    }

    if (rhs.empty()) {
        bVector.resize(0);
    } else {
        bVector.setData(rhs.data(), static_cast<int>(rhs.size()));
    }
}
