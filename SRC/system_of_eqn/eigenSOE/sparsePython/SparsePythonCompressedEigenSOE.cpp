#include "SparsePythonCompressedEigenSOE.h"

#include "SparsePythonCompressedEigenSolver.h"

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

SparsePythonCompressedEigenSOE::SparsePythonCompressedEigenSOE(
    SparsePythonCompressedEigenSolver &theSolver,
    StorageScheme scheme)
    : EigenSOE(theSolver, EigenSOE_TAGS_SparsePythonCompressedEigenSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED),
      storageScheme(scheme)
{
    if (storageScheme == StorageScheme::COO) {
        opserr << "WARNING: SparsePythonCompressedEigenSOE does not support COO; defaulting to CSR" << endln;
        storageScheme = StorageScheme::CSR;
    }
    theSolver.setEigenSOE(*this);
}

SparsePythonCompressedEigenSOE::SparsePythonCompressedEigenSOE(StorageScheme scheme)
    : EigenSOE(EigenSOE_TAGS_SparsePythonCompressedEigenSOE),
      matrixStatus(MatrixStatus::STRUCTURE_CHANGED),
      storageScheme(scheme)
{
    if (storageScheme == StorageScheme::COO) {
        opserr << "WARNING: SparsePythonCompressedEigenSOE does not support COO; defaulting to CSR" << endln;
        storageScheme = StorageScheme::CSR;
    }
}

SparsePythonCompressedEigenSOE::~SparsePythonCompressedEigenSOE() = default;

int
SparsePythonCompressedEigenSOE::getNumEqn(void) const
{
    return static_cast<int>(indexPtr.size() > 0 ? indexPtr.size() - 1 : 0);
}

void
SparsePythonCompressedEigenSOE::buildCSRPattern(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    const std::size_t nnz = static_cast<std::size_t>(size) + 2 * theGraph.getNumEdge();

    // Build shared sparsity pattern
    indexPtr.assign(static_cast<std::size_t>(size) + 1u, 0);
    indices.resize(nnz);
    kValues.assign(nnz, 0.0);
    mValues.assign(nnz, 0.0);

    int last = 0;
    indexPtr[0] = 0;
    Vertex *theVertex = nullptr;
    for (int row = 0; row < size; ++row) {
        theVertex = theGraph.getVertexPtr(row);
        if (theVertex == nullptr) {
            opserr << "WARNING: SparsePythonCompressedEigenSOE::buildCSRPattern - "
                   << "vertex " << row << " not found in graph, "
                   << "size set to 0" << endln;
            indexPtr.assign(1u, 0);
            indices.clear();
            kValues.clear();
            mValues.clear();
            return;
        }

        const ID &theAdjacency = theVertex->getAdjacency();
        ID columns(0, theAdjacency.Size() + 1);

        columns.insert(theVertex->getTag());
        for (int i = 0; i < theAdjacency.Size(); ++i) {
            columns.insert(theAdjacency(i));
        }

        for (int i = 0; i < columns.Size(); ++i) {
            indices[static_cast<std::size_t>(last)] = columns(i);
            last++;
        }

        indexPtr[static_cast<std::size_t>(row + 1)] = last;
    }
}

void
SparsePythonCompressedEigenSOE::buildCSCPattern(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    const std::size_t nnz = static_cast<std::size_t>(size) + 2 * theGraph.getNumEdge();

    // Build shared sparsity pattern
    indexPtr.assign(static_cast<std::size_t>(size) + 1u, 0);
    indices.resize(nnz);
    kValues.assign(nnz, 0.0);
    mValues.assign(nnz, 0.0);

    int last = 0;
    indexPtr[0] = 0;
    Vertex *theVertex = nullptr;
    for (int col = 0; col < size; ++col) {
        theVertex = theGraph.getVertexPtr(col);
        if (theVertex == nullptr) {
            opserr << "WARNING: SparsePythonCompressedEigenSOE::buildCSCPattern - "
                   << "vertex " << col << " not found in graph, "
                   << "size set to 0" << endln;
            indexPtr.assign(1u, 0);
            indices.clear();
            kValues.clear();
            mValues.clear();
            return;
        }

        const ID &theAdjacency = theVertex->getAdjacency();
        ID rows(0, theAdjacency.Size() + 1);

        rows.insert(theVertex->getTag());
        for (int i = 0; i < theAdjacency.Size(); ++i) {
            rows.insert(theAdjacency(i));
        }

        for (int i = 0; i < rows.Size(); ++i) {
            indices[static_cast<std::size_t>(last)] = rows(i);
            last++;
        }

        indexPtr[static_cast<std::size_t>(col + 1)] = last;
    }
}

int
SparsePythonCompressedEigenSOE::setSize(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    if (size < 0) {
        opserr << "WARNING: SparsePythonCompressedEigenSOE::setSize - "
               << "size of SOE < 0" << endln;
        return -1;
    }

    if (storageScheme == StorageScheme::CSR) {
        buildCSRPattern(theGraph);
    } else if (storageScheme == StorageScheme::CSC) {
        buildCSCPattern(theGraph);
    } else {
        opserr << "WARNING: SparsePythonCompressedEigenSOE::setSize - unknown storage scheme\n";
        return -1;
    }

    EigenSolver *theSolver = this->getSolver();
    if (theSolver != nullptr) {
        const int solverOK = theSolver->setSize();
        if (solverOK < 0) {
            opserr << "WARNING: SparsePythonCompressedEigenSOE::setSize -"
                   << "solver failed setSize()" << endln;
            return solverOK;
        }
    }

    matrixStatus = MatrixStatus::STRUCTURE_CHANGED;

    return 0;
}

SparsePythonCompressedEigenSOE::MatrixStatus
SparsePythonCompressedEigenSOE::getMatrixStatus() const
{
    return matrixStatus;
}

void
SparsePythonCompressedEigenSOE::updateMatrixStatus()
{
    if (matrixStatus != MatrixStatus::STRUCTURE_CHANGED) {
        matrixStatus = MatrixStatus::COEFFICIENTS_CHANGED;
    }
}

int
SparsePythonCompressedEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "SparsePythonCompressedEigenSOE::addA - Matrix and ID not of similar sizes\n";
        return -1;
    }

    const int size = getNumEqn();

    if (storageScheme == StorageScheme::CSR) {
        if (fact == 1.0) {
            for (int i = 0; i < idSize; ++i) {
                const int row = id(i);
                if (row < 0 || row >= size) {
                    continue;
                }
    
                const int start = indexPtr[static_cast<std::size_t>(row)];
                const int end = indexPtr[static_cast<std::size_t>(row + 1)];
    
                for (int j = 0; j < idSize; ++j) {
                    const int col = id(j);
                    if (col < 0 || col >= size) {
                        continue;
                    }
    
                    for (int k = start; k < end; ++k) {
                        if (indices[static_cast<std::size_t>(k)] == col) {
                            kValues[static_cast<std::size_t>(k)] += m(i, j);
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

                const int start = indexPtr[static_cast<std::size_t>(row)];
                const int end = indexPtr[static_cast<std::size_t>(row + 1)];

                for (int j = 0; j < idSize; ++j) {
                    const int col = id(j);
                    if (col < 0 || col >= size) {
                        continue;
                    }

                    for (int k = start; k < end; ++k) {
                        if (indices[static_cast<std::size_t>(k)] == col) {
                            kValues[static_cast<std::size_t>(k)] += fact * m(i, j);
                            break;
                        }
                    }
                }
            }
        }
    } else { // CSC
        if (fact == 1.0) {
            for (int j = 0; j < idSize; ++j) {
                const int col = id(j);
                if (col < 0 || col >= size) {
                    continue;
                }

                const int start = indexPtr[static_cast<std::size_t>(col)];
                const int end = indexPtr[static_cast<std::size_t>(col + 1)];

                auto beginIt = indices.begin() + start;
                auto endIt = indices.begin() + end;

                for (int i = 0; i < idSize; ++i) {
                    const int row = id(i);
                    if (row < 0 || row >= size) {
                        continue;
                    }

                    auto it = std::lower_bound(beginIt, endIt, row);
                    if (it != endIt && *it == row) {
                        const std::size_t offset = static_cast<std::size_t>(it - indices.begin());
                        kValues[offset] += m(i, j);
                    }
                }
            }
        } else {
            for (int j = 0; j < idSize; ++j) {
                const int col = id(j);
                if (col < 0 || col >= size) {
                    continue;
                }

                const int start = indexPtr[static_cast<std::size_t>(col)];
                const int end = indexPtr[static_cast<std::size_t>(col + 1)];

                auto beginIt = indices.begin() + start;
                auto endIt = indices.begin() + end;

                for (int i = 0; i < idSize; ++i) {
                    const int row = id(i);
                    if (row < 0 || row >= size) {
                        continue;
                    }

                    auto it = std::lower_bound(beginIt, endIt, row);
                    if (it != endIt && *it == row) {
                        const std::size_t offset = static_cast<std::size_t>(it - indices.begin());
                        kValues[offset] += fact * m(i, j);
                    }
                }
            }
        }
    }

    updateMatrixStatus();
    return 0;
}

int
SparsePythonCompressedEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "SparsePythonCompressedEigenSOE::addM - Matrix and ID not of similar sizes\n";
        return -1;
    }

    const int size = getNumEqn();

    if (storageScheme == StorageScheme::CSR) {
        if (fact == 1.0) {
            for (int i = 0; i < idSize; ++i) {
                const int row = id(i);
                if (row < 0 || row >= size) {
                    continue;
                }
    
                const int start = indexPtr[static_cast<std::size_t>(row)];
                const int end = indexPtr[static_cast<std::size_t>(row + 1)];
    
                for (int j = 0; j < idSize; ++j) {
                    const int col = id(j);
                    if (col < 0 || col >= size) {
                        continue;
                    }
    
                    for (int k = start; k < end; ++k) {
                        if (indices[static_cast<std::size_t>(k)] == col) {
                            mValues[static_cast<std::size_t>(k)] += m(i, j);
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

                const int start = indexPtr[static_cast<std::size_t>(row)];
                const int end = indexPtr[static_cast<std::size_t>(row + 1)];

                for (int j = 0; j < idSize; ++j) {
                    const int col = id(j);
                    if (col < 0 || col >= size) {
                        continue;
                    }

                    for (int k = start; k < end; ++k) {
                        if (indices[static_cast<std::size_t>(k)] == col) {
                            mValues[static_cast<std::size_t>(k)] += fact * m(i, j);
                            break;
                        }
                    }
                }
            }
        }
    } else { // CSC
        if (fact == 1.0) {
            for (int j = 0; j < idSize; ++j) {
                const int col = id(j);
                if (col < 0 || col >= size) {
                    continue;
                }

                const int start = indexPtr[static_cast<std::size_t>(col)];
                const int end = indexPtr[static_cast<std::size_t>(col + 1)];

                auto beginIt = indices.begin() + start;
                auto endIt = indices.begin() + end;

                for (int i = 0; i < idSize; ++i) {
                    const int row = id(i);
                    if (row < 0 || row >= size) {
                        continue;
                    }

                    auto it = std::lower_bound(beginIt, endIt, row);
                    if (it != endIt && *it == row) {
                        const std::size_t offset = static_cast<std::size_t>(it - indices.begin());
                        mValues[offset] += m(i, j);
                    }
                }
            }
        } else {
            for (int j = 0; j < idSize; ++j) {
                const int col = id(j);
                if (col < 0 || col >= size) {
                    continue;
                }

                const int start = indexPtr[static_cast<std::size_t>(col)];
                const int end = indexPtr[static_cast<std::size_t>(col + 1)];

                auto beginIt = indices.begin() + start;
                auto endIt = indices.begin() + end;

                for (int i = 0; i < idSize; ++i) {
                    const int row = id(i);
                    if (row < 0 || row >= size) {
                        continue;
                    }

                    auto it = std::lower_bound(beginIt, endIt, row);
                    if (it != endIt && *it == row) {
                        const std::size_t offset = static_cast<std::size_t>(it - indices.begin());
                        mValues[offset] += fact * m(i, j);
                    }
                }
            }
        }
    }

    updateMatrixStatus();
    return 0;
}

void
SparsePythonCompressedEigenSOE::zeroA(void)
{
    std::fill(kValues.begin(), kValues.end(), 0.0);
    updateMatrixStatus();
}

void
SparsePythonCompressedEigenSOE::zeroM(void)
{
    std::fill(mValues.begin(), mValues.end(), 0.0);
    updateMatrixStatus();
}

int
SparsePythonCompressedEigenSOE::setPythonSolver(SparsePythonCompressedEigenSolver &newSolver)
{
    return setSolver(newSolver);
}

int
SparsePythonCompressedEigenSOE::sendSelf(int commitTag, Channel &theChannel)
{
    // Not implemented yet
    return 0;
}

int
SparsePythonCompressedEigenSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // Not implemented yet
    return 0;
}

