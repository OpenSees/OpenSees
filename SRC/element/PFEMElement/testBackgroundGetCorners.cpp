#include <iostream>
#include <vector>

using VInt = std::vector<int>;
using VVInt = std::vector<VInt>;

int ndm = 2;
int OPS_GetNDM() { return ndm; }

void getCorners(const VInt& index, int num, int dim, VVInt& indices) {
    int ndm = OPS_GetNDM();
    int counter = 0;

    if (ndm == 2) {
        indices.resize(2 * num + 1);
        int dim2 = dim + 1;
        if (dim2 >= ndm) {
            dim2 -= ndm;
        }
        for (int j = -num; j <= num; ++j) {
            indices[counter] = index;
            indices[counter][dim2] += j;
            ++counter;
        }
    } else if (ndm == 3) {
        indices.resize((2 * num + 1) * (2 * num + 1));
        int dim2 = dim + 1;
        if (dim2 >= ndm) {
            dim2 -= ndm;
        }
        int dim3 = dim + 2;
        if (dim3 >= ndm) {
            dim3 -= ndm;
        }
        for (int j = -num; j <= num; ++j) {
            for (int k = -num; k <= num; ++k) {
                indices[counter] = index;
                indices[counter][dim2] += j;
                indices[counter][dim3] += k;
                ++counter;
            }
        }
    }
}

std::ostream& operator<<(std::ostream& os, const VInt& v) {
    for (VInt::size_type i = 0; i < v.size(); i++) {
        os << v[i] << " ";
    }
    os << "\n";

    return os;
}

int main() {
    ndm = 2;
    VInt index(ndm);
    VVInt indices;

    std::cout << "2D:\n";
    std::cout << "index = " << index;
    VInt nums = {1, 2};
    VInt dims = {0, 1};
    for (int num : nums) {
        for (int dim : dims) {
            std::cout << "dim = " << dim << ", num = " << num << "\n";
            getCorners(index, num, dim, indices);
            std::cout << "num of indices = " << indices.size() << ": \n";
            for (const auto& ind: indices) {
                std::cout << "ind: " << ind;
            }
        }
    }

    ndm = 3;
    index.assign(ndm, 0);
    std::cout << "3D:\n";
    std::cout << "index = " << index;
    nums = {1, 2};
    dims = {0, 1, 2};
    for (int num : nums) {
        for (int dim : dims) {
            std::cout << "dim = " << dim << ", num = " << num << "\n";
            getCorners(index, num, dim, indices);
            std::cout << "num of indices = " << indices.size() << ": \n";
            for (const auto& ind: indices) {
                std::cout << "ind: " << ind;
            }
        }
    }

    return 0;
}