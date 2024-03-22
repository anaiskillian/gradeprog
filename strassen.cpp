#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>

using namespace std;

// addition function
vector<vector<int>> addMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return C;
}

// subtraction function
vector<vector<int>> subtractMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }

    return C;
}

// normal multiplication
vector<vector<int>> normalMult(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}
// function to pad matrix
vector<vector<int>> padMatrix(const vector<vector<int>>& A) {
    if (A.size()% 2 == 0){
        return A;
    }
    else {
        int n = A.size();
        int newSize = n + 1;
        vector<vector<int>> paddedMatrix(newSize, vector<int>(newSize, 0));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                paddedMatrix[i][j] = A[i][j];

        return paddedMatrix;
    }
}

// function to remove padding
vector<vector<int>> removePadding(const vector<vector<int>>& C, int originalSize) {
    //int n = C.size();
    vector<vector<int>> result(originalSize, vector<int>(originalSize));

    for (int i = 0; i < originalSize; i++)
        for (int j = 0; j < originalSize; j++)
            result[i][j] = C[i][j];

    return result;
}

// Strassen's algorithm
vector<vector<int>> strassen(const vector<vector<int>>& X, const vector<vector<int>>& Y) {
    // save original size to remove padding later
    int origSize = X.size();

    // Cross over point
    if (origSize <= 64)
        return normalMult(X, Y);

    // Pad matrices
    vector<vector<int>> paddedA = padMatrix(X);
    vector<vector<int>> paddedB = padMatrix(Y);

    // Continue with the padded matrices
    int newSize = paddedA.size();

    // Splitting matrices into submatrices
    int halfSize = newSize / 2;
    vector<vector<int>> A(halfSize, vector<int>(halfSize));
    vector<vector<int>> B(halfSize, vector<int>(halfSize));
    vector<vector<int>> C(halfSize, vector<int>(halfSize));
    vector<vector<int>> D(halfSize, vector<int>(halfSize));
    vector<vector<int>> E(halfSize, vector<int>(halfSize));
    vector<vector<int>> F(halfSize, vector<int>(halfSize));
    vector<vector<int>> G(halfSize, vector<int>(halfSize));
    vector<vector<int>> H(halfSize, vector<int>(halfSize));

    for (int i = 0; i < halfSize; i++) {
        for (int j = 0; j < halfSize; j++) {
          A[i][j] = paddedA[i][j];
          B[i][j] = paddedA[i][j + halfSize];
          C[i][j] = paddedA[i + halfSize][j];
          D[i][j] = paddedA[i + halfSize][j + halfSize];
          E[i][j] = paddedB[i][j];
          F[i][j] = paddedB[i][j + halfSize];
          G[i][j] = paddedB[i + halfSize][j];
          H[i][j] = paddedB[i + halfSize][j + halfSize];
        }
    }

    // intermediate matrices
    vector<vector<int>> P1 = strassen(A, subtractMatrices(F, H));
    vector<vector<int>> P2 = strassen(addMatrices(A, B), H);
    vector<vector<int>> P3 = strassen(addMatrices(C, D), E);
    vector<vector<int>> P4 = strassen(D, subtractMatrices(G, E));
    vector<vector<int>> P5 = strassen(addMatrices(A, D), addMatrices(E, H));
    vector<vector<int>> P6 = strassen(subtractMatrices(B, D), addMatrices(G, H));
    vector<vector<int>> P7 = strassen(subtractMatrices(C, A), addMatrices(E, F));

    // result submatrices
    vector<vector<int>> C11 = addMatrices(subtractMatrices(P4, P2), addMatrices(P5, P6));
    vector<vector<int>> C12 = addMatrices(P1, P2);
    vector<vector<int>> C21 = addMatrices(P3, P4);
    vector<vector<int>> C22 = addMatrices(subtractMatrices(P1, P3), addMatrices(P5, P7));

    // Combine into result matrix
    vector<vector<int>> C(newSize, vector<int>(newSize));
    for (int i = 0; i < halfSize; i++) {
        for (int j = 0; j < halfSize; j++) {
            C[i][j] = C11[i][j];
            C[i][j + halfSize] = C12[i][j];
            C[i + halfSize][j] = C21[i][j];
            C[i + halfSize][j + halfSize] = C22[i][j];
        }
    }
    // Remove padding from output
    return removePadding(C, origSize);
}


// Function to read matrices from an input file
pair<vector<vector<int>>, vector<vector<int>>> readFile(const string& filename, int dim) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(1);
    }

    vector<vector<int>> A(dim, vector<int>(dim));
    vector<vector<int>> B(dim, vector<int>(dim));

    // read matrix A 
   for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            file >> A[i][j];
        }
    }

    // Read matrix B
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            file >> B[i][j];
        }
    }

    file.close();
    return make_pair(A, B);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " flag dimension inputfile" << endl;
        return 1;
    }

    
    int flag = atoi(argv[1]);
    int dim = atoi(argv[2]);
    string file = argv[3];

    pair<vector<vector<int>>, vector<vector<int>>> matrices = readFile(file, dim);
    vector<vector<int>> A = matrices.first;
    vector<vector<int>> B = matrices.second;

    vector<vector<int>> C = strassen(A, B);
    for (int i = 0; i < dim; i++){
        cout << C[i][i] << endl;
    }

    return 0;
}