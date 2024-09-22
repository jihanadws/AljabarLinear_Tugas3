#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <cmath>

using namespace std;

const double EPSILON = 1e-10;

void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << setw(10) << fixed << setprecision(4) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

vector<vector<double>> gaussJordanElimination(vector<vector<double>> matrix) {
    int n = matrix.size();
    int m = matrix[0].size();

    for (int i = 0; i < n && i < m; i++) {
        // Mencari pivot
        int pivotRow = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(matrix[j][i]) > abs(matrix[pivotRow][i])) {
                pivotRow = j;
            }
        }

        // Menukar baris
        if (pivotRow != i) {
            swap(matrix[i], matrix[pivotRow]);
        }

        // Membuat pivot menjadi 1
        double pivot = matrix[i][i];
        if (abs(pivot) < EPSILON) {
            continue; // Skip jika pivot mendekati nol
        }
        for (int j = i; j < m; j++) {
            matrix[i][j] /= pivot;
        }

        // Eliminasi
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int j = i; j < m; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }
    }

    return matrix;
}

vector<vector<double>> findInverse(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    if (n != matrix[0].size()) {
        throw runtime_error("Matriks harus persegi untuk menghitung invers");
    }

    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n));

    // Membuat matriks augmented [A|I]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];
            augmentedMatrix[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Melakukan eliminasi Gauss-Jordan
    augmentedMatrix = gaussJordanElimination(augmentedMatrix);

    // Mengecek apakah matriks dapat diinverskan
    for (int i = 0; i < n; i++) {
        if (abs(augmentedMatrix[i][i]) < EPSILON) {
            throw runtime_error("Matriks singular, tidak dapat diinverskan");
        }
    }

    // Mengekstrak matriks invers
    vector<vector<double>> inverse(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    return inverse;
}

bool isSingular(const vector<vector<double>>& matrix) {
    try {
        findInverse(matrix);
        return false;
    } catch (const runtime_error&) {
        return true;
    }
}

int main() {
    int rows, cols;

    cout << "Masukkan jumlah baris: ";
    cin >> rows;
    cout << "Masukkan jumlah kolom: ";
    cin >> cols;

    vector<vector<double>> matrix(rows, vector<double>(cols));

    cout << "Masukkan elemen-elemen matriks:" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << "Matriks[" << i + 1 << "][" << j + 1 << "]: ";
            cin >> matrix[i][j];
        }
    }

    cout << "\nMatriks awal:" << endl;
    printMatrix(matrix);

    // Hasil eliminasi Gauss-Jordan
    cout << "Hasil eliminasi Gauss-Jordan (SPL):" << endl;
    vector<vector<double>> result = gaussJordanElimination(matrix);
    printMatrix(result);

    // Jika matriks persegi, hitung invers
    if (rows == cols) {
        if (!isSingular(matrix)) {
            vector<vector<double>> inverse = findInverse(matrix);
            cout << "Matriks invers:" << endl;
            printMatrix(inverse);
            
            cout << "Verifikasi: A * A^(-1) = I" << endl;
            vector<vector<double>> product(rows, vector<double>(cols, 0.0));
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    for (int k = 0; k < cols; k++) {
                        product[i][j] += matrix[i][k] * inverse[k][j];
                    }
                }
            }
            printMatrix(product);

            cout << "Untuk SPL Ax = b, solusi unik dapat ditemukan dengan x = A^(-1)b" << endl;
        } else {
            cout << "Matriks singular, tidak memiliki invers." << endl;
        }
    } else {
        cout << "Matriks tidak persegi, tidak dapat dihitung inversnya." << endl;
    }

    return 0;
}
