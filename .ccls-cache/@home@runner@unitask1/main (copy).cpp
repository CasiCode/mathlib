#include<bits/stdc++.h>
#include <fstream>

#include "Matrix.hpp"

// Сводит матрицу к ступенчатой форме. Вощвращает число, чтобы определить,
// вырожденная матрица или нет
int forwardElimination(Matrix<double> mat);

// Обратный ход
void backwardSubstitution(Matrix<double> mat);

// Прямой ход
void gaussianElimination(Matrix<double> mat) {
    /* Сведение к ступенчатому виду */
    int singularFlag = forwardElimination(mat);

    /* Если вырожденная */
    if (singularFlag != -1) {
        printf("Singular Matrix.\n");
        /* Если справа от "равно" у нулевой строки стоит нуль,
        то система имеет бесконечно много решений, иначе решений нет*/
        if (mat.getData()[singularFlag][mat.getRows()])
            printf("Inconsistent System.");
        else
            printf("May have infinitely many "
                   "solutions.");
        return;
    }

    /* По приведении к ступенчатому виду начинаем обратный ход */
    backwardSubstitution(mat);
}

// Операция поэлементной смены двух строк местами
void swapRows(Matrix<double> mat, int i, int j) {
    double** matrix = new double*[mat.getRows()];
    for (int m = 0; m < mat.getRows(); m++) {
        matrix[m] = new double[mat.getCols()];
        for (int n = 0; n < mat.getCols(); n++) {
            matrix[m][n] = mat.getData()[m][n];
        }
    }
    for (int k = 0; k <= mat.getRows(); k++) {
        double temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
    // Free memory
    for (int m = 0; m < mat.getRows(); m++) {
        delete[] matrix[m];
    }
    delete[] matrix;
}
/*void swapRows(Matrix<double> mat, int i, int j) {
    //printf("Swapped rows %d and %d\n", i, j);
    double matrix[mat.getRows()][mat.getCols()] = mat.getData();
    for (int k = 0; k <= mat.getRows(); k++) {
        double temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
}*/

// выводит содержимое матрицы
void print(Matrix<double> mat) {
    double** matrix = new double*[mat.getRows()];
    for (int m = 0; m < mat.getRows(); m++) {
        matrix[m] = new double[mat.getCols()];
        for (int n = 0; n < mat.getCols(); n++) {
            matrix[m][n] = mat.getData()[m][n];
        }
    }
    for (int i = 0; i < mat.getRows(); i++, printf("\n"))
        for (int j = 0; j <= mat.getRows(); j++)
            printf("%lf ", matrix[i][j]);
    printf("\n");
    // Free memory
    for (int m = 0; m < mat.getRows(); m++) {
        delete[] matrix[m];
    }
    delete[] matrix;
}

// Сводит матрицу к ступенчатому виду
int forwardElimination(Matrix<double> mat) {
    for (int k = 0; k < mat.getRows(); k++) {
        // Инициализируем значение и индекс главного элемента
        int iMax = k;
        int vMax = mat.getData()[iMax][k];

        /* ищем элемент с большей амплитудой */
        for (int i = k + 1; i < mat.getRows(); i++)
            if (fabs(mat.getData()[i][k]) > vMax)
                vMax = mat.getData()[i][k], iMax = i;

        /* Если элемент главной диагонали равен нулю,
         * это значит, что матрица вырожденная, а
         * это позже приведет к делению на нуль. */
        if (!mat.getData()[k][iMax])
            return k; // Матрица вырожденная

        /* меняем местами главную строку с текущей */
        if (iMax != k)
            swapRows(mat, k, iMax);

        for (int i = k + 1; i < mat.getRows(); i++) {
            /* Коэффициент f который установит k-й элемент текущей строки равным нулю,
             * а в последствии и k-й столбец целиком */
            double f = mat.getData()[i][k]/mat.getData()[k][k];

            /* вычитаем соответствующий элемент k-ой строки, умноженный на f */
            for (int j=k+1; j<=mat.getRows(); j++)
                mat.getData()[i][j] -= mat.getData()[k][j]*f;

            /* Заполняем нижнюю треугольную матрицу нулями */
            mat.getData()[i][k] = 0;
        }
        //print(mat);
    }
    //print(mat);
    return -1;
}

// Обратный ход
void backwardSubstitution(Matrix<double> mat) {
    double x[mat.getRows()];  // массив, хранящий решение

    /* вычисляем неизвестные начиная с последней строки и заканчивая первой */
    for (int i = mat.getRows()-1; i >= 0; i--) {
        /* начинаем с правой части уравнения */
        x[i] = mat.getData()[i][mat.getRows()];

        /* Инициализируем j как i+1 так как матрица треугольная */
        for (int j=i+1; j<mat.getRows(); j++) {
            /* вычитаем все значения с левой части,
             * кроме коэффициента при неизвестной, значение
             * которой вычисляется на этой итерации */
            x[i] -= mat.getData()[i][j]*x[j];
        }

        /* делим правую часть на оставшийся коэффициент */
        x[i] = x[i]/mat.getData()[i][i];
    }

    printf("\nSolution for the system:\n");
    for (int i=0; i<mat.getRows(); i++)
        printf("%lf\n", x[i]);
}

// Ввод из файла
void inputFromFile(std::fstream& fin, const std::string& strFilename, double) {

}

// Драйвер
int main() {
    /* Матрица на ввод */
    double mat[6][7] = {
        {0.337, 0.001, 0.015, -0.092, 0.965, 0.119, -0.306},
        {0.675, 1.203, 0.692, -1.479, 0.822, 1.358, -0.841},
        {-0.739, 0.245, 0.099, 1.152, -0.235, 0.482, 1.009},
        {0.834, -0.480, 1.429, -0.519, 0.562, -0.348, 1.345},
        {-0.902, 1.051, -0.350, 0.129, -1.470, 1.380, 0.639},
        {-0.684, 1.359, -1.010, -1.167, -1.196, -0.296, -0.488}
    };
    /*double mat[N][N+1] {
        {0.0, 0.0, 0.0, 1.0, 2.0},
        {0.0, 0.0, 1.0, 5.0, 1.0},
        {0.0, 1.0, 3.0, 1.0, 1.0},
        {1.0, 12.0, 2.0, 6.0, 3.0}
    };*/

    //gaussianElimination(readDoubleMatrix(std::fstream("input.txt"), "input.txt"));
    gaussianElimination(Matrix<double>(6, 7, mat));

    return 0;
}