#include <iostream>
#include <fstream>
#include <cstdlib>   // For calloc and free
#include <iomanip>   // For setprecision
#include <clocale>   // For setlocale
double M = 1000000;
using namespace std;
double** pmat, ** pmat2;
int rivn, zm = 2, kilk = 1, imin2, flag = 0, type, maxSize = 3;
void Print();
void FillZ();
int MinZ();
int MaxZ();
int FindRiad(int st);
void zmina(int xmin, int ymin);
void Znaku();
bool CheckSt(int i);
int main()
{
    setlocale(LC_ALL, ""); // Enable UTF-8 for Cyrillic output
    int k, i, j, imin, jmin;
    double nmin, xmax = 0, ymax = 0, xmax2 = 0, ymax2 = 0, yg1;
    cout << fixed << setprecision(2); // Set precision for output
    cout << "Лабораторна робота №6\nЗ предмету математичні методи дослідження операцій " << endl;
    cout << "Нелінійне програмування : задачі дробово - лінійного програмування" << endl;
    cout << "Кількість рівнянь: ";
    cin >> rivn;
    rivn++;
    cout << "Кількість змінних: ";
    cin >> zm;
    maxSize++;
    pmat = (double**)calloc(rivn + 2, sizeof(double*)); // Allocate memory for array of pointers
    if (!pmat) {
        cout << "Memory allocation failed!" << endl;
        return 1;
    }
    for (k = 0; k < rivn + 2; k++) { // Allocate memory for each row
        pmat[k] = (double*)calloc(zm + rivn + maxSize, sizeof(double));
        if (!pmat[k]) {
            cout << "Memory allocation failed!" << endl;
            return 1;
        }
    }
    pmat2 = (double**)calloc(rivn + 2, sizeof(double*)); // Allocate memory for second matrix
    if (!pmat2) {
        cout << "Memory allocation failed!" << endl;
        return 1;
    }
    for (k = 0; k < rivn + 2; k++) {
        pmat2[k] = (double*)calloc(zm + rivn + maxSize, sizeof(double));
        if (!pmat2[k]) {
            cout << "Memory allocation failed!" << endl;
            return 1;
        }
    }
    for (i = 0; i < rivn + 2; i++)
        for (j = 0; j < zm + rivn + 2; j++)
            pmat[i][j] = 0;
    cout << "Введіть замість а,b,c (ax1+bx2 ? c) коефіцієнти рівняння " << endl;
    for (i = 1; i < rivn; i++) {
        cout << "Нерівність №" << i << endl;
        cin >> pmat[i][2] >> pmat[i][3] >> pmat[i][1];
        pmat[i][1] *= -1;
        cout << "Введіть знак нерівності (1- \"<=\", 2- \">=\", 3- \"=\")" << endl;
        cin >> pmat[i][rivn + zm + maxSize - 1];
        if (pmat[i][rivn + zm + maxSize - 1] == 2) {
            for (j = 1; j < rivn + zm + maxSize - 2; j++)
                if (pmat[i][j] != 0)
                    pmat[i][j] *= -1;
        }
    }
    cout << "Введіть замість q,w,e,r (qx1+wx2/ex1+rx2) коефіцієнти ф-ї мети" << endl;
    for (i = 2; i < zm + 2; i++) {
        cin >> pmat[0][i];
    }
    cin >> pmat[rivn][2] >> pmat[rivn][3];
    pmat[rivn][rivn + zm + maxSize - 2] = 1;
    pmat[rivn][0] = -M;
    pmat[0][rivn + zm + maxSize - 3] = -M;
    cout << "Виберіть тип задачі(1-на максимум, 2-на мінімум)\n";
    cin >> type;
    cout << "Отримані результати:" << endl;
    for (i = 1, j = 4; i < rivn + 1; i++, j++) {
        pmat[i][j] = 1;
    }
    for (i = 2; i < zm + 2; i++)
        if (type == 2)
            pmat[0][i] *= -1;
    FillZ();
    Print();
    while (1) {
        imin = MinZ();
        if (!imin) {
            cout << "Цей план є оптимальним\n";
            break;
        }
        if (kilk == 2 && flag == 0) {
            for (i = 0; i < rivn + 2; i++)
                for (j = 0; j < rivn + zm + maxSize; j++)
                    pmat2[i][j] = pmat[i][j];
            flag = 1;
        }
        jmin = FindRiad(imin);
        if (!jmin) {
            cout << "Неможливо застосувати правило В\n" << endl;
            cin.get(); // Pause program
            return 0;
        }
        nmin = pmat[jmin][imin];
        for (i = 1; i < rivn + zm + maxSize - 1; i++)
            pmat[jmin][i] /= nmin;
        cout << "\nВибрано ключовий елемент (" << jmin << ',' << imin << ')';
        zmina(imin, jmin);
        FillZ();
        Print();
    }
    for (i = 1; i < rivn + 1; i++) {
        if (pmat[i][0] != 0)
            for (j = 2; j < zm + 2; j++)
                if (pmat[i][j] == 1) {
                    if (CheckSt(j)) {
                        if (j == zm)
                            xmax = pmat[i][rivn + zm + maxSize - 2];
                        if (j == zm + 1)
                            ymax = pmat[i][rivn + zm + maxSize - 2];
                    }
                }
    }
    for (i = 1; i < rivn + 1; i++) {
        if (pmat[i][1] == 1)
            yg1 = pmat[i][rivn + zm + maxSize - 2];
        else if (pmat[i][1] != 0)
            return 0;
    }
    ymax /= yg1;
    xmax /= yg1;
    if (type == 1) {
        cout << "Максимум";
    }
    else {
        cout << "Мінімум";
    }
    cout << " ф-ї знаходиться в точці (" << xmax << ',' << ymax << ')' << "=";
    if (type == 1) {
        cout << pmat[rivn + 1][rivn + zm + maxSize - 2] << endl;
    }
    else {
        cout << -pmat[rivn + 1][rivn + zm + maxSize - 2] << endl;
    }
    cin.get(); // Pause program
    // Free allocated memory
    for (i = 0; i < rivn + 2; i++) {
        free(pmat[i]);
        free(pmat2[i]);
    }
    free(pmat);
    free(pmat2);
    return 0;
}
void Print()
{
    int i, l, size2 = maxSize;
    cout << endl << endl << '\t';
    for (i = 1; i < rivn + zm + 1; i++) {
        cout << pmat[0][i] << '\t';
    }
    if (pmat[0][i] == -M) {
        cout << "-M" << '\t';
    }
    cout << endl << "Ci\t";
    for (i = 0; i < rivn + zm; i++) {
        cout << "p" << i + 1 << '\t';
    }
    while (size2-- > 3) {
        cout << "z" << ++i << '\t';
    }
    cout << 'b' << endl;
    for (int i = 1; i < rivn + 1; i++) {
        for (int j = 0; j < rivn + zm + maxSize - 1; j++) {
            if (j == 0 && pmat[i][j] == M) {
                cout << "M\t";
            }
            else {
                if (j == 0 && pmat[i][j] == -M) {
                    cout << "-M\t";
                }
                else {
                    cout << pmat[i][j] << "\t";
                }
            }
        }
        cout << endl;
    }
    cout << "Z-C\t";
    for (int j = 1; j < rivn + zm + maxSize - 1; j++) {
        if (pmat[rivn + 1][j] > M / 1000) {
            for (l = 1; l < rivn + 1; l++)
                if (pmat[l][0] == M || pmat[l][0] == -M)
                    break;
            cout << "M";
            if (l < rivn + 1) {
                cout << '*' << pmat[rivn + 1][j] / M << '\t';
            }
            else {
                cout << '\t';
            }
        }
        else {
            if (pmat[rivn + 1][j] < (-M / 1000)) {
                for (l = 1; l < rivn + 1; l++)
                    if (pmat[l][0] == M || pmat[l][0] == -M)
                        break;
                cout << "-M";
                if (l < rivn + 1) {
                    cout << '*' << pmat[rivn + 1][j] / -M << '\t';
                }
                else {
                    cout << '\t';
                }
            }
            else {
                cout << pmat[rivn + 1][j] << "\t";
            }
        }
    }
    cout << endl;
}
void FillZ()
{
    int i, j;
    double znach = 0;
    for (i = 1; i < rivn + zm + maxSize - 1; i++) {
        for (j = 1; j < rivn + 1; j++)
            znach += pmat[j][0] * pmat[j][i];
        pmat[rivn + 1][i] = znach - pmat[0][i];
        znach = 0;
    }
}
bool CheckSt(int i)
{
    int k, kst = 0;
    for (k = 1; k < rivn + 1; k++) {
        if (pmat[k][i] == 1) {
            kst++;
            continue;
        }
        if (pmat[k][i] != 0)
            return false;
    }
    if (kst < 2)
        return true;
    return false;
}
int MinZ()
{
    int imin = 0;
    double min = 0;
    for (int i = 1; i < rivn + zm + maxSize - 2; i++)
        if (pmat[rivn + 1][i] < min) {
            imin = i;
            min = pmat[rivn + 1][i];
        }
    for (int i = 1; i < rivn + zm + maxSize - 2; i++) {
        if (i == imin)
            continue;
        if (pmat[rivn + 1][i] == min && min != 0) {
            imin2 = i;
            kilk = 2;
            break;
        }
    }
    return imin;
}
int MaxZ()
{
    int imax = 0;
    double max = 0;
    for (int i = 1; i < rivn + zm + maxSize - 2; i++)
        if (pmat[rivn + 1][i] > max) {
            imax = i;
            max = pmat[rivn + 1][i];
        }
    for (int i = 1; i < rivn + zm + maxSize - 2; i++) {
        if (i == imax)
            continue;
        if (pmat[rivn + 1][i] == max && max != 0) {
            imin2 = i;
            kilk = 2;
            break;
        }
    }
    return imax;
}
int FindRiad(int st) {
    int imin = 0;
    double min = 100000000;
    for (int i = 1; i < rivn + 1; i++)
        if (pmat[i][rivn + zm + maxSize - 2] / pmat[i][st] < min && pmat[i][rivn + zm + maxSize - 2] / pmat[i][st] >= 0 && pmat[i][st] > 0) {
            imin = i;
            min = pmat[i][rivn + zm + maxSize - 2] / pmat[i][st];
        }
    return imin;
}
void zmina(int xmin, int ymin) {
    int x, y;
    double mnozh;
    for (y = 1; y < rivn + 1; y++) {
        if (y == ymin)
            continue;
        mnozh = -pmat[y][xmin];
        for (x = 1; x < rivn + zm + maxSize - 1; x++) {
            pmat[y][x] = pmat[y][x] + pmat[ymin][x] * mnozh;
        }
    }
    pmat[ymin][0] = pmat[0][xmin];
}
void Znaku()
{
    int i;
    for (i = 1; i < rivn + 1; i++)
        if (pmat[i][rivn + zm + maxSize - 1] == 2) {
            for (int j = 1; j < rivn + zm + 2; j++)
                if (pmat[i][j])
                    pmat[i][j] *= -1;
        }
}