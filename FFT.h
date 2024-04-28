#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

class FastFourierTransform
{
    private:
    //Уникальные значения матриц ДПФ
    std::vector<std::complex<double>> DFT3Matrix = {{1.0, 0.0}, {-0.5, -0.866}, {-0.5, 0.866}}; //Уникальные значения матрицы ДПФ N = 3
    std::vector<std::complex<double>> DFT5Matrix = {{1.0, 0.0}, {0.309, -0.951},{-0.809, -0.588},{-0.809, 0.588},{0.309, 0.951}}; //Уникальные значения матрицы ДПФ N = 5
    
    public:
    std::vector<std::complex<double>> DFT2(std::vector<std::complex<double>> signal, bool flag); //Функция ДПФ для N = 2, flag указывает прямое (0) или обратное (1) преобразование
    std::vector<std::complex<double>> DFT3(std::vector<std::complex<double>> signal, bool flag); //Функция ДПФ для N = 3
    std::vector<std::complex<double>> DFT5(std::vector<std::complex<double>> signal, bool flag); //Функция ДПФ для N = 5
    std::vector<std::complex<double>> Trasform(std::vector<std::complex<double>> signal, bool flag); //Функция БПФ для N кратного 2, 3 и 5
};