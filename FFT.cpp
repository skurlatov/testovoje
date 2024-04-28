#include "FFT.h"
#pragma once

std::vector<std::complex<double>> FastFourierTransform::DFT2(std::vector<std::complex<double>> signal, bool flag) //Функция ДПФ-2
{
    std::vector<std::complex<double>>transformed(2); //Инициализация Фурье-образа изначальной последовательности
    transformed[0] = signal[0] + signal[1]; //ДПФ-2 выполняется простыми арифметическими действиями 
    transformed[1] = signal[0] - signal[1];
    if (flag == 0)
        return {transformed};
    else
        {
            for (int i = 0; i<2; i++)
                transformed[i] = 0.5*transformed[i]; //Обратное ДПФ-2
            return {transformed};
        }

}
std::vector<std::complex<double>> FastFourierTransform::DFT3(std::vector<std::complex<double>> signal, bool flag)//Функция ДПФ-3
{
    std::vector<std::complex<double>>transformed(3);
    transformed[0] = signal[0]*DFT3Matrix[0] + signal[1]*DFT3Matrix[0] + signal[2]*DFT3Matrix[0]; //ДПФ-3 также можно представить в виде небольшого количества операций
    transformed[1] = signal[0]*DFT3Matrix[0] + signal[1]*DFT3Matrix[1+int(flag)] + signal[2]*DFT3Matrix[2-int(flag)]; //Значение flag влияет на индексацию таким образом,
    transformed[2] = signal[0]*DFT3Matrix[0] + signal[1]*DFT3Matrix[2-int(flag)] + signal[2]*DFT3Matrix[1+int(flag)]; //чтобы получилось обратное преобразование (без множителя 1/N) при flag=1 и прямое при flag=0
    if (flag == 0)
        return {transformed};
    else
        {
            for (int i = 0; i<3; i++)
                transformed[i] = 0.333*transformed[i];
            return {transformed};
        }
}
std::vector<std::complex<double>> FastFourierTransform::DFT5(std::vector<std::complex<double>> signal, bool flag)//Функция ДПФ-5
{
    std::vector<std::complex<double>> transformed(5);
    transformed[0] = signal[0] + signal[1] + signal[2] + signal[3] + signal[4]; //Нулевой отчет формируется простым суммированием
    for(int i = 1; i<5; i++)
    {
        transformed[i] = signal[0] + signal[1]*DFT5Matrix[(i*(1+int(flag)*3))%5] + signal[2]*DFT5Matrix[(i*2*(1+int(flag)*3))%5] + signal[3]*DFT5Matrix[(i*3*(1+int(flag)*3))%5] + signal[4]*DFT5Matrix[(i*4*(1+int(flag)*3))%5];
    }
    if (flag == 0)
        return {transformed};
     else
        {
            for (int i = 0; i<5; i++)
                transformed[i] = 0.2*transformed[i];
            return {transformed};
        }
}
std::vector<std::complex<double>> FastFourierTransform::Trasform(std::vector<std::complex<double>> signal, bool flag)
{


    return {(0,0)};

}