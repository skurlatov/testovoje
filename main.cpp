#include "FFT.h"

int main()
{
    FastFourierTransform F1;
    std::vector<std::complex<double>> signal(300); //длина может быть любой, но кратной 2, 3 или 5
    std::vector<std::complex<double>> transformedSignal; //трансформированная последовательность
    std::vector<std::complex<double>> error(signal); //вектор ошибки
    for(int i = 0; i<signal.size(); i++)
        signal[i] = {rand()%10, rand()%10};
    transformedSignal = F1.Trasform(signal, 0); //прямое БПФ
    transformedSignal = F1.Trasform(transformedSignal, 1); //обратное БПФ
    for(int i = 0; i<signal.size(); i++)
    {
        error[i] = signal[i] - transformedSignal[i]; //формирование вектора ошибки
        std::cout << "Re: " << error[i].real() << " Im: " << error[i].imag() << std::endl; //вывод значений вектора ошибки
    }
    return 0;
}