#include "FFT.h"

int main()
{
    FastFourierTransform F1;
    std::vector<std::complex<double>> signal {{1, 4}, {5,10}, {11, 10}, {1, 1}, {1.0, 0.0}}; //, {11, 10}, {1, 1}, {1.0, 0.0}
    std::vector<std::complex<double>> FFTSignal = F1.DFT5(signal, 0);
    std::vector<std::complex<double>> FFTSignal2 = F1.DFT5(FFTSignal, 1);
    for (int i = 0; i<5; i++)
    {
        std::cout << "Re: " << FFTSignal[i].real() << " Im: " << FFTSignal[i].imag() << std::endl;
    }
    std::cout <<"______________________________________________"<< std::endl;
    for (int i = 0; i<5; i++)
    {
        std::cout << "Re: " << FFTSignal2[i].real() << " Im: " << FFTSignal2[i].imag() << std::endl;
    }
    return 0;
}