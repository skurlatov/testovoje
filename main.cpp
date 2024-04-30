#include "FFT.h"

int main()
{
    FastFourierTransform F1;
    std::vector<std::complex<double>> signal(300);
    std::vector<std::complex<double>> transformedSignal;
    std::vector<std::complex<double>> error(signal);
    for(int i = 0; i<signal.size(); i++)
        signal[i] = {rand()%10, rand()%10};
    transformedSignal = F1.Trasform(signal, 0);
    transformedSignal = F1.Trasform(transformedSignal, 1);
    for(int i = 0; i<signal.size(); i++)
    {
        error[i] = signal[i] - transformedSignal[i];
        std::cout << "Re: " << error[i].real() << " Im: " << error[i].imag() << std::endl;
    }
    return 0;
}