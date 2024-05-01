#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>
#define PI 3.141592654

class FastFourierTransform
{
    private:
    //Уникальные значения матриц ДПФ
    std::vector<std::complex<double>> DFT3Matrix = {{1.0, 0.0}, {-0.5, -0.866}, {-0.5, 0.866}}; //Уникальные значения матрицы ДПФ N = 3
    std::vector<std::complex<double>> DFT5Matrix = {{1.0, 0.0}, {0.309, -0.951},{-0.809, -0.588},{-0.809, 0.588},{0.309, 0.951}}; //Уникальные значения матрицы ДПФ N = 5
        std::vector<std::complex<double>> DFT2(std::vector<std::complex<double>> signal, bool flag)
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
    } //Функция ДПФ для N = 2, flag указывает прямое (0) или обратное (1) преобразование
    std::vector<std::complex<double>> DFT3(std::vector<std::complex<double>> signal, bool flag){
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
    } //Функция ДПФ для N = 3
    std::vector<std::complex<double>> DFT5(std::vector<std::complex<double>> signal, bool flag){
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
    } //Функция ДПФ для N = 5
    
    
    public:
    std::vector<std::complex<double>> Trasform(std::vector<std::complex<double>> signal, bool flag)
    {
        int N = signal.size();
        std::vector<std::complex<double>> copySignal = signal;
        int marker = 0;
        if (N % 5 == 0 && marker == 0)
        {
        marker++;
        std::vector<std::vector<std::complex<double>>> transformed(N/5);
        for (int i = 0; i<N/5; i++) //проход по числу разбиений
        {
            for (int j = 0; j<N; j++) //проход по всем элементам сигнала
            {
                if(j%(N/5) == i)
                transformed[i].push_back(copySignal[j]);
            }
            transformed[i] = DFT5(transformed[i], flag);
        }
        if (N/5 != 1)
        {
            std::vector<std::vector<std::complex<double>>> transformedRec(5);
            for(int i = 0; i < 5; i++)
            {
                for(int j = 0; j < N/5; j++)
                transformedRec[i].push_back(transformed[j][i]*std::polar(1.0, (-1+2*int(flag))*2*PI*j*i/N));
                transformedRec[i] = Trasform(transformedRec[i], flag);
            }
            
            for(int i = 0; i < 5; i++)
            {
                for(int j = 0; j < N/5; j++)
                transformed[j][i] = transformedRec[i][j];
            }
            
        }
        std::vector<std::complex<double>> TEMP;
        for(int i = 0; i < N/5; i++)
            for(int j = 0; j < 5; j++)
                TEMP.push_back(transformed[i][j]);
        copySignal = TEMP;
        return copySignal;
        }

        else if (N % 3 == 0 && marker == 0)
        {
        marker++;
        std::vector<std::vector<std::complex<double>>> transformed(N/3);
        for (int i = 0; i<N/3; i++) //проход по числу разбиений
        {
            for (int j = 0; j<N; j++) //проход по всем элементам сигнала
            {
                if(j%(N/3) == i)
                transformed[i].push_back(copySignal[j]);
            }
            transformed[i] = DFT3(transformed[i], flag);
        }
        if (N/3 != 1)
        {
            std::vector<std::vector<std::complex<double>>> transformedRec(3);
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < N/3; j++)
                transformedRec[i].push_back(transformed[j][i]*std::polar(1.0, (-1+2*int(flag))*2*PI*j*i/N));
                transformedRec[i] = Trasform(transformedRec[i], flag);
            }
            
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < N/3; j++)
                transformed[j][i] = transformedRec[i][j];
            }
            
        }
        std::vector<std::complex<double>> TEMP;
        for(int i = 0; i < N/3; i++)
            for(int j = 0; j < 3; j++)
                TEMP.push_back(transformed[i][j]);
        copySignal = TEMP;
        return copySignal;
    }
    else if (N % 2 == 0 && marker == 0)
    {
        marker++;
        std::vector<std::vector<std::complex<double>>> transformed(N/2);
        for (int i = 0; i<N/2; i++) //проход по числу разбиений
        {
            for (int j = 0; j<N; j++) //проход по всем элементам сигнала
            {
                if(j%(N/2) == i)
                transformed[i].push_back(copySignal[j]);
            }
            transformed[i] = DFT2(transformed[i], flag);
        }
        if (N/2 != 1)
        {
            std::vector<std::vector<std::complex<double>>> transformedRec(2);
            for(int i = 0; i < N/2; i++)
            {
                for(int j = 0; j < 2; j++)
                transformed[i][j] = transformed[i][j]*std::polar(1.0, (-1+2*int(flag))*2*PI*j*i/N);
            }
            for(int i = 0; i < 2; i++)
            {
                for(int j = 0; j < N/2; j++)
                transformedRec[i].push_back(transformed[j][i]);
                transformedRec[i] = Trasform(transformedRec[i], flag);
            }
            for(int i = 0; i < 2; i++)
            {
                for(int j = 0; j < N/2; j++)
                transformed[j][i] = transformedRec[i][j];
            }
            
        }
        std::vector<std::complex<double>> TEMP;
        for(int i = 0; i < N/2; i++)
            for(int j = 0; j < 2; j++)
                TEMP.push_back(transformed[i][j]);
        copySignal = TEMP;
        return copySignal;
    }
    return {(0,0)};
    }


    }; //Функция БПФ для N кратного 2, 3 и 5




int main()
{
    FastFourierTransform F1;
    std::vector<std::complex<double>> signal{{0, 0}, {0, 1}, {1,0}, {1, 1}};
    std::vector<std::complex<double>> FFTS = F1.Trasform(signal, 0);
    for (int i = 0; i < signal.size(); i++)
        std::cout << '#' << i+1 << " Re: " << FFTS[i].real() << " Im: " << FFTS[i].imag()<<std::endl;
    std::cout << "__________________________________"<<std::endl;
    FFTS = F1.Trasform(FFTS, 1);
    for (int i = 0; i < signal.size(); i++)
        std::cout << '#' << i+1 << " Re: " << FFTS[i].real() << " Im: " << FFTS[i].imag()<<std::endl;
    return 0;
}