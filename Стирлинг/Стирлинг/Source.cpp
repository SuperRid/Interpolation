#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>


using namespace std;

double f(double x)
{
	return tan(x);
}

void KonRaz(double* y, int n, double** KR)
{
	for (int i = -n / 2; i <= n / 2; i++)
	{
		KR[i][-n/2] = y[i];
	}

	for (int j = (-n / 2) + 1; j <= n / 2; j++)
	{
		for (int i = -n/2; i <= -j; i++)
		{
			KR[i][j] = KR[i + 1][j - 1] - KR[i][j - 1];
		}
	} 
} 

double Stirling(double* x_, double* y, int n, double& x, double h,double** KR)
{
	double pq1 = (x - x_[0]) / h, pq2 = pow((x - x_[0]) / h, 2) / 2, rq1 = 0, rq2 = 0, P = 0, q = 0;

	KonRaz(y, n, KR);
	
	for (int i = 2; i <= n/2; i++)
	{
		q = pow(((x - x_[0]) / h),2) - pow((i - 1),2);
		pq1 *= q / ((2+2*(i-2)) * (3 + 2 * (i - 2)));
		rq1 += pq1 * (KR[-i][((-n/2)+1) + 2*(i-1)] + KR[-i + 1][((-n / 2) + 1) + 2 * (i - 1)]) / 2;
		pq2 *= q /((3 + 2 * (i - 2)) * (4 + 2 * (i - 2)));
		rq2 += pq2 * KR[-i][((-n / 2) + 2) + 2 * (i - 1)];
		
		
	}
	P = y[0] + ((x-x_[0]) / h)* ((KR[-1][(-n / 2) + 1] + KR[0][(-n / 2) + 1]) / 2)  + rq1 + ((pow((x - x_[0]) / h,2)/2) * KR[-1][(-n / 2) + 2]) + rq2;
	return P;
}



int main()
{
	setlocale(LC_ALL, "ru");

	double a = 0, b = 1.5;
	int n;

	cout <<"Введите n, n - четное" << endl << "n = ";
	cin >> n;

	while (n % 2 != 0)
	{
		cout << "Вы ввели нечетное n, введите четное n" << endl << "n = ";
		cin >> n;
	}
	cout << endl; 

	
	double* x = new double[n + 1];
	double* y = new double[n + 1];
	double** KR = new double* [n + 1];
	for (int i = -n / 2; i <= n / 2; i++)
	{
		KR[i] = new double[n + 1];
	}


	cout << "Таблица значений функции yi = f(xi)" << endl;
	cout << " " << endl;
	cout << setw(5) << "x" << " | " << setw(10) << "f(x)" << endl;
	double h = (b - a) / n;
	for (int i = -n/2; i <= n/2; i++)
	{
		x[i] = 0.5*(b-a) + i * h;
		y[i] = f(x[i]);
		cout << setw(5) << x[i] << " | " << setw(10) << y[i] << " | " << endl;

	}
	cout << " " << endl;


	cout << "Таблица с реализованным методом " << endl;
	cout << " " << endl;
	cout << setw(5) << "xj" << " | " << setw(10) << "f(xj)" << " | " << setw(15) << "F(xj)" << " | " << setw(20) << "f(xj) - F(xj)" << endl;
	for (double xj = 0; xj < 1.55; xj += 0.05)
	{
		double fi = f(xj);
		double F = Stirling(x, y, n, xj, h, KR);
		cout << setw(5) << xj << " | " << setw(10) << fi << " | " << setw(15) << F << " | " << setw(20) << fi - F << endl;
	}

	cin >> n;
	for (int i = -n / 2; i <= n / 2; i++)
	{
		delete[] KR[i];
	}

	delete[] KR;
	delete[] x;
	delete[] y;

	return 0;
}
