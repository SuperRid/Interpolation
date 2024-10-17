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
	// Инициализируем массив разделенных разностей с нулевыми значениями
	for (int i = 0; i <= n; i++)
	{
		KR[i][0] = y[i];
	}

	// Вычисляем разделенные разности
	for (int j = 1; j <= n; j++)
	{
		for (int i = 0; i <= n - j; i++)
		{
			KR[i][j] = KR[i + 1][j - 1] - KR[i][j - 1];
		}
	}
}

double Newton(double* x_, double* y, int n, double& x,double h)
{
	double pq = 1, rq = 0, P = 0,q = 0;
	

	double** KR = new double* [n + 1];
	for (int i = 0; i <= n; i++)
	{
		KR[i] = new double[n + 1];
	}

	KonRaz(y, n, KR);

	for (int i = 1; i <= n; i++)
	{
		q = ((x - x_[0]) / h) - (i-1) ;
		pq *= q / i;
		rq += pq * KR[0][i];
	}
	P = y[0] + rq;
	
	for (int i = 0; i <= n; i++)
	{
		delete[] KR[i];
	}
	
	delete[] KR;
	return P;
}



int main()
{
	setlocale(LC_ALL, "ru");

	double a = 0.5, b = 1.5;
	int n;
	cout << " n = ";
	cin >> n;
	cout << endl;

	double x, g;
	double* x_ = new double[n + 1];
	double* y = new double[n + 1];
	

	cout << "Таблица значений функции yi = f(xi)" << endl;
	cout << " " << endl;
	cout << setw(5) << "x" << " | " << setw(10) << "f(x)" << endl;

	double h = (b - a) / n;
	for (int i = 0; i <= n; i++)
	{
		x_[i] = a + i * h;
		y[i] = f(x_[i]);
		cout << setw(5) << x_[i] << " | " << setw(10) << y[i] << " | " << endl;

	}
	cout << " " << endl;


	cout << "Таблица с реализованным методом " << endl;
	cout << " " << endl;
	cout << setw(5) << "xj" << " | " << setw(10) << "f(xj)" << " | " << setw(15) << "F(xj)" << " | " << setw(20) << "f(xj) - F(xj)" << endl;
	for (double xj = 0.2; xj <1; xj += 0.1)
	{
		double fi = f(xj);
		double F = Newton(x_, y, n, xj,h);

		cout << setw(5) << xj << " | " << setw(10) << fi << " | " << setw(15) << F << " | " << setw(20) << fi - F << endl;
	}

	cout << " " << endl;
	cout << "Введите нужный угол в градусах: " << endl;
	cin >> g;
	x = g * M_PI / 180;
	cout << Newton(x_, y, n, x,h) << endl;


	delete[] x_;
	delete[] y;


	
} 
