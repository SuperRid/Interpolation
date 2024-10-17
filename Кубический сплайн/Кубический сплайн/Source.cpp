#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>


using namespace std;

double f(double x)
{
	return 1/(x*x);
}

double Splain (double* x_, double* y, int n, double h, double& x)
{
	double S = 0, half1 = 0, half2 = 0, half3 = 0, half4 = 0;
	double* m = new double[n + 1];

	for (int i = 0; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			m[0] = (4 * y[1] - y[2] - 3 * y[0]) / (2 * h);
			m[n] = (3 * y[n] + y[n - 2] - 4 * y[n - 1]) / (2 * h);
			m[j] = (y[j + 1] - y[j - 1]) / (2 * h);
		}

		half1 = (pow(x_[i + 1] - x, 2) * (2 * (x - x_[i]) + h) * y[i]) / pow(h, 3);
		half2 = (pow(x - x_[i], 2) * (2 * (x_[i + 1] - x) + h) * y[i + 1]) / pow(h, 3);
		half3 = (pow(x_[i + 1] - x, 2) * (x - x_[i]) * m[i]) / pow(h, 2);
		half4 = (pow(x - x_[i], 2) * (x - x_[i + 1]) * m[i + 1]) / pow(h, 2);

		S = half1 + half2 + half3 + half4;
	}
	delete[] m;
	return S;
}


int main()
{
	setlocale(LC_ALL, "ru");

	double a = 0.1, b = 1;
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
	for (double xj = 0.1; xj < 1.2; xj += 0.05)
	{
		double fi = f(xj);
		double F = Splain(x_, y, n, h, xj);

		cout << setw(5) << xj << " | " << setw(10) << fi << " | " << setw(15) << F << " | " << setw(20) << fi - F << endl;
	}

	cout << " " << endl;
	cout << "Введите нужный угол в градусах: " << endl;
	cin >> g;
	x = g * M_PI / 180;
	cout << Splain(x_, y, n, h, x) << endl;

	delete[] x_;
	delete[] y;


}
