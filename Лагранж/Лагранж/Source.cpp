#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>


using namespace std;

double f(double x) 
{
	return sin(x*x*x);
}

double Lagr(double* x_, double* y, double n, double& x )
{
	double L = 0;

	for (int i = 0; i <= n; i++)
	{
		double p = 1;
		for (int j = 0; j <= n; j++)
		{
			if(j!=i) p *= (x - x_[j]) / (x_[i] - x_[j]);
		}
		L += y[i] * p;

	}
	return L;
}



int main()
{
	setlocale(LC_ALL, "ru");
	
	double a = 0, b = 2;
	int n;
	cout <<" n = ";
	cin >> n;
	cout << endl;

	double x, g;
	double* x_ = new double[n + 1];
	double* y = new double[n + 1];

	cout << "Таблица значений функции yi = f(xi)" << endl;
	cout << " " << endl;


	
	cout << setw(5) << "x" << " | " << setw(10) << "f(x)" << endl;
	for (int i = 0; i <= n; i++)
	{
		double h = (b - a) / n;
		x_[i] = a + i * h;
		y[i] = f(x_[i]);
		cout << setw(5) << x_[i]<< " | " << setw(10) << y[i] << " | " << endl;
		
	}
	cout << " " << endl;

	cout << "Таблица с реализованным методом " << endl;
	cout << " " << endl;
	cout << setw(5) << "xj" << " | " << setw(10) << "f(xj)" << " | " << setw(15) << "F(xj)" << " | " << setw(20) << "f(xj) - F(xj)" << endl;
	for (double xj = 0; xj < 2.08; xj += 0.08)
	{
		double fi = f(xj);
		double F = Lagr(x_, y, n, xj);

		cout << setw(5) << xj << " | " << setw(10) << fi << " | " << setw(15) << F << " | " << setw(20) << fi - F << endl;
	}
	
	cout << " "  << endl;
	cout << "Введите нужный угол в градусах: " << endl;
	cin >> g;
	x = g * M_PI / 180;
	cout << Lagr(x_, y, n, x)  << endl;
	
	delete[] x_;
	delete[] y;


}
