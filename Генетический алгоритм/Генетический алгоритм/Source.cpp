#define _USE_MATH_DEFINES
#define slon
#define population_size 100
#ifdef slon

#include <iostream>
#include <cmath>
#include <ctime>
//#include "Vector.h"
#include <valarray>

using namespace std;

double rosenbrock_func(valarray<double> arg)
{
	if (arg.size() == 2)
	{
		double x = arg[0];
		double y = arg[1];
		return pow((1 - x), 2) + 100 * pow((y - x * x), 2);
	}
	else
		throw - 1;
}

double j_func(valarray<double> arg)
{
	if (arg.size() == 2)
	{
		double x = arg[0];
		double y = arg[1];
		return (1+sin(10*x)+cos(2*x)+cos(2*x+2*y)+cos(2*y)+sin(20*y)+y*y);
	}
	else
		throw - 1;
}

valarray<double> mutation(const valarray<double> &x0,const valarray<double> &x1)  // мутация: генерация случайной величины
{
	valarray<double> child(2);
	const int NUM= 100000000;	
	child[0] = fabs((double)((NUM * rand()) % (int)((x1[0] - x0[0]) * NUM) + 1) / NUM) + x0[0];
	child[1] = fabs((double)((NUM * rand()) % (int)((x1[1] - x0[1]) * NUM) + 1) / NUM) + x0[1];
	return child;
}

double mutation_double(double x0, double x1)  // мутация: генерация случайной величины
{
	double child;
	const int NUM = 100000000;
	child = fabs((double)((NUM * rand()) % (int)((x1 - x0) * NUM) + 1) / NUM) + x0;	
	return child;
}

valarray<double> inversion(const valarray<double> &x, double eps)  // инверсия: поиск в окрестностях точки
{
	valarray<double> x_invers;
	static int sign = 0;
	sign++;
	sign %= 2;
	for (int i = 0; i < 2; i++)
	{
		if (sign == 0)
			 x_invers[i]=x[i] - eps;
		else
			x_invers[i] = x[i] + eps;
	}	
	return x_invers;
}

void crossover(valarray<double> *population, double eps,const valarray <double> x0, const valarray <double> x1)  // кроссовер: среднее арифметическое
{
	int k = 99;
	for (int i = 0; i < 8; i++)
		for (int j = i + 1; j < 8; j++)
		{
			population[k] = (population[i] + population[j]) * 0.5;
			k--;
		}
	for (int i = 0; i < 8; i++)
	{
		population[k] = inversion(population[i], eps); k--;
		population[k] = inversion(population[i], eps); k--;
	}
	for (int i = 8; i < k; i++)
		population[i] = mutation(x0, x1);
}

void crossover_population(valarray<double> *population, double eps,const valarray<double> &x0,const valarray<double> &x1)
{
	for (int i = 0; i < 2; i++)
	{
		crossover(population, eps, x0[i], x1[i]);
	}
}

void sort(valarray<double> *x, double *y)  // сортировка
{
	for (int i = 0; i < population_size; i++)
		for (int j = i + 1; j < population_size; j++)	
			if (y[j] < y[i])
			{
				double temp = y[i];
				valarray<double> temp_vector=x[i];
				y[i] = y[j];
				y[j] = temp;
				x[i] = x[j];
				x[j] = temp_vector;
			}
}

valarray<double> genetic(const valarray<double> &x0,const valarray<double> &x1, double eps)  // поиск решения с использованием ГА
{
	valarray <double> *population = new valarray<double>[population_size];
	double f[population_size];
	int iter = 0;
	for (int i = 0; i < population_size; i++)   // Формирование начальной популяции
	{
		population[i] = mutation(x0, x1);
		f[i] = rosenbrock_func(population[i]);
	}
	sort(population, f);
	do {
		iter++;
		crossover_population(population, eps, x0, x1);
		for (int i = 0; i < population_size; i++)
			f[i] = rosenbrock_func(population[i]);
		sort(population, f);
	} while (fabs(f[0]) > eps && iter < 200000);
	cout << iter << " iterations" << endl;
	return population[0];
}
int main()
{
	srand(time(NULL));
	valarray<double> individ_1(2);	
	valarray<double> individ_2(2);
	for (int i = 0; i < 2; i++)
	{
		individ_1[i] = 3+i;
		individ_2[i] = 2+i;
	}
	cout << genetic(individ_1, individ_2, 0.000001);
	cin.get();
	return 0;
}

#else //RRRRRRRRRRRRRRRRRRRRRRRRRR
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;
double func(double x)
{
	return sin(M_PI * x / 180) + 1 / x;
}
double mutation(double x0, double x1)  // мутация: генерация случайной величины
{
	const int NUM = 100000000;
	return fabs((double)((rand() * NUM) % (int)((x1 - x0)*NUM) + 1) / NUM) + x0;
}
double inversion(double x, double eps)  // инверсия: поиск в окрестностях точки
{
	static int sign = 0;
	sign++;
	sign %= 2;
	if (sign == 0) return x - eps;
	else return x + eps;
}
void crossover(double *x, double eps, double x0, double x1)  // кроссовер: среднее арифметическое
{
	int k = 99;
	for (int i = 0; i < 8; i++)
		for (int j = i + 1; j < 8; j++)
		{
			x[k] = (x[i] + x[j]) / 2;
			k--;
		}
	for (int i = 0; i < 8; i++)
	{
		x[k] = inversion(x[i], eps); k--;
		x[k] = inversion(x[i], eps); k--;
	}
	for (int i = 8; i < k; i++)
		x[i] = mutation(x0, x1);
}
void sort(double *x, double *y)  // сортировка
{
	for (int i = 0; i < 100; i++)
		for (int j = i + 1; j < 100; j++)
			if (fabs(y[j]) < fabs(y[i])) {
				double temp = y[i];
				y[i] = y[j];
				y[j] = temp;
				temp = x[i];
				x[i] = x[j];
				x[j] = temp;
			}
}
double genetic(double x0, double x1, double eps)  // поиск решения с использованием ГА
{
	double population[100];
	double f[100];
	int iter = 0;
	for (int i = 0; i < 100; i++)   // Формирование начальной популяции
	{
		population[i] = mutation(x0, x1);
		f[i] = func(population[i]);
	}
	sort(population, f);
	do {
		iter++;
		crossover(population, eps, x0, x1);
		for (int i = 0; i < 100; i++)
			f[i] = func(population[i]);
		sort(population, f);
	} while (fabs(f[0]) > eps && iter < 20000);
	cout << iter << " iterations" << endl;
	return population[0];
}
int main()
{
	srand(time(NULL));
	cout << genetic(1.0, 10.0, 0.000001);
	cin.get();
	return 0;
}
#endif // slon