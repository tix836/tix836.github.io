#define _CRT_SECURE_NO_WARNINGS
#include "stdio.h"
#include<iostream>
#include<fstream>
#include<vector>
#include "LOS.h"


using namespace std;
// 192 общие формулы
// 67 на шаблонных базисных функциях
void read_grid(string grid_file, int& grid_size, vector<double>& grid )
{
	fstream file_input;

	file_input.open(grid_file);
	if (file_input.is_open())
	{
		file_input >> grid_size;
		grid.resize(grid_size);
		for (int i = 0; i < grid_size; i++)
		{
			file_input >> grid[i];
		}
	}
	else
	{		
		perror("error");
	}
	file_input.close();
}

void read_params(double& chi, double& sigma, double& lamda)
{
	fstream file_input;
	file_input.open("params.txt");
	if (file_input.is_open())
	{
		file_input >> chi;
		file_input >> sigma;
		file_input >> lamda;
	}
	file_input.close();
}

double F_function(double x, double t)
{
	
	return x;
}
 
double U_function(double x, double t)
{
	
	return x * t;
}

vector<vector<double>> Local_M(double K, double x0, double x1)
{
	vector<vector<double>> local_matrix_M =
	{
		{2,1},
		{1,2}
	};
	double h = x1 - x0;

	for (int i = 0; i < local_matrix_M.size(); i++)
	{
		for (int j = 0; j < local_matrix_M.size(); j++)
		{
			local_matrix_M[i][j] = ((K * h) / 6)* local_matrix_M[i][j];
		}
	}

	return local_matrix_M;
}

vector<vector<double>> Global_M(vector<double> grid_space, double K)
{
	int grid_space_size = grid_space.size();
	int num_of_fin_elems = grid_space.size() - 1;

	vector<vector<double>> global_M;

	global_M.resize(grid_space_size);
	for (int i = 0; i < grid_space_size; i++)
	{
		global_M[i].resize(grid_space_size);
	}

	for (int i = 0; i < num_of_fin_elems; i++)
	{
		vector<vector<double>> local_M = Local_M(K, grid_space[i], grid_space[i + 1]);
		global_M[i][i] += local_M[0][0];
		global_M[i][i + 1] += local_M[0][1];
		global_M[i + 1][i] += local_M[1][0];
		global_M[i + 1][i + 1] += local_M[1][1];
	}
	return global_M;
}

vector<vector<double>> Local_G(double lamda, double x0, double x1)
{
	vector<vector<double>> local_matrix_G =
	{
		{1, -1},
		{-1, 1}
	};
	double h = x1 - x0;

	for (int i = 0; i < local_matrix_G.size(); i++)
	{
		for (int j = 0; j < local_matrix_G.size(); j++)
		{
			local_matrix_G[i][j] = (lamda / h) * local_matrix_G[i][j];
		}
	}

	return local_matrix_G;
}

vector<vector<double>> Global_G(vector<double> grid_space, double lamda)
{
	int grid_space_size = grid_space.size();
	int num_of_fin_elems = grid_space.size() - 1;

	vector<vector<double>> global_G;

	global_G.resize(grid_space_size);
	for (int i = 0; i < grid_space_size; i++)
	{
		global_G[i].resize(grid_space_size);
	}

	for (int i = 0; i < num_of_fin_elems; i++) 
	{
		vector<vector<double>> local_G = Local_G(lamda, grid_space[i], grid_space[i + 1]);
		global_G[i][i] += local_G[0][0];
		global_G[i][i + 1] += local_G[0][1];
		global_G[i + 1][i] += local_G[1][0];
		global_G[i + 1][i + 1] += local_G[1][1];
	}
	return global_G;
}

vector<double> Local_B(double x0, double x1, double t)
{
	vector<double> local_B(2);
	double h = x1 - x0;

	local_B[0] = ( h / 6 ) * (2 * F_function(x0, t) + F_function(x1, t));
	local_B[1] = ( h / 6 ) * (F_function(x0, t) + 2 * F_function(x1, t));

	return local_B;
}

vector<double> Global_B(vector<double> grid_space, double t)
{
	int num_of_fin_elems = grid_space.size() - 1;
	vector<double> global_B(grid_space.size());
	for (int i = 0; i < num_of_fin_elems; i++)
	{
		vector<double> local_B = Local_B(grid_space[i], grid_space[i + 1], t);
		global_B[i] += local_B[0];
		global_B[i + 1] += local_B[1];
	}
	return global_B;
}

void Conditions(vector<vector<double>>& matrix, vector<double>& vect, vector<double> grid_space, double t) {
	int n = matrix.size();
	for (int i = 0; i < n; ++i) {
		matrix[0][i] = 0;
		matrix[n - 1][i] = 0;
	}
	matrix[0][0] = 1;
	matrix[n - 1][n - 1] = 1;

	vect[0] = U_function(grid_space[0], t);
	vect[n - 1] = U_function(grid_space[n - 1], t);
}

// это для теста. Потому что я мог накосячить в своем лосе, но вроде все ок
void LU(int n, vector<vector<double>>& A, vector<double>& b) {
	for (int i = 1; i < n; ++i) {
		A[i][0] /= A[0][0];
	}

	for (int i = 1; i < n; ++i) {
		for (int j = 1; j < n; ++j) {
			if (i > j) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum += A[i][k] * A[k][j];
				}
				A[i][j] = (A[i][j] - sum) / A[j][j];
			}
			else {
				double sum = 0;
				for (int k = 0; k < i; k++) {
					sum += A[i][k] * A[k][j];
				}
				A[i][j] -= sum;
			}
		}
	}

	for (int i = 1; i < n; ++i) {
		double sum = 0;
		for (int k = 0; k < i; ++k) {
			sum += b[k] * A[i][k];
		}
		b[i] -= sum;
	}

	b[n - 1] /= A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; --i) {
		double sum = 0;
		for (int k = i + 1; k < n; ++k) {
			sum += A[i][k] * b[k];
		}
		b[i] = (b[i] - sum) / A[i][i];
	}
}

void output(vector<double> grid_time, vector<double> grid_space, vector<vector<double>> U)
{

	FILE* fout = fopen("output.txt","w");
	for (int i = 0; i < grid_time.size(); i++)
	{
		
		fprintf(fout, "t = %lf \n", grid_time[i]);
		fprintf(fout, "U\t\t\t\t\t\t\t\t\t\tU0\t\t\t\t\t\t\t\tErr\n");
		for (int j = 0; j < grid_space.size(); j++)
		{
			fprintf(fout, "%2.18le\t\t", U[i][j]);
			fprintf(fout, "%2.18le\t\t", U_function(grid_space[j], grid_time[i]));
			fprintf(fout, "%2.18le\n", U[i][j] - U_function(grid_space[j], grid_time[i]));
		}
	}
	fclose(fout);
}
void main()
{
	vector<double> grid_space;
	vector<double> grid_time;
	int grid_space_size;
	int grid_time_size;
	double chi;
	double sigma;
	double lamda;
	fstream file_input;

	read_grid("grid_space.txt", grid_space_size, grid_space);
	read_grid("grid_time.txt", grid_time_size, grid_time);
	read_params(chi, sigma, lamda);

	vector<vector<double>> U(grid_time_size);

	// вот такая реализация позволяет хранить значения q на всех слоях, но это затратно по памяти
	// я бы мог хранить только 2 вектора для 3-слойной схемы(как пример) и получать на их основе следкющий вектор
	// выводить рез в файл
	// потом переопределять векторы

	for (int i = 0; i < U.size(); i++)
	{
		U[i].resize(grid_space_size);
	}

	// вместо 2 можно указать n-слойность схемы и тогда логично вывести все в функцию.
	// Это норм идея если в будущем я захочу сделать комлекс решателей со множеством схем
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < grid_space_size; j++)
		{
			U[i][j] = U_function(grid_space[j], grid_time[i]);
		}
	}
	for (int i = 2; i < grid_time_size; i++) {

		// это я игрался

		// схема  для равномерного шага по времени равномерная сетка
		// вроде работает
		//double delta_t = 1 / (grid_time[i] - grid_time[i - 1]);

		//vector<vector<double>> G = Global_G(grid_space, lamda);
		//vector<vector<double>> M_chi = Global_M(grid_space, chi);
		//vector<vector<double>> M_sigma = Global_M(grid_space, sigma);
		//vector<vector<double>> global_A = (delta_t * delta_t) * M_chi + (0.5 * delta_t) * M_sigma;

		//vector<double> B(grid_space_size);
		//B = Global_B(grid_space, grid_time[i - 1]);

		//vector<double> D(grid_space_size);
		//D = B + (2 * delta_t * delta_t) * M_chi * U[i - 1] - (delta_t * delta_t) * M_chi * U[i - 2] + (0.5 * delta_t) * M_sigma * U[i - 2] - G * U[i - 1];
		// схема  для равномерного шага по времени КОН

		// неравномерная сетка НАЧ	
		double delta_t = (grid_time[i] - grid_time[i - 2]);
		double delta_t0 = (grid_time[i - 1] - grid_time[i - 2]);
		double delta_t1 = (grid_time[i] - grid_time[i - 1]);

		vector<vector<double>> G = Global_G(grid_space, lamda);
		vector<vector<double>> M_chi = Global_M(grid_space, chi);
		vector<vector<double>> M_sigma = Global_M(grid_space, sigma);

		// в учебнике опечатка см 7.42 7.43 там delta_t1 при М_сигма должна быть
		vector<vector<double>> global_A = (2 / (delta_t * delta_t0)) * M_chi + (delta_t1 / (delta_t * delta_t0)) * M_sigma;

		vector<double> B(grid_space_size);
		B = Global_B(grid_space, grid_time[i - 1]); // опечатка. в учебнике j слой а должен быть j-1

		vector<double> D(grid_space_size);

		D = B - (2 / (delta_t * delta_t1)) * M_chi * U[i - 2] + (2 / (delta_t0 * delta_t1)) * M_chi * U[i - 1] + (delta_t0 / (delta_t * delta_t1)) * M_sigma * U[i - 2] - ((delta_t0 - delta_t1) / (delta_t1 * delta_t0)) * M_sigma * U[i - 1] - G * U[i - 1];
		//// неравномерная сетка КОН

		// Краевые условия
		Conditions(global_A, D, grid_space, grid_time[i]);
		// Решатель Слау
		U[i] = LOS(global_A, D, grid_space_size);
		//LU(grid_space_size, global_A, D);
		//U[i] = D;		
	}
	output(grid_time, grid_space, U);
}