#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <clocale>
#include <fstream>
#include <iostream>
#include<vector>

	using namespace std;
	int Size = 3;
	double epsilon = 1e-13;

	vector<vector<double>> operator+(const vector<vector<double>>& matrix, const vector<vector<double>>& matrix2) {
		vector<vector<double>> c;
		int size = matrix.size();
		c.resize(size);
		for (auto& vect : c) {
			vect.resize(size);
		}
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				c[i][j] = matrix[i][j] + matrix2[i][j];
			}
		}
		return c;
	}

	vector<vector<double>> operator*(const double& a, const vector<vector<double>>& matrix) {
		vector<vector<double>> c;
		int size = matrix.size();
		c.resize(size);
		for (auto& vect : c) {
			vect.resize(size);
		}
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				c[i][j] = matrix[i][j] * a;
			}
		}
		return c;
	}

	vector<double> operator*(const vector<vector<double>>& matrix, const vector<double>& vect) {
		vector<double> c;
		int size = matrix.size();
		c.resize(size);

		for (int i = 0; i < size; ++i) {
			double sum = 0;
			for (int j = 0; j < size; ++j) {
				sum += matrix[i][j] * vect[j];
			}
			c[i] = sum;
		}
		return c;
	}

	vector<double> operator* (const double& c, const vector<double>& vect)
	{
		vector<double> result = vect;
		for (int i = 0; i < result.size(); i++)
		{
			result[i] = vect[i] * c;
		}
		return result;
	}

	vector<double> operator+(const vector<double>& vect1, const vector<double>& vect2) {
		vector<double> c;
		int size = vect1.size();
		c.resize(size);

		for (int i = 0; i < size; ++i) {
			c[i] = vect1[i] + vect2[i];
		}
		return c;
	}

	vector<double> operator-(const vector<double>& vect1, const vector<double>& vect2) {
		vector<double> c;
		int size = vect1.size();
		c.resize(size);

		for (int i = 0; i < size; ++i) {
			c[i] = vect1[i] - vect2[i];
		}
		return c;
	}

	void Read_matrix(string matrix_filename, vector<vector<double>>& matrix_A, int Size)
	{
		fstream file_start_cons;
		file_start_cons.open(matrix_filename);
		if (file_start_cons.is_open())
		{
			for (int i = 0; i < Size; i++)
			{
				for (int j = 0; j < Size; j++)
				{
					file_start_cons >> matrix_A[i][j];
				}
			}
		}
		file_start_cons.close();
	}
	void Read_vector(string vector_filename, vector<double>& vector, int Size)
	{
		fstream file_start_cons;
		file_start_cons.open(vector_filename);
		if (file_start_cons.is_open())
		{
			for (int i = 0; i < Size; i++)
			{
				file_start_cons >> vector[i];
			}
		}
		file_start_cons.close();
	}

	double ScalVectOnVect(vector<double> vector1, vector<double> vector2, int Size)
	{
		int N = vector1.size();
		double sum = 0;
		for (int i = 0; i < N; i++)
		{
			sum += (vector1[i] * vector2[i]);
		}
		return sum;
	}


	vector<double> LOS(vector<vector<double>> matrix_A, vector<double> vector_B, int Size) // Для любой матрицы A
	{
		vector<double> r0(Size);
		vector<double> r_next(Size);
		vector<double> z0(Size);
		vector<double> z_next(Size);
		vector<double> x_next(Size);
		vector<double> p0(Size);
		vector<double> p_next(Size);
		vector<double> x0(Size);
		vector<double> result = x0;

		double normF;
		double sum = 0;
		double alpha = 0;
		int iter_count = 0;

		int max_iter = 10000;

		double normRkkvad = 0;
		r0 = vector_B - matrix_A * result; // (1 2 3) - (14 32 50) = (-13 - 30 - 47)
		z0 = r0;
		p0 = matrix_A * z0;

		do
		{
			iter_count++;
			//alpha nach
			double alphaChisl = 0;
			double alphaZnam = 0;
			alphaChisl = ScalVectOnVect(p0, r0, Size);
			alphaZnam = ScalVectOnVect(p0, p0, Size);

			alpha = alphaChisl / alphaZnam;
			//alpha kon

			x_next = result + alpha * z0;

			//rk
			r_next = r0 - alpha * p0;

			normRkkvad = ScalVectOnVect(r_next, r_next, Size);


			//betta nach
			double bettaChisl = 0;
			double bettaZnam = 0;
			double betta = 0;
			vector<double> temp;
			temp = matrix_A * r_next;

			bettaZnam = ScalVectOnVect(p0, p0, Size);
			bettaChisl = ScalVectOnVect(p0, temp, Size);

			betta = -(bettaChisl / bettaZnam);
			//betta kon

			//zk
			z_next = r_next + betta * z0;


			p_next = matrix_A * r_next + betta * p0;


			result = x_next;
			r0 = r_next;
			z0 = z_next;
			p0 = p_next;

		} while (normRkkvad > epsilon && iter_count < max_iter);

		return result;
	}
