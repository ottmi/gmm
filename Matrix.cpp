#include "Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(unsigned int dim)
{
	_dim = dim;
	vector<double> v = vector<double>(_dim, .0);
	_m = vector<vector<double> >(_dim, v);
}

Matrix::~Matrix()
{
	// TODO Auto-generated destructor stub
}

void Matrix::zero()
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = .0;
}

void Matrix::set(vector<vector<double> > x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = x[i][j];
}

void Matrix::set(vector<double> x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = x[i * _dim + j];
}

void Matrix::setDiag(double x)
{
	for (unsigned int i = 0; i < _dim; i++)
		_m[i][i] = x;
}

void Matrix::setOffDiag(double x)
{
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			if (i != j) _m[i][j] = x;
}

void Matrix::setRow(unsigned int row, vector<double> x)
{
	for (unsigned int j = 0; j < _dim; j++)
		_m[row][j] = x[j];
}

void Matrix::setCol(unsigned int col, vector<double> x)
{
	for (unsigned int i = 0; i < _dim; i++)
		_m[i][col] = x[i];
}

void Matrix::setEntry(unsigned int row, unsigned int col, double x)
{
	_m[row][col] = x;
}

void Matrix::update(Matrix& x)
{
	double sum = 0;
	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			sum += x.getEntry(i, j);

	for (unsigned int i = 0; i < _dim; i++)
		for (unsigned int j = 0; j < _dim; j++)
			_m[i][j] = x.getEntry(i, j) / sum;
}

double Matrix::getEntry(unsigned int row, unsigned int col)
{
	return _m[row][col];
}

double Matrix::getRowSum(unsigned int row)
{
	double sum = .0;
	for (unsigned int i = 0; i < _dim; i++)
		sum += _m[row][i];
	return sum;
}

double Matrix::getColSum(unsigned int col)
{
	double sum = .0;
	for (unsigned int i = 0; i < _dim; i++)
		sum += _m[i][col];
	return sum;
}

double Matrix::determinant()
{
	double det = 0;
	if (_dim == 2)
	{
		det = _m[1][1] * _m[2][2] - _m[1][2] * _m[2][1];
	} else if (_dim == 3)
	{
		det = _m[0][0] * _m[1][1] * _m[2][2] + _m[0][1] * _m[1][2] * _m[2][0] + _m[0][2] * _m[1][0] * _m[2][1] - _m[0][0] * _m[1][2] * _m[2][1]
				- _m[0][1] * _m[1][0] * _m[2][2] - _m[0][2] * _m[1][1] * _m[2][0];
	} else if (_dim == 4)
	{
		det = _m[0][3] * _m[1][2] * _m[2][1] * _m[3][0] - _m[0][2] * _m[1][3] * _m[2][1] * _m[3][0] - _m[0][3] * _m[1][1] * _m[2][2] * _m[3][0]
				+ _m[0][1] * _m[1][3] * _m[2][2] * _m[3][0] + _m[0][2] * _m[1][1] * _m[2][3] * _m[3][0] - _m[0][1] * _m[1][2] * _m[2][3] * _m[3][0]
				- _m[0][3] * _m[1][2] * _m[2][0] * _m[3][1] + _m[0][2] * _m[1][3] * _m[2][0] * _m[3][1] + _m[0][3] * _m[1][0] * _m[2][2] * _m[3][1]
				- _m[0][0] * _m[1][3] * _m[2][2] * _m[3][1] - _m[0][2] * _m[1][0] * _m[2][3] * _m[3][1] + _m[0][0] * _m[1][2] * _m[2][3] * _m[3][1]
				+ _m[0][3] * _m[1][1] * _m[2][0] * _m[3][2] - _m[0][1] * _m[1][3] * _m[2][0] * _m[3][2] - _m[0][3] * _m[1][0] * _m[2][1] * _m[3][2]
				+ _m[0][0] * _m[1][3] * _m[2][1] * _m[3][2] + _m[0][1] * _m[1][0] * _m[2][3] * _m[3][2] - _m[0][0] * _m[1][1] * _m[2][3] * _m[3][2]
				- _m[0][2] * _m[1][1] * _m[2][0] * _m[3][3] + _m[0][1] * _m[1][2] * _m[2][0] * _m[3][3] + _m[0][2] * _m[1][0] * _m[2][1] * _m[3][3]
				- _m[0][0] * _m[1][2] * _m[2][1] * _m[3][3] - _m[0][1] * _m[1][0] * _m[2][2] * _m[3][3] + _m[0][0] * _m[1][1] * _m[2][2] * _m[3][3];
	} else
	{
		for (int c = 0; c < (int) _dim; c++)
		{
			Matrix m = minor(0, c);
			double d = m.determinant();
			det += (1 - c % 2 - c % 2) * _m[0][c] * d;
		}
	}

	return det;
}

Matrix Matrix::minor(const unsigned int row, const unsigned int col)
{
	Matrix ret(_dim - 1);
	vector<vector<double> > m(_dim - 1, vector<double>(_dim - 1));

	if (row < _dim && col < _dim)
	{
		unsigned int r2 = 0;
		for (unsigned int r1 = 0; r1 < _dim; r1++)
		{
			if (r1 != row)
			{
				unsigned int c2 = 0;
				for (unsigned int c1 = 0; c1 < _dim; c1++)
				{
					if (c1 != col)
					{
						m[r2][c2] = _m[r1][c1];
						c2++;
					}
				}
				r2++;
			}
		}
	} else
		throw("Matrix::minor(): index out of range");

	ret.set(m);
	return ret;
}

void Matrix::print()
{
	cout.precision(8);
	for (unsigned int i = 0; i < _dim; i++)
	{
		for (unsigned int j = 0; j < _dim; j++)
			cout << fixed << _m[i][j] << " ";
		cout << endl;
	}
}
