#include "Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(unsigned int dim)
{
	_dim = dim;
	vector<double> v = vector<double>(_dim, .0);
	_m = vector< vector<double> >(_dim, v);
}


Matrix::~Matrix()
{
	// TODO Auto-generated destructor stub
}


void Matrix::zero()
{
	for (unsigned int i=0; i<_dim; i++)
		for (unsigned int j=0; j<_dim; j++)
			_m[i][j] = .0;
}


void Matrix::set(vector< vector<double> > x)
{
	for (unsigned int i=0; i<_dim; i++)
		for (unsigned int j=0; j<_dim; j++)
			_m[i][j] = x[i][j];
}


void Matrix::set(vector<double> x)
{
	for (unsigned int i=0; i<_dim; i++)
		for (unsigned int j=0; j<_dim; j++)
			_m[i][j] = x[i*_dim + j];
}


void Matrix::setDiag(double x)
{
	for (unsigned int i=0; i<_dim; i++)
		_m[i][i] = x;
}


void Matrix::setOffDiag(double x)
{
	for (unsigned int i=0; i<_dim; i++)
		for (unsigned int j=0; j<_dim; j++)
			if (i != j)
				_m[i][j] = x;
}


void Matrix::setRow(unsigned int row, vector<double> x)
{
	for (unsigned int j=0; j<_dim; j++)
		_m[row][j] = x[j];
}


void Matrix::setCol(unsigned int col, vector<double> x)
{
	for (unsigned int i=0; i<_dim; i++)
		_m[i][col] = x[i];
}


void Matrix::print()
{
	cout.precision(8);
	for (unsigned int i=0; i<_dim; i++)
	{
		for (unsigned int j=0; j<_dim; j++)
			cout << fixed << _m[i][j] << " ";
		cout << endl;
	}
}
