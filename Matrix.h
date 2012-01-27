#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
using namespace std;

class Matrix
{
public:
	Matrix(unsigned int dim);
	virtual ~Matrix();
	void zero();
	void set(vector< vector<double> > m);
	void set(vector<double> m);
	void setDiag(double x);
	void setOffDiag(double x);
	void setRow(unsigned int row, vector<double> m);
	void setCol(unsigned int col, vector<double> m);

	double getEntry(unsigned int row, unsigned int col);
	double getRowSum(unsigned int row);
	double getColSum(unsigned int col);

	void print();


private:
	unsigned int _dim;
	vector< vector<double> > _m;
};

#endif /* MATRIX_H_ */
