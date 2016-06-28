#ifndef VECTORFUNCTIONS_H
#define VECTORFUNCTIONS_H

#include "El.hpp"
using namespace El;

void GetAutocovarianceVector(vector<double>& v, vector<double>&v_out);
void Autocovariance1DReal(vector<double>& v,vector<double>& autocovariance);
void GetRow(DistMatrix<double>& A,int row,vector<double>& v);
void GetCirculantGeneratingVector(vector<double>& v_original, vector<double>& v_out);
void VectorToDistMatrix(vector<double>& v, DistMatrix<double>& M);
void DistMatrixToVector(DistMatrix<double>& A, vector<double>& v_out);

#endif
