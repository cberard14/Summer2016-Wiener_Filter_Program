#ifndef MATRIXSETUPFUNCTIONS_H
#define MATRIXSETUPFUNCTIONS_H

#include "El.hpp"
using namespace El;

void SetSubMatrix(DistMatrix<double>& A,DistMatrix<double>& Asub,int istart,int jstart);
void GetAutocovarianceMatrix(DistMatrix<double>& v_data, DistMatrix<double>&  CovarMatrix);
void FillBlockToeplitz(DistMatrix<double>& A,DistMatrix<double>& X);

#endif
