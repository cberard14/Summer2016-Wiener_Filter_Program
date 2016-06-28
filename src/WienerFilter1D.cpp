#include "../include/VectorFunctions.hpp"
#include "../include/MatrixSetUpFunctions.hpp"
#include "El.hpp"

// This file contains the source code of the 1D Wiener Filter function

void WienerFilter1D(DistMatrix<double>& unlensed, DistMatrix<double>& lense, 
			DistMatrix<double>& lensed, DistMatrix<double>& noise, DistMatrix<double>& recovered)
{
	// this function accepts unlensed, lense, lensed matrices
	// assumes unlensed, lensed, and noise data are column matrices (that is, vectors in DistMatrix form)
	// lense must be an n-by-n matrix
	int n = unlensed.Height(); // determines size of square matrices to be formed

	// initialize matrices needed for Wiener Filter                                                                              
	DistMatrix<double> Co_xx(n,n),Co_nn(n,n);
	DistMatrix<double> xxh(n,n),hxxh(n,n),Co_xx_plus_nn(n,n),filter(n,n);
	Zeros(xxh,n,n); // this matrix holds Co_xx * lensematrix adjoint
	Zeros(hxxh,n,n); // this matrix holds lense * Co_xx * lense adjoint
	Zeros(Co_xx,n,n); // unlensed autocovariance matrix
	Zeros(Co_nn,n,n); // noise autocovariance matrix -- should be diagonal
	Zeros(Co_xx_plus_nn,n,n); 

	// fill noise, unlensed signal autocovariance matrices
	GetAutocovarianceMatrix(noise,Co_nn);
	GetAutocovarianceMatrix(unlensed,Co_xx);

	// Form Co_xx * lens adjoint and lense * Co_xx * lense adjoint 
	Gemm(NORMAL,ADJOINT,Real(1),Co_xx,lense,Real(0),xxh);
	Gemm(NORMAL,NORMAL,Real(1),lense,xxh,Real(0),hxxh);

	// Form the Wiener Filter Co_xx*xxh*[(hxxh+Co_nn)^{-1}]y
	Co_xx_plus_nn = hxxh;
	Co_xx_plus_nn += Co_nn;

	// **** KNOWN ISSUE HERE: Co_xx, Co_nn not numerically HPD so (hxxh+Co_nn) isn't either! But it should be.
	SymmetricInverse(LOWER,Co_xx_plus_nn,false); 
	Gemm(NORMAL,NORMAL,Real(1),xxh,Co_xx_plus_nn,Real(0),filter);
	Gemv(NORMAL,Real(1),filter,lensed,Real(0),recovered);
}
