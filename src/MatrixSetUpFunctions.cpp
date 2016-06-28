#include "El.hpp"
#include "../include/VectorFunctions.hpp"
using namespace El;
typedef double Real;


void SetSubMatrix(DistMatrix<double>& A,DistMatrix<double>& Asub,int istart,int jstart)
{
	// set a portion of matrix A, beginning at indices istart,jstart, to be 
	// equal to the smaller matrix Asub
	int i,j,iloc,jloc,iend,jend;
	iend=Asub.Height();
	jend=Asub.Width();
	iloc=0;

	for (i=istart;i<istart+iend;i++)
	{
		jloc=0;
		for (j=jstart;j<jstart+jend;j++)
		{
			A.Set(i,j,Asub.Get(iloc,jloc));
			jloc+=1;
		}
		iloc+=1; 
	}  
}

void GetAutocovarianceMatrix(DistMatrix<double>& v_data, DistMatrix<double>&  CovarMatrix)
{
	// take read matrix (column matrix) to serve as vector, 
	// output covariance matrix CovarMatrix    
																																						
	int n = v_data.Height();
	vector<double> v_data_vector(n),v_covariance(n),v_covar_full(2*n-1);

	DistMatrixToVector(v_data, v_data_vector);
	Autocovariance1DReal(v_data_vector,v_covariance);
	GetAutocovarianceVector(v_covariance,v_covar_full);
	Toeplitz(CovarMatrix,n,n,v_covar_full);
}

void FillBlockToeplitz(DistMatrix<double>& A,DistMatrix<double>& X)
{
	// WHAT IF X WAS SIZE 2N-1 x 2N-1? -- X is a matrix of correlation vectors...                                                                                                                                                       
	// function overwrites A, which becomes block Toeplitz                                                                                                                
	// each generating vector for each submatrix C, Cdiag, is a row of X                                                                                                  

	int nsubmatrices = X.Height(); // 2*n-1 distinct blocks
	int blocksize = X.Width();
	int nskip = X.Width();   // submatrices same size as row vectors of X
	int sizeA = A.Height();
	int i,j,k,iloc,jloc;
	DistMatrix<double> C,Cdiag;
	vector<double> vrow(blocksize);
	vector<double> vrow2(blocksize*2-1);
	
	Zeros(C,blocksize,blocksize);
	Zeros(Cdiag,blocksize,blocksize);
	
	iloc=0;
	jloc=0;

	// handle diagonal matrices separately         
	GetRow(X,0,vrow);
	GetCirculantGeneratingVector(vrow,vrow2);
	Toeplitz(Cdiag,blocksize,blocksize,vrow2);
		
	for (k=0;k<sizeA;k+=nskip)
	{
		SetSubMatrix(A,Cdiag,k,k);
	}
	
	// fill lower off-diagonal blocks              
	int rownum=0;
	for (i=nskip;i<pow(nskip,2);i+=nskip)
	{
		rownum+=1;
		for (j=0;j<(pow(nskip,2)-i);j+=nskip)
		{
			GetRow(X,rownum,vrow);
			GetCirculantGeneratingVector(vrow,vrow2);
			Toeplitz(C,blocksize,blocksize,vrow2);
			SetSubMatrix(A,C,i+j,j);
		}
	} 

	// fill upper off-diagonal blocks                                         
	rownum=nskip*2-1; // nskip-1;
	for (i=nskip;i<pow(nskip,2);i+=nskip)
	{
		rownum-=1;
		for (j=0;j<(pow(nskip,2)-i);j+=nskip)
		{
			GetRow(X,rownum,vrow);
			GetCirculantGeneratingVector(vrow,vrow2);
			Toeplitz(C,blocksize,blocksize,vrow2);
			SetSubMatrix(A,C,j,j+i);
		}
	}
}
