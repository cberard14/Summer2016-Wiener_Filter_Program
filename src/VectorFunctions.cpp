#include "El.hpp"
using namespace El;
typedef double Real;

void GetAutocovarianceVector(vector<double>& v, vector<double>&v_out)
{
	// Accepts correlation vector of length n
	// Returns full correlation vector (length 2n-1)
	// length 2n-1 in order to generate Toeplitz matrix
	// leads to symmetric matrix
	int n = v.size();
	for (int i=0;i<2*n-1;i++)
	{
		v_out[i]=v[i+1];
		if (i==(n-1))
		{
			v_out[i]=v[0];
		}
		else if (i>=n) 
		{
			v_out[i]=v[(2*n-1)-i];
		}
	}
}


void Autocovariance1DReal(vector<double>& v,vector<double>& autocovariance)
{
	// Accepts vector of length n
	// Returns autocovariance vector of length n
	int n = v.size();
	int sum,k,kk;
	double mean=0;

	for (int i=0;i<n;i++)
	{
		mean+=(1.0/n)*v[i];
	}

	for (k=0;k<n;k++)
	{
		sum=0;
		for (kk=0;kk<n;kk++)
		{
			sum+=(1.0/(n-1.0))*(v[kk]-mean)*(v[kk-k]-mean);
		}
		autocovariance[k]=sum;
	}
}

void GetRow(DistMatrix<double>& A,int row,vector<double>& v)
{
	// Get a row from a matrix A, fill and update vector v  
	int j,jend;
	jend = A.Width();

	for (j=0;j<jend;j++)
	{
		v[j] = A.Get(row,j);
	}
}


void GetCirculantGeneratingVector(vector<double>& v_original, vector<double>& v_out)
{
	int n = v_original.size();

	for (int i=0;i<2*n-1;i++)
	{
		v_out[i]=v_original[i+1];
		if (i==(n-1)) 
		{
			v_out[i]=v_original[0];
		}
		else if (i>=n) 
		{
			v_out[i]=v_original[i-n+1];
		}
	}
}

void VectorToDistMatrix(vector<double>& v, DistMatrix<double>& M)
{
	// convert a vector to a column matrix
	int n = v.size();
	for (int k=0;k<n;k++)
	{
		M.Set(k,0,v[k]);
	}
}

void DistMatrixToVector(DistMatrix<double>& A, vector<double>& v_out)
{
	int n = A.Height();

	for (int k=0;k<n;k++)
	{
		v_out[k] = A.Get(k,0);
	}
}
