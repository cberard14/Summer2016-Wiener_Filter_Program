#include "El.hpp"
#include "./include/VectorFunctions.hpp"
#include "./include/MatrixSetUpFunctions.hpp"
#include "./include/WienerFilter1D.hpp"
using namespace El;
typedef double Real;


int main( int argc, char* argv[] )
{

	Environment env( argc, argv );
	mpi::Comm comm = mpi::COMM_WORLD;
	const int commrank = mpi::Rank(comm);
	const int commsize = mpi::Size(comm);

	try
	{
		// Define flags
		const int n = Input("--n","n",1000);
		const std::string unlensed =
			Input("--unlensed","unlensed",std::string(""));
		const std::string lensed =
			Input("--lensed","lensed",std::string(""));
		const std::string noise  =
			Input("--noise","noise",std::string(""));
		const std::string lense =
			Input("--lense","lense",std::string(""));
		ProcessInput();


		DistMatrix<double> lensecolmatrix(n,1);
		Read(lensecolmatrix, lense);
		vector<double> lensevector(n);
		DistMatrixToVector(lensecolmatrix,lensevector);
		DistMatrix<double> lensematrix(n,n);
		Circulant(lensematrix,lensevector); // periodic convolution with lense vector

		DistMatrix<double> S(n,n),N(n,1),X(n,n), H(n,n);
		DistMatrix<double> xcolumnmatrix(n,1), noisecolumnmatrix(n,1); 
		DistMatrix<double> scolumnmatrix(n,1), ycolumnmatrix(n,1),xhat(n,1);
		vector<double> xvector(n),svector(n),noisevector(n);
		Read( X, unlensed); // read unlensed data
		GetRow(X,0,xvector);
		VectorToDistMatrix(xvector,xcolumnmatrix);

		// form ytest -- check if ytest = y from daniel's simulation      
		DistMatrix<double> signal(n,1),ysimmatrix(n,1);
		Gemv(NORMAL,Real(1),lensematrix,xcolumnmatrix,Real(0),signal);

		Read( S, lensed); // read lensed 
		GetRow(S,0,svector);
		VectorToDistMatrix(svector,scolumnmatrix);

		Read( H, noise); // read noise vector (in DistMatrix form)
		GetRow(H,0,noisevector);
		VectorToDistMatrix(noisevector,noisecolumnmatrix);
	
		ysimmatrix = signal;
		ysimmatrix += noisecolumnmatrix;

		// add signal + noise to get y
		vector<double> yvector(n);
		for (int k=0;k<n;k++)
	  	{
			yvector[k] = svector[k]+noisevector[k];
	  	}
		VectorToDistMatrix(yvector,ycolumnmatrix);

		// OVERWRITE ycolumnmatrix         
		ycolumnmatrix = ysimmatrix;
		Print(ycolumnmatrix,"y = Hx+n simulated");

		WienerFilter1D(xcolumnmatrix,lensematrix,ycolumnmatrix,noisecolumnmatrix,xhat);
		Print(xhat, "result");
	}
	catch( std::exception& e ) { ReportException(e); }
	return 0;
}
