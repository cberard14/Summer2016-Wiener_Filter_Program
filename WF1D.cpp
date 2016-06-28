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
		const std::string raw =
			Input("--raw","raw",std::string(""));
		const std::string noise  =
			Input("--noise","noise",std::string(""));
		const std::string transformation =
			Input("--transformation","transformation",std::string(""));
		const bool print_intermediate = 
			Input("--print_intermediate","print_intermediate",true);
		ProcessInput();

		// Form the convolution matrix
		DistMatrix<double> transcolmatrix(n,1);
		Read(transcolmatrix, transformation);
		vector<double> transvector(n);
		DistMatrixToVector(transcolmatrix,transvector);
		DistMatrix<double> transmatrix(n,n);
		Circulant(transmatrix,transvector); // periodic convolution with lense vector

		DistMatrix<double> N(n,1),X(n,n);
		DistMatrix<double> xcolumnmatrix(n,1), noisecolumnmatrix(n,1); 
		DistMatrix<double> ycolumnmatrix(n,1),xhat(n,1);
		vector<double> xvector(n),noisevector(n);

		Read(X,raw); // read raw data
		GetRow(X,0,xvector);
		VectorToDistMatrix(xvector,xcolumnmatrix);
		if (print_intermediate)
		{
			Print(xcolumnmatrix,"x");
		}

		// form the observable      
		DistMatrix<double> rawsignal(n,1);
		Gemv(NORMAL,Real(1),transmatrix,xcolumnmatrix,Real(0),rawsignal);
		if (print_intermediate)
		{
			Print(rawsignal,"h*x");
		}
		Read(N,noise); // read noise vector (in DistMatrix form)
		GetRow(N,0,noisevector);
		VectorToDistMatrix(noisevector,noisecolumnmatrix);
		// noisecolumnmatrix *= 1.0e4; // dialed up noise for example plot
		if (print_intermediate)
		{
			Print(noisecolumnmatrix,"n");
		}

		ycolumnmatrix = rawsignal;
		ycolumnmatrix += noisecolumnmatrix;

		if (print_intermediate)
		{
			Print(ycolumnmatrix,"y = h*x+n simulated");
		}

		WienerFilter1D(xcolumnmatrix,transmatrix,ycolumnmatrix,noisecolumnmatrix,xhat);
		Print(xhat, "result");
	}
	catch( std::exception& e ) { ReportException(e); }
	return 0;
}
