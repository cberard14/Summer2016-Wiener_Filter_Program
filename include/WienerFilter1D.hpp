#ifndef WIENERFILTER_H
#define WIENERFILTER_H

#include "El.hpp"
using namespace El;

void WienerFilter1D(DistMatrix<double>& unlensed, DistMatrix<double>& lense,
                    DistMatrix<double>& lensed, DistMatrix<double>& noise, DistMatrix<double>& recovered);

#endif
