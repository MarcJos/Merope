//! Copyright : see license.txt
//!
//! \briefScalar FFT fields with FFTW
//!


#include "FFTW3/FFTScalarField.hxx"

#include <cmath>
#include <cstdlib>
#include "../../Geometry/include/GeomTypes.hxx"
#include "Field/CovSum.hxx"
#include "Grid/Geostat.hxx"
#include "Grid/Grid.hxx"
#include "VTKinout/VTKStream.hxx"


namespace merope {

// 2 Variables aléatoires normales indépendantes
// X,Y : Variables aléatoires normales indépendantes
// U1,U2: Variables aléatoires uniformes indépendantes
inline void VAGauss(double& X, double& Y, double& U1, double& U2) {
	double R = sqrt(-2 * log(U1));
	X = cos(2 * m_PI * U2) * R;
	Y = sin(2 * m_PI * U2) * R;
}

FFTScalarField::FFTScalarField(const Grid& grid, const bool IP,
	const unsigned flags_i) :
	FFTField(grid, flags_i) {
	// Variables number
	nv = 1;
	alloc(IP);
}

void FFTScalarField::prepareCovarianceInFourier() {
	auto function = [&](bool, bool, bool,
		auto& rF) {
			double Dk = max(0., rF[0]);
			rF[0] = Dk;
			rF[1] = 0;
		};
	this->loopOnFrequencies(function);
}

void FFTScalarField::randFunc(int seed) {
	srand(seed);
	// Uniform Random values
	vector<double> U(2 * FSize);
	for (auto& iU : U) {
		iU = rand() / (RAND_MAX + 1.);
	}
	size_t i_courant = 0;
	auto function = [&](bool fx, bool fy, bool fz,
		auto& rF) {
			double Dk = max(0., rF[0]);
			double G1, G2;
			size_t I = 2 * i_courant;
			VAGauss(G1, G2, U[I], U[I + 1]);
			// Qualify frequencies
			if (fx || fy || fz) {
				// Ordinary frequency
				Dk = sqrt(0.5 * Dk);
				rF[0] = Dk * G1;
				rF[1] = Dk * G2;
			} else {
				// Fréquences réelles
				Dk = sqrt(Dk);
				rF[0] = Dk * G1;
				rF[1] = 0;
			}
			i_courant++;
		};
	this->loopOnFrequencies(function, false);
}

FFTScalarField::FFTScalarField(const Grid& grid, const gaussianField::CovSum& cs,
	const bool IP, int seed, const unsigned flags_i) :
	FFTField(grid, flags_i) {
	build(&grid, cs, IP, seed, flags_i, false);
}

FFTScalarField::FFTScalarField(const Grid& grid,
	const std::function<double(array<double, 3>)>& cs, bool IP, int seed,
	bool showCovariance, unsigned flags_i) :
	FFTField(grid, flags_i) {
	build(&grid, cs, IP, seed, flags_i, showCovariance);
}

FFTScalarField::FFTScalarField(const Grid& grid,
	const std::function<double(array<double, 2>)>& cs, bool IP, int seed,
	bool showCovariance, unsigned flags_i) :
	FFTField(grid, flags_i) {
	build(&grid, cs, IP, seed, flags_i, showCovariance);
}


void FFTScalarField::toVTKCELL(VTKstream& fvtk, const char* const cname) const {
	checkSpatial("toVTKCELL");

	size_t nx, ny;
	getDim(nx, ny);
	fvtk.setCELL(fSize);

	fvtk << "SCALARS " << cname << " float\n";
	fvtk << "LOOKUP_TABLE default" << endl;

	// Header
	float x;
	const rfloat* f2;
	unsigned i, j, k;
	switch (d) {
	case 1:
		for (i = 0, f2 = f; i < nz; ++i, ++f2) {
			x = (float)(*f2);
			fvtk.write(x);
		}
		break;
	case 2:
		for (j = 0; j < nz; ++j) {
			for (i = 0; i < ny; ++i) {
				x = (float)(f[j + nzb * i]);
				fvtk.write(x);
			}
		}
		break;
	case 3:
		for (k = 0; k < nz; k++) {
			for (j = 0; j < ny; j++) {
				for (i = 0; i < nx; i++) {
					x = (float)(f[k + nzb * (j + ny * i)]);
					fvtk.write(x);
				}
			}
		}
	}
	fvtk << endl;
}

}  // namespace merope
