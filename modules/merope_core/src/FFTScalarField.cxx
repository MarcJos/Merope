//! Copyright : see license.txt
//!
//! \brief Scalar FFT fields with FFTW
//!


#include "FFTW3/FFTScalarField.hxx"

#include <cmath>
#include <cstdlib>
#include "../../AlgoPacking/src/Geometry/GeomTypes.hxx"
#include "Field/CovSum.hxx"
#include "Grid/Geostat.hxx"
#include "Grid/Grid.hxx"
#include "VTKinout/VTKStream.hxx"


#include "MeropeNamespace.hxx"


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
	const unsigned flags_i):
	FFTField(grid, flags_i) {
	// Variables number
	nv = 1;
	alloc(IP);
}

void FFTScalarField::RandFunc(int seed) {
	srand(seed);
	checkSpectral("RandFunc");
	size_t nx, ny;
	getDim(nx, ny);
	// Uniform Random values
	vector<double> U(2 * FSize);
	for (auto& iU : U) {
		iU = rand() / (RAND_MAX + 1.);
	}

#pragma omp parallel default(shared)
	{
		// Limits and indexes for loops
		int i2, j2, ncx, ncy, ncx2, ncy2;
		unsigned i, j, ncz, ncz2;
		// Dimension 1,2 and 3 constants
		ncx = nx / 2 + 1;
		if (nx % 2) {
			ncx2 = ncx;
		}
		else {
			ncx2 = ncx - 1;
		}

		ncy = ny / 2 + 1;
		if (ny % 2) {
			ncy2 = ncy;
		}
		else {
			ncy2 = ncy - 1;
		}

		ncz = nz / 2 + 1;
		if (nz % 2) {
			ncz2 = ncz;
		}
		else {
			ncz2 = ncz - 1;
		}
#pragma omp for collapse(2)
		// Loop on x frequencies
		for (i = 0; i < nx; ++i) {
			// Loop on y frequencies
			for (j = 0; j < ny; ++j) {

				cfloat* pF = &F[(i * ny + j) * ncz];
				i2 = i;
				if (i >= (unsigned)ncx) {
					i2 = i - nx;
				}
				bool fx = (i2 != ncx2 && i2);
				j2 = j;
				if (j >= (unsigned)ncy) {
					j2 = j - ny;
				}
				bool fy = (j2 != ncy2 && j2);

				// Loop on z frequencies
				for (unsigned k = 0; k < ncz; ++k) {
					cfloat& rF = pF[k];
					bool fz = (k != ncz2 && k);
					double Dk = 0;
					if (rF[0] > Dk)
						Dk = rF[0];
					double G1, G2;
					size_t I = 2 * ((i * ny + j) * ncz + k);
					VAGauss(G1, G2, U[I], U[I + 1]);
					// Qualify frequencies
					if (fx || fy || fz) {
						// Ordinary frequency
						Dk = sqrt(0.5 * Dk);
						rF[0] = Dk * G1;
						rF[1] = Dk * G2;
					}
					else {
						// Fréquences réelles
						Dk = sqrt(Dk);
						rF[0] = Dk * G1;
						rF[1] = 0;
					}
				}
			}
		}
	}
}

FFTScalarField::FFTScalarField(const Grid& grid, const gaussianField::CovSum& cs,
	const bool IP, int seed, const unsigned flags_i):
	FFTField(grid, flags_i) {
	// Variables number
	nv = 1;
	alloc(IP);
	// Fill the grid with a covariance function
	setCov(grid.getLx(), grid.getLy(), grid.getLz(), cs);
	forward();
	// Generate the random field
	RandFunc(seed);
	backward();
}

FFTScalarField::FFTScalarField(const Grid& grid,
	const std::function<double(array<double, 3>)>& cs, bool IP, int seed, unsigned flags):
	FFTField(grid, flags) {
	// Variables number
	nv = 1;
	alloc(IP);
	// Fill the grid with a covariance function
	setCov(grid.getLx(), grid.getLy(), grid.getLz(), cs);
	forward();
	// Generate the random field
	RandFunc(seed);
	backward();
}

FFTScalarField::FFTScalarField(const Grid& grid,
	const std::function<double(array<double, 2>)>& cs, bool IP, int seed, unsigned flags):
	FFTField(grid, flags) {
	// Variables number
	nv = 1;
	alloc(IP);
	// Fill the grid with a covariance function
	setCov(grid.getLx(), grid.getLy(), grid.getLz(), cs);
	forward();
	// Generate the random field
	RandFunc(seed);
	backward();
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

} // namespace merope
