//! Copyright : see license.txt
//!
//! \brief
//!


#include "../../Geometry/include/GeomConstants.hxx"
#include "../../merope_core/include/geneOrientations.hxx"


namespace merope {
namespace geneOrientation {

inline void test() {
	// create source of randomness, and initialize it with non-deterministic seed
	random_device r;
	seed_seq seed{ r() };
	mt19937 eng{ seed };
	// a distribution that takes randomness and produces values in specified range
	uniform_real_distribution<> dist(-6, 6);

	// Same seed -> each run produce the same numbers
	default_random_engine generator(2);
	uniform_real_distribution<> distribution(-6, 6);

	// Bruit aleatoire entre [-b,b] = 2b(rand()/intervalTotalDeRetourDeLaFonctionRand - 0.5)
	double b = 6;
	double aa = 2 * b / ((double)RAND_MAX + 1.);
	double Br;
	srand(4);

	cout << "\n# 1:production aleatoire de l'aleatoire 2:C++11 default (generator) 3:srand\n";
	for (size_t i = 0; i < 10; ++i) {
		Br = aa * rand() - b;
		cout << dist(eng) << " " << distribution(generator) << " (" << generator
			<< ") " << Br << endl;
	}
	cout << endl;
}

inline void EulerToAxis(const double* const ang, double* const vec) {
	// Premiere colonne de la matrice de passage a
	vec[0] = cos(ang[2]) * cos(ang[0])
		- sin(ang[2]) * cos(ang[1]) * sin(ang[0]);
	vec[1] = cos(ang[2]) * sin(ang[0])
		+ sin(ang[2]) * cos(ang[1]) * cos(ang[0]);
	vec[2] = sin(ang[2]) * sin(ang[1]);
	// Deuxieme colonne de la matrice de passage a
	vec[3] = -cos(ang[0]) * sin(ang[2])
		- cos(ang[1]) * cos(ang[2]) * sin(ang[0]);
	vec[4] = cos(ang[0]) * cos(ang[1]) * cos(ang[2])
		- sin(ang[0]) * sin(ang[2]);
	vec[5] = cos(ang[2]) * sin(ang[1]);
	// Troisieme colonne de la matrice de passage a
	vec[6] = sin(ang[1]) * sin(ang[0]);
	vec[7] = -1 * sin(ang[1]) * cos(ang[0]);
	vec[8] = cos(ang[1]);
}

inline void checkOrth(const double* const vec) {
	double p12 = vec[0] * vec[3] + vec[1] * vec[4] + vec[2] * vec[5];
	double p13 = vec[0] * vec[6] + vec[1] * vec[7] + vec[2] * vec[8];
	double p23 = vec[3] * vec[6] + vec[4] * vec[7] + vec[5] * vec[8];
	if (p12 > 1e-14 or p13 > 1e-14 or p23 > 1e-14)
		throw runtime_error(
			"Generation d'orientations aleatoires : le produit scalaire des vecteurs de la base n'est pas nul\n");

	double n1 = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	double n2 = sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]);
	double n3 = sqrt(vec[6] * vec[6] + vec[7] * vec[7] + vec[8] * vec[8]);
	if (abs(n1 - 1.) > 1e-14 or abs(n2 - 1.) > 1e-14 or abs(n3 - 1.) > 1e-14)
		throw runtime_error(
			"Generation d'orientations aleatoires : la norme des vecteurs de la base n'est pas egale a 1\n");
}

inline void PremCad(const double* const vec, double* const v) {
	// Normalisation du vecteur
	double n = sqrt(vec[2] * vec[2] + vec[5] * vec[5] + vec[8] * vec[8]);
	if (vec[2] > 0)
		v[0] = vec[2] / n;
	else
		v[0] = -vec[2] / n;

	if (vec[5] > 0)
		v[1] = vec[5] / n;
	else
		v[1] = -vec[5] / n;

	if (vec[8] > 0)
		v[2] = vec[8] / n;
	else
		v[2] = -vec[8] / n;

	// Inversion des axes pour realiser les inegalites 0 <= v[1] <= v[0] <= v[2]
	double tmp;
	if (v[2] < v[0]) {
		tmp = v[0];
		v[0] = v[2];
		v[2] = tmp;
	}
	if (v[2] < v[1]) {
		tmp = v[1];
		v[1] = v[2];
		v[2] = tmp;
	}
	if (v[0] < v[1]) {
		tmp = v[0];
		v[0] = v[1];
		v[1] = tmp;
	}
}

inline void printPremCad(const double* const ang, const double* const vec,
	ofstream& f) {
	double v[3];
	PremCad(vec, v);
	f << ang[0] << " " << ang[1] << " " << ang[2] << " " << v[0] << " " << v[1]
		<< " " << v[2] << endl;
}

inline void printAllCad(const double* const ang, const double* const vec,
	ofstream& f) {
	f << ang[0] << " " << ang[1] << " " << ang[2] << " " << vec[2] << " "
		<< vec[5] << " " << vec[8] << endl;
}

inline void printLicos(const unsigned short n, const double* const vec, ofstream& f) {
	f << "Material '" << n << "'\n";
	f << "AnisotropicAxises {{'" << vec[0] << "','" << vec[1] << "','" << vec[2]
		<< "'},{'" << vec[3] << "','" << vec[4] << "','" << vec[5]
		<< "'}}\n";
	f << "EndOfMaterial\n\n";
}

inline void printTmfft(const unsigned short n, const double* const ang, ofstream& f) {
	f << "Material " << n << "\n";
	f << "AnisotropicAxises<Angles> {" << ang[0] * 180 / m_PI << ","
		<< ang[1] * 180 / m_PI << "," << ang[2] * 180 / m_PI << "}\n";
	f << "EndOfMaterial\n\n";
}

inline void printNoFormatting(const unsigned short n, const double* const vec,
	ofstream& f) {
	f << n << " " << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3]
		<< " " << vec[4] << " " << vec[5] << endl;
}

inline void randomOrient0(const unsigned short N, const unsigned short seed_i,
	ofstream& f) {
	double ang[3], vec[9];
	default_random_engine gen(seed_i);
	uniform_real_distribution<> dist(0, 2 * m_PI), dist2(0, 1);
	for (unsigned short i = 0; i < N; ++i) {
		// Angles D'Euler : ang[0] = psi, ang[1] = theta, ang[2] = phi
		// en notation de Bunge : ang[0] = phi_1, ang[1] = phi, ang[2] = phi_2
		ang[0] = dist(gen);
		ang[1] = acos(1 - 2 * dist2(gen));
		ang[2] = dist(gen);

		EulerToAxis(ang, vec);
		checkOrth(vec);
		printNoFormatting(i + 1, vec, f);
		// printLicos(i+1,vec,f);
		// printTmfft(i,ang,g);
		// printPremCad(ang,vec,f);
		// printAllCad(ang,vec,f);
	}
}

inline void readOrient0(const char* const outEBSD, ofstream& f, ofstream& g) {
	ifstream febsd(outEBSD);
	if (!febsd)
		throw logic_error("Le fichier sortie d'EBSD n'existe pas");

	string lin;
	unsigned short pixOK, n = 0;
	double ang[3], vec[9], tmp;
	while (febsd.good()) {
		getline(febsd, lin);
		istringstream ist(lin);
		ist >> pixOK;
		if (ist.good() and pixOK > 0) {
			ist >> tmp >> tmp >> ang[0] >> ang[1] >> ang[2];
			ang[0] = ang[0] * m_PI / 180;
			ang[1] = ang[1] * m_PI / 180;
			ang[2] = ang[2] * m_PI / 180;
			EulerToAxis(ang, vec);
			checkOrth(vec);
			++n;
			printLicos(n, vec, f);
			printTmfft(n - 1, ang, g);
			// printPremCad(ang,vec,f);
		}
	}
	febsd.close();
}

void printRandomOrient_3D(const unsigned short N,
	const unsigned short seed_i, std::string nameFile) {
	ofstream file(nameFile);
	randomOrient0(N, seed_i, file);
}

}  // namespace geneOrientation
}  // namespace merope

