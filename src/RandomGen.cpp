/**********************************************************
RandomGen.cpp
Last modified: 11/28/2023
***********************************************************/

#include "RandomGen.h"
#include <iomanip>      // std::setprecision

/************************ DistrParams ************************/

const char* DistrParams::Equel = " = ";
const char* DistrParams::SepCl = ": ";
const char* DistrParams::SepSCl = "; ";

const char* DistrParams::sMean = "mean";
const char* DistrParams::sSigma = "sigma";
const char* DistrParams::sDistrib = "distribution";

float DistrParams::lnMean;
float DistrParams::lnSigma;
float DistrParams::ssMean = 0;
int	  DistrParams::ssSigma;
int	  DistrParams::rdMean = 0;
int	  DistrParams::rdSigma;

void DistrParams::Init(const pairVal& ln, const pairVal& ss, bool isSS, const pairVal& rd, bool isRD)
{
	lnMean = ln.first,
	lnSigma = ln.second;
	if (IsLn()) {
		if (isSS)
			ssMean = ss.first == -1 ? LnMean() : ss.first,
			ssSigma = int(ss.second);
		if (isRD)
			rdMean = int(rd.first),
			rdSigma = int(rd.second);
	}
}

#define LF	'\n'

void DistrParams::PrintFragDistr(std::ostream& s, const char* title, bool bothDistrs)
{
	const char* sFragment[]{ "fragment","Fragment" };

	if (bothDistrs || !IsSS()) {
		s << title << sFragment[bothDistrs] << " lognorm size " << sDistrib << SepCl;
		if (IsLn()) {
			s << sMean << Equel << lnMean << SepSCl << sSigma << Equel << lnSigma;
			if (bothDistrs)
				s << SepSCl << "Mean" << Equel << std::setprecision(5) << LnMean()
				<< SepSCl << "Mode" << Equel << LnMode();
		}
		else {
			s << "OFF" << SepSCl << sFragment[0] << " length" << Equel << lnMean;
		}
		s << LF;
	}
	if (bothDistrs || IsSS()) {
		s << title << sFragment[bothDistrs] << " size selection " << sDistrib << SepCl;
		if (IsSS())
			s << sMean << Equel << std::setprecision(5) << ssMean
			  << SepSCl << sSigma << Equel << ssSigma << LF;
		else
			s << "OFF\n";
	}
}

void DistrParams::PrintReadDistr(std::ostream& s, const char* title, const char* rTitle)
{
	if (IsRVL())
		s << title << rTitle << " length norm " << sDistrib
		<< SepCl << sMean << Equel << rdMean
		<< SepSCl << sSigma << Equel << rdSigma << LF;
}

/************************ end of DistrParams ************************/

/************************ class Random ************************/

const float PI = 3.14159265f;
int Random::Seed;
float Random::ExpLambda = 1;		// Average rate of occurrence in std::exponential_distribution

int Random::SetSeed(int seed, randval expoBase)
{
	ExpLambda = 10.f / expoBase;
	if (seed)
		Seed = 12345678 + (seed << 24);	// any number with capacity 6-8 (needed RAND_XORSHIFT constructor)
	else {
		time_t tm;
		time(&tm);	// get current time; same as: timer = time(NULL)
		struct tm y2k = { 0 };
		y2k.tm_year = 117; y2k.tm_mday = 1;
		Seed = int(difftime(tm, mktime(&y2k)));	// seconds since January 1, 2017
	}
	return Seed;
}

float Random::ssFactor0;
const float Random::ssFactor1 = 2.5f / (float)sqrt(2 * PI);
/*
Current release (2 vars defined lognormal distribution):
X_lognorm = exp( Normal() * Sigma + Mean );
where Mean = 5.46, Sigma = 0.4,
Normal(): mean = 0, SD = 1
*/

Random::Random() : _normal_x2(DRand())//, distribution(ExpLambda)
{
#ifdef RAND_STD
	srand((unsigned)Seed);
	_seed = Seed;
#elif defined RAND_MT
	mt[0] = Seed;
	for (mti = 1; mti < MERS_N; mti++)
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
	for (int i = 0; i < 37; i++) rand();		// Randomize some more
#elif defined RAND_XORSHIFT
	x = Seed;
	// initialize to fix random generator. Any initialization of y, w, z in fact
	y = x >> 1;	 w = y + 1000;  z = w >> 1;
#endif
}

#ifdef RAND_MT
// Generates 32 random bits
uint32_t Random::rand()
{
	if (mti >= MERS_N) {
		// Generate MERS_N words at one time
		const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
		const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
		static const uint32_t mag01[2] = { 0, MERS_A };
		int kk;
		for (kk = 0; kk < MERS_N - MERS_M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MERS_M] ^ (y >> 1) ^ mag01[y & 1];
		}

		for (; kk < MERS_N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MERS_M - MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
		}

		y = (mt[MERS_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MERS_N - 1] = mt[MERS_M - 1] ^ (y >> 1) ^ mag01[y & 1];
		mti = 0;
	}
	y = mt[mti++];

	// Tempering (May be omitted):
	y ^= y >> MERS_U;
	y ^= (y << MERS_S) & MERS_B;
	y ^= (y << MERS_T) & MERS_C;
	y ^= y >> MERS_L;

	return y;
}
#elif defined RAND_XORSHIFT

uint32_t Random::rand()
{
	uint32_t t = x ^ (x << 11);
	x = y; y = z; z = w;
	return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}
#endif	

int Random::Range(int max)
{
#ifdef RAND_MT
	if (max == 1)	return 1;
	return min(int(DRand() * (max - 1) + 1), max);
#else
	return int(DRand() * --max + 1);
#endif
}

bool Random::Sample(float sample)
{
	return sample >= 1.0 ? true : (sample == 0.0 ? false : (DRand() <= sample));
}

// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
//	Used Agner Fog method
//	return:  value with Gaussian likelihood between about -5 and +5
//	(from -6 to 6 in 1000000000 cycles) 
double Random::Normal()
{
	double normal_x1;	// first random coordinate
	double w;			// radius

	if (_phase) {		// we have a valid result from last call
		_phase = 0;
		return _normal_x2;
	}
	do {				// make two normally distributed variates by Box-Muller transformation
		normal_x1 = 2. * DRand() - 1.;
		_normal_x2 = 2. * DRand() - 1.;
		w = normal_x1 * normal_x1 + _normal_x2 * _normal_x2;
	} while (w >= 1. || w < 1E-30);

	w = sqrt(log(w) * (-2. / w));
	_normal_x2 *= w;	// normal_x1 and normal_x2 are independent normally distributed variates
	_phase = 1;
	return normal_x1 * w;	// return normal distributed value
}

/************************ end of class Random ************************/