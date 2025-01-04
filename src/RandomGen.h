/**********************************************************
RandomGen.h
Provides random number, normal and lognormal distribution generation functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/28/2023
***********************************************************/
#pragma once

#include <cmath>
#include <utility>		// std::pair
#include <iostream>

//#include <random>	// std::exponential_distribution

// 'DistrParams' keeps fragment initial lognormal and normal 'size selection' (SS) distributions params
static struct DistrParams
{
private:
	static const char* Equel;
	static const char* SepCl;		// colon separator
	static const char* SepSCl;		// semicolon separator

	static const char* sMean;
	static const char* sSigma;
	static const char* sDistrib;

public:
	typedef std::pair<float, float> pairVal;
	typedef std::pair<int, int> rdParams;

	static float lnMean;	// mean of lognormal distribution
	static float lnSigma;	// sigma of lognormal distribution
	static float ssMean;	// mean of 'size selection' normal distribution,
							// or 0 if SS is off, or mean of lognormal distribution by default
	static int	 ssSigma;	// sigma of 'size selection' normal distribution
	static int	 rdMean;	// mean of Read length normal distribution, or 0 if RD is off
	static int	 rdSigma;	// sigma of Read length normal distribution

	// Initializes distribution values;
	//	@param ln: expectation & standard deviation of frag lognormal distribution
	//	@param ss: expectation & standard deviation of size sel normal distribution
	//	@param isSS: true if size sel is applied
	//	@param rd: mean & sigma of Read length distriburion
	//	@param isRD: true if Read variable length mode is set
	static void Init(const pairVal& ln, const pairVal& ss, bool isSS, const pairVal& rd, bool isRD);

	// Returns mean of lognormal distribution
	static float LnMean() { return IsLn() ? exp(lnMean + lnSigma * lnSigma / 2) : 1.5f * lnMean; }

	// Returns mode of lognormal distribution
	static float LnMode() { return exp(lnMean - lnSigma * lnSigma); }

	// Returns true if lognormal distribution is ON
	static bool IsLn() { return lnSigma; }

	// Returns true if size selection is ON
	static bool IsSS() { return bool(ssMean); }

	// Returns true if read variable length is ON
	static bool IsRVL() { return bool(rdMean); }

	// Prints title and the fragment distribution parameters
	//	@param s: outstream to print
	//	@param title: pre-title to print first
	//	@param bothDistrs: if true then print both frag dist and size selection, otherwise only one of them
	static void PrintFragDistr(std::ostream& s, const char* title, bool bothDistrs);
	
	// Prints title and the read distribution parameters
	//	@param s: outstream to print
	//	@param title: pre-title to print first
	//	@param rTitle: read title
	static void PrintReadDistr(std::ostream& s, const char* title, const char* rTitle);

} distrParams;


//#define RAND_STD			// rand() (Windows) or rand_r(int *seed) (Linux)
// Mersenne Twister by Agner Fog, 2008-11-16 http://www.agner.org/random/
//#define RAND_MT
// Xorshift by George Marsaglia http://en.wikipedia.org/wiki/Xorshift
#define RAND_XORSHIFT

#ifdef RAND_STD
#define RAND_MAX_	RAND_MAX
#ifdef __unix__
#define rand()	rand_r(&_seed)
#endif
#else
#define RAND_MAX_	0xFFFFFFFF	//(65536.*65536.)	// 2^32
#endif

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
  // Compilers supporting C99 or C++0x have stdint.h defining these integer types
#include <stdint.h>
#define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
  // 16 bit systems use long int for 32 bit integer.
typedef   signed long int int32_t;
typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Older Microsoft compilers have their own definition
typedef   signed __int32  int32_t;
typedef unsigned __int32 uint32_t;
typedef   signed __int64  int64_t;
typedef unsigned __int64 uint64_t;
#define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
  // This works with most compilers
typedef signed int          int32_t;
typedef unsigned int       uint32_t;
typedef long long           int64_t;
typedef unsigned long long uint64_t;
#define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif

// 'Random' encapsulates random number generator.
class Random
{
#ifdef RAND_MT
	// Version of Mersenne Twister: MT11213A or MT19937
#if 0	// Constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else   // Constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
#endif

private:
	static int Seed;
	static float ExpLambda;		// Average rate of occurrence in std::exponential_distribution
	//std::default_random_engine generator;
	//std::exponential_distribution<float> distribution;

#ifdef RAND_STD
	int _seed;
#elif defined RAND_MT
	int mti;					// Index into mt
	uint32_t y, mt[MERS_N];		// State vector

	// Generates 32 random bits
	uint32_t rand();
#elif defined RAND_XORSHIFT
	uint32_t x, y, z, w;

	// Generates 32 random bits
	uint32_t rand();
#endif

	double	_normal_x2;		// second random coordinate (for normal())
	short	_phase;			// phase (for normal RNG)

protected:
	// Generates random double number within interval 0 <= x < 1
	double DRand() { return (double)rand() / RAND_MAX_; }

	static float ssFactor0;			// factor sigma*sqrt(2) in the size sel norm distr
	const static float ssFactor1;	// factor 2.5/sqrt(2PI) in the size sel norm distr

public:
	using randval = uint32_t;

	// Sets and returns real seed
	//	@param seed: if 0, random seed
	//	@param expoBase: the distance from the site boundary at which the probability increases exponentially
	//		(exonuclease 'headroom' length)
	static	int SetSeed(int seed, randval expoBase = 10);

	Random();

	// Returns random integer within interval [1, max]
	int	Range(int max);

	// Returns true with given likelihood
	//	@param sample: probability of returning true; from 0.0. to 1.0
	bool Sample(float sample);

	// Returns random true ot false with probability 0.5
	bool Boolean() { return rand() & 0x1; }

	// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
	//	Used Agner Fog method
	//	@returns: value with Gaussian likelihood between about -5 and +5 (from -6 to 6 in 1000000000 cycles) 
	double Normal();

	// Random number distrib that produces values according to a lognormal distrib.
	// canonical form: exp( (Normal()*Sigma + Mean) / LnFactor + LnTerm )
	// About 1.5 times faster then std::lognormal_distribution (<random>)
	randval Lognormal()
	{
		return randval(exp(Normal() * DistrParams::lnSigma + DistrParams::lnMean));
	}

	// Generates number within interval 0 <= x < max with exponentially decreasing probability (in the max direction)
	randval Expo() { return randval(-log(1 - rand() / (RAND_MAX_ + 1.0)) / ExpLambda); }
	//randval Expo() { return fraglen(distribution(generator)); }
};
