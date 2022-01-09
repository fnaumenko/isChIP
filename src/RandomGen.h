/**********************************************************
OutTxtFile.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 16.11.2020
-------------------------
Provides random generation functionality
***********************************************************/

#pragma once
#include "common.h"
#include <cmath>
//#include <random>	// std::exponential_distribution

// 'DistrParams' keeps fragment initial lognormal and normal size selection (SS) distributions params
// It is defined here but not in Imitator.h becuase of using in 
static struct DistrParams
{
public:
	typedef pair<int, int> rdParams;

	static float lnMean;	// mean of initial lognormal distribution
	static float lnSigma;	// sigma of initial lognormal distribution
	static float ssMean;	// mean of size selection normal distribution,
							// or 0 if SS is off, or mean of lognormal distribution by default
	static int	 ssSigma;	// sigma of size selection normal distribution
	//static rdParams RDParams;	// mean & sigma of Read normal distribution
	static int	 rdMean;	// mean of size selection normal distribution,
							// or 0 if SS is off, or mean of lognormal distribution by default
	static int	 rdSigma;	// sigma of size selection normal distribution

	// Initializes lognorm fragment and size selection distribution values;
	//	@ln: expectation & standard deviation of frag lognormal distribution
	//	@isFDset: true if paramters od lognorm distr were assigned by user
	//	@ss: expectation & standard deviation of size sel normal distribution
	//	@isSSset: true if size sel is applied
	//	@rd: mean & sigma of Read length distriburion
	//	@isRVLset: true if Read variable length mode is set
	static void Init(const pairVal& ln, bool isFDset, const pairVal& ss, bool isSSset, const pairVal& rd, bool isRVLset);

	// Returns mean of lognormal distribtion
	static float LnMean() { return exp(lnMean + lnSigma * lnSigma / 2); }

	// Returns mode of lognormal distribtion
	static float LnMode() { return exp(lnMean - lnSigma * lnSigma); }

	// Returns true if size selection is ON
	inline static bool IsSS() { return bool(ssMean); }

	// Returns true if read variable length is ON
	inline static bool IsRVL() { return bool(rdMean); }

	// Prints title and the set fragment distribution parameters
	//	@s: outstream to print
	//	@startWord: substring to print first
	//	@all: if true then print both frag dist and size selection, otherwise only one of them
	static void PrintFragDistr(ostream& s, const char* startWord, bool all);
	
	// Prints title and the set read distribution parameters
	//	@s: outstream to print
	//	@startWord: substring to print first
	//	@rTitle: read title
	static void PrintReadDistr(ostream& s, const char* startWord, const char* rTitle);

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
	inline double DRand() { return (double)rand() / RAND_MAX_; }

	static float ssFactor0;			// factor sigma*sqrt(2) in the size sel norm distr
	const static float ssFactor1;	// factor 2.5/sqrt(2PI) in the size sel norm distr

public:
	// Sets and returns real seed
	//	@seed: if 0, random seed
	//	@expoBase: the distance from the site boundary at which the probability increases exponentially
	//		(exonuclease 'headroom' length)
	static	int SetSeed(UINT seed, readlen expoBase = 10);

	Random();

	// Returns random integer within interval [1, max]
	int	Range(int max);

	// Returns true with given likelihood
	//	@sample: probability of returning true; from 0.0. to 1.0
	bool Sample(float sample);

	// Returns random true ot false with probability 0.5
	inline bool Boolean() { return rand() & 0x1; }

	// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
	//	Used Agner Fog method
	//	return:  value with Gaussian likelihood between about -5 and +5
	//	(from -6 to 6 in 1000000000 cycles) 
	double Normal();

	// Random number distrib that produces values according to a lognormal distrib.
	// canonical form:  exp( (Normal()*Sigma + Mean) / LnFactor + LnTerm )
	// About 1.5 times faster then std::lognormal_distribution (<random>)
	inline fraglen Lognormal() {
		return fraglen(exp(Normal() * DistrParams::lnSigma + DistrParams::lnMean));
	}

	// Generates number within interval 0 <= x < max with exponentially decreasing probability (in the max direction)
	fraglen Expo() { return fraglen(-log(1 - rand() / (RAND_MAX_ + 1.0)) / ExpLambda); }
	//inline fraglen Expo() { return fraglen(distribution(generator)); }

};
