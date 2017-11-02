#include "isChIP.h"
#include "Imitator.h"
#include <algorithm>    // std::sort

static const char* SignDbg = "## ";	// Marker of output debug info
static const char* cAverage = "AVERAGE";

/************************  class Average ************************/
void Average::operator+=(UINT val)
{
	if( _keep )
		if( _summator < (ULLONG_MAX-val) ) {
			_summator += val;
			_count++;
			_keep = _count < ULONG_MAX;		//  cut off treatment due to counter overflow
			if( !_keep && Imitator::Verbose(vDEBUG) )
				Err("block count", cAverage).Throw();
		}
		else {
			_keep = false;					//  cut off treatment due to summator overflow
			if( Imitator::Verbose(vDEBUG) )	
				Err("block summation", cAverage).Throw();
		}
}
/************************  end of class Average ************************/

/************************ class Random ************************/
int Random::Seed = 123456;		// any number with capacity 6-8 (needed RAND_XORSHIFT constructor)

// Sets  and returns seed
//	@random: if true, random seed
int Random::SetSeed(bool random)
{
	if( random ) {
		time_t tm;
		time(&tm);	// get current time; same as: timer = time(NULL)
		struct tm y2k = {0};
		y2k.tm_year = 117; y2k.tm_mday = 1;
		Seed = int(difftime(tm, mktime(&y2k)));	// seconds since January 1, 2017
	}
	return Seed;
}

Random::Random()
{
#ifdef RAND_STD
	srand( (unsigned)Seed );
	_seed = Seed;
#elif defined(RAND_MT)
	Init0(Seed);
	for (int i = 0; i < 37; i++) BRandom();		// Randomize some more
#elif defined(RAND_XORSHIFT)
	x = Seed;
	// initialize to fix random generator. Any initialization of y, w, z in fact
	y = x >> 1;	 w = y + 1000;  z = w >> 1;
#endif

	_phase = 0;
}

#ifdef RAND_MT
void Random::Init0(int seed)
{
	// Seed generator
	const uint32_t factor = 1812433253UL;
	mt[0]= seed;
	for (mti=1; mti < MERS_N; mti++) {
		mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
	}
}

// Generates 32 random bits
uint32_t Random::BRandom()
{
	if (mti >= MERS_N) {
		// Generate MERS_N words at one time
		const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
		const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
		static const uint32_t mag01[2] = {0, MERS_A};
		int kk;
		for (kk=0; kk < MERS_N-MERS_M; kk++) {    
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

		for (; kk < MERS_N-1; kk++) {    
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

		y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
		mti = 0;
	}
	y = mt[mti++];

	// Tempering (May be omitted):
	y ^=  y >> MERS_U;
	y ^= (y << MERS_S) & MERS_B;
	y ^= (y << MERS_T) & MERS_C;
	y ^=  y >> MERS_L;

	return y;
}
#endif

#ifdef RAND_XORSHIFT
// Generates 32 random bits
uint32_t Random::IRand()
{
    uint32_t t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}
#endif	

// Generates random double number in the interval 0 <= x < 1
inline double	Random::DRand()
{
#ifdef RAND_STD
	return (double)rand() / RAND_MAX;
#elif defined(RAND_MT)
	return (float)BRandom() * (1./(65536.*65536.));		// Multiply by 2^(-32)
#elif defined(RAND_XORSHIFT)
	return (double)IRand() / 0xFFFFFFFF;
#endif
}

// Returns random integer within interval [1, max]
#ifdef RAND_MT
int Random::Range(int max) {
	if( max == 1 )	return 1;
	int res = int(DRand() * (max-1) + 1);
	return min(res, max);
}
#else
inline int Random::Range(int max) {
	return int(DRand() * (max - 1) + 1);
}
#endif


inline bool Random::Boolean()
{ 
#ifdef RAND_STD
	return rand() & 0x1;
#elif defined RAND_MT
	return DRand() >= 0.5F;
#elif defined RAND_XORSHIFT
	return IRand() & 0x1;
#endif
}

// Returns true with given likelihood
//	@sample: probability of returning true; from 0.0. to 1.0
bool Random::RequestSample(float sample)
{
	return sample == 1.0 ? true : (sample == 0.0 ? false : DRand() <= sample);
}

// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
//	return: value with gaussian likelihood between about -5 and +5
double Random::Normal() {
	double normal_x1;		// first random coordinate (normal_x2 is member of class)
	double w;				// radius
	
	if (_phase) {			// we have a valid result from last call
		_phase = 0;
		return normal_x2;
	}    
	do {					// make two normally distributed variates by Box-Muller transformation
		normal_x1 = 2. * DRand() - 1.;
		normal_x2 = 2. * DRand() - 1.;
		w = normal_x1 * normal_x1 + normal_x2 * normal_x2;
	} while (w >= 1. || w < 1E-30);

	w = sqrt( log(w) * (-2./w) );
	_phase = 1;
	normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
	return normal_x1 * w;
}

// Normal distribution method discussed in Knuth and due originally to Marsaglia
//double Random::Normal0()
// The fastest on http://c-faq.com/lib/gaussian.html. Modified to non-static to use in different threads.
// It's approximately equal to std::normal_distribution (<random>) by using standart rand(),
// but 1.3 times faster by using fastRand().
//{
//	if(_phase == 0) {
//		do {
//			_V1 = 2. * DRand() - 1;
//			_V2 = 2. * DRand() - 1;
//			_S = _V1 * _V1 + _V2 * _V2;
//		} while(_S >= 1 || _S == 0);
//		_phase = 1;
//		return _V1 * sqrt(-2 * log(_S) / _S);
//	} 
//	_phase = 0;
//	return _V2 * sqrt(-2 * log(_S) / _S);
//}

/************************ end of class Random ************************/

/************************  class LognormDistribution ************************/

float LognormDistribution::_RelSigma;
float LognormDistribution::_RelMean;
float LognormDistribution::_szSelSigma;
Average LognormDistribution::_Average;

// Random number distribution that produces floating-point values according to a lognormal distribution,
// with or without output accumulation to calculate average
// About 1.5 times faster then std::lognormal_distribution (<random>)
fraglen LognormDistribution::NextWithAccum() {
	if( !_saveAverage )	return Next();
	short res = Next();
	_Average += res;
	return res;
}

/*
float LognormDistribution::LognormalFactorCalc()
{
	if( _lognFactor == 0 )
	{
		int i=0, k=0;
		float y;
		for(i=0; i<96; i++)
		{
			y = exp( GaussStrict()/Options::LnFactor + Options::LnTerm );
			if( y >= Options::FragLenMin && y <= Options::FragLenMax )
				k++;
			_lognMean += y;
		}
		_lognFactor = (float)k/i;
		_lognMean = _lognMean/i;
	}
	return _lognFactor;
}

float LognormDistribution::LognormalMeanCalc()
{
	float s = 1.8;
	float m = 4.0;
	float x = exp( (m + s*s/2)/1 + 1 );
	for(int x=0; x<7; x++)
	{
		cout<< "x: " << x << endl;
		for(int i=1; i<7; i++)
			cout << i*100 << '\t' << exp( (m + s*s)/i + (float)x ) << endl;
	}
	cout  << endl;
	return 1.0;
	if( _lognMean == 0 )
	{
		double sqrt2PI = sqrt(2*PI);
		float y, z;
		float half;
		float m=0, s=1;
		int i=0, limit = 40;
		for(int x=-limit; x<=limit; x++)
		{
			z = (float)x/10;
			//y =	exp( (-log(z)*log(z)/2) ) / (z*sqrt(2*PI));
			//y =	exp(-(z*z/2)/sqrt2PI);
			y =	exp(-((log(z)-m)*(log(z)-m)/(2*s*s)) ) / (z*sqrt(2*PI*s));
			cout << y << endl;
			_lognMean += y;
			i++;
		}
		_lognMean = _lognMean/i;
	}
	return _lognMean;
}
*/
/************************  end of class LognormDistribution ************************/

/************************  class Amplification ************************/
short Amplification::Coefficient = 1;
Average Amplification::_Average;

inline Amplification::Amplification() //: _memFactor(3), calcAverage(false) {
{
	_memFactor = 3;
	calcAverage = false;
	_fractions = Coefficient > 1 ? new Fraction[_memFactor*Coefficient] : NULL;
}

void Amplification::Reset(int fragLen)
{
	_fragLen = fragLen;
	_initCnt = -(Coefficient/2);
	if( _fractions ) {
		// fill the array by amplificated fractions
		_currCnt = -_initCnt;
		memset(_fractions, 0, FractionsSize());	// clear shifts
		Split(true);
		_currCnt = -(Coefficient/2) - 1;		// clear current index to prepare for reading
	}
	else {
		_currCnt = -2;
		_count = 0;
	}
}

 // Recursively fills array _fractions by splitted fragments (segments).
 //  firstCall: true if it is a first call
void Amplification::Split(bool firstCall)
{
	short leftFrac, rightFrac;//, lenMin = Imitator::FragLenMin;
	_count = _currCnt;
	for(short i=(_initCnt+_currCnt)/2; i<_count; i++)	// start from second half of added fractions
	{	// temporary use rightFrac as a Current Fragment Length
		if( firstCall )
			rightFrac = _fragLen;
		else if((rightFrac=_fractions[i].length) < 0)	// skip "hole"
			continue;
		leftFrac = short(random->Range(rightFrac));
		rightFrac -= leftFrac;							// use rightFrac as Right Fraction Length
		if(leftFrac >= Imitator::FragLenMin) {			// save left fraction on the same position
			_fractions[i].length = leftFrac;			// shift doesn't change for left fraction
			if(rightFrac >= Imitator::FragLenMin) {		// save right fraction on the new position
				if(_currCnt+1 >= _memFactor*Coefficient) {
					// increase memory for array if necessary
					Fraction* tmp = _fractions;
					_fractions = new Fraction[(++_memFactor) * Coefficient];
					memset(_fractions, 0, FractionsSize());		// clear new block
					memcpy(_fractions, tmp, sizeof(struct Fraction)*(_memFactor-1)*Coefficient);
					delete [] tmp;
					if( Imitator::Verbose(vDEBUG) ) {
						Mutex::Lock(Mutex::OUTPUT);
						cout << "Amplification's buffer was increased; new memFactor = "
							 << (int)_memFactor << endl;
						Mutex::Unlock(Mutex::OUTPUT);
					}
				}
				// save right fraction on the new position
				_fractions[_currCnt++].Set(rightFrac, _fractions[i].shift + leftFrac);	
			}
		}
		else			// discard left fraction
			if( rightFrac >= Imitator::FragLenMin)
				// save right fraction on the same position
				_fractions[i].Set(rightFrac, _fractions[i].shift + leftFrac);	
			else
				// discard right fraction: position should be free (a "hole" is created)
				// mark "hole" and clear shift for case saving right fraction in the future
				_fractions[i].Set(-1, 0);	
	}
	if( !firstCall )
		_initCnt = _count;
	if( _currCnt > _initCnt )
		Split(false);
}

 // Gets length and shift of amplificated fraction (MDA)
 //  shift: returned shift of fraction
 //  return: length of fraction or 0 if fractions are ended
short Amplification::GetFraction	(short* shift)
{
	if( ++_currCnt < 0 ) {				// "virtual" part of array: return initial fragment
		_count++;
		*shift = 0;
		return _fragLen;
	}
	if( _currCnt == _initCnt ) {				// end of records
		if( calcAverage )	_Average += _count;	// correct Amplification Mean
		return 0;
	}
	while( _fractions[_currCnt].length < 0 )	// skip "holes"
		_currCnt++;
	_count++;
	*shift = _fractions[_currCnt].shift;
	return _fractions[_currCnt].length;			// return current fraction
}

/************************  end of class Amplification ************************/

/************************  class ChromsThreads ************************/
struct ChrSize {
	chrid	ID;
	chrlen	Size;

	ChrSize(chrid cID, chrlen size) : ID(cID), Size(size) {}
	// for sorting by descent
	inline bool operator < (const ChrSize& chrSize) const {	return (chrSize.Size < Size); }
};

// Distributes chroms among threads possibly according equally runtime,
// mining runtime is in proportion to vhroms treated length
ChromsThreads::ChromsThreads(threadnumb thrCnt, const ChromFiles& chrFiles)
{
	chrid	cCnt = chrFiles.TreatedCount(),
			cCntInTread = cCnt / thrCnt + 1;	// +1 for case unzero remainder of the division
	threadnumb i;
	vector<ChromsThread>::iterator itThr;

	// initialize threads
	_threads.reserve(UINT(thrCnt));
	for(i=0; i<thrCnt; i++)
		_threads.push_back( ChromsThread(cCntInTread) );
	
	if( thrCnt == 1 ) {			// one thread
		itThr=_threads.begin();
		itThr->Numb = 1;
		itThr->chrIDs.reserve(cCnt);
		// insert treated chromosomes
		for(ChromFiles::cIter it=chrFiles.cBegin(); it!=chrFiles.cEnd(); it++)
			if( chrFiles.IsTreated(it) )
				itThr->chrIDs.push_back(CID(it));
	}
	else {
		chrlen chrLen;
		ULLONG sumSize = 0;
		vector<ChrSize> sizes;		// temporary vector of treated chroms to sort by descent
		
		sizes.reserve(cCnt);
		// fill temporary vector of treated chroms
		for(ChromFiles::cIter it=chrFiles.cBegin(); it!=chrFiles.cEnd(); it++)
			if( Imitator::All || chrFiles.IsTreated(it) ) {
				chrLen = chrFiles.ChromTreatLength(it, 1);
				sizes.push_back( ChrSize(CID(it), chrLen ) );
				sumSize += chrLen;
			}
		// in genome each next chromosome has decreasing size, but this isn't always true.
		// So do sort it for any case
		sort(sizes.begin(), sizes.end());	// by descent

		i=0;
		char step = 1;
		for(vector<ChrSize>::iterator it=sizes.begin(); it!=sizes.end(); it++) {
			_threads[i].sumSize += it->Size;
			_threads[i].chrIDs.push_back(it->ID);
			// go through threads, then in reverse order, and again
			if( (i+step)/thrCnt )	step = -1;
			else if( i+step < 0 )	step = 1;
			else					i += step;
		}
		// sort by sumSize ascending to set the minimal first. First is the main.
		sort(_threads.begin(), _threads.end());
		// set thread's numbers: 1 first
		i=1;
		for(itThr=_threads.begin(); itThr!=_threads.end(); itThr++)
			itThr->Numb = i++;
	}
}

void ChromsThreads::Print()
{
	vector<chrid>::iterator it1;
	vector<ChromsThread>::iterator it = _threads.begin();
	//cout << "THREADS:  min weight " << it->sumSize << EOL;	// size of the first (minimal) thread

	for(; it!=_threads.end(); it++) {
		cout << "thr " << int(it->Numb)
			 //<< ":\tweight " << setprecision(3) << double(it->sumSize)/minSize << SepClTab;
			 << ":\tweight " << (it->sumSize) << SepClTab;
		for(it1=it->chrIDs.begin(); it1!=it->chrIDs.end(); it1++) {
			cout << Chrom::AbbrName(*it1) << BLANK;
			if( Chrom::NameLength(*it1) == 1 )	cout << BLANK;	// padding
		}
		cout << endl;
	}
}
/************************  end of class ChromsThreads ************************/

/************************ class AvrFragLengths ************************/
#define SEP	"-"

AvrFragLengths::AvrFragLengths	(const string path)
{
	_isChanged = false;
	ostringstream oss;
	oss << path
		<< Options::GetIVal(oFRAG_LEN)	<< SEP
		<< Options::GetIVal(oFRAG_DEV)	<< SEP
		<< Options::GetIVal(oMEAN)		<< SEP
		<< Options::GetIVal(oSIGMA)		<< SEP
		<< Options::GetIVal(oLN_FACTOR)	<< SEP
		<< Options::GetDVal(oLN_TERM)	<< ".txt";

	TabFile file(_fileName = oss.str(), TxtFile::ALL, 2);
	const char* currLine;
	while( currLine = file.GetLine() )
		_avrLens.push_back(AvrFragLen(
			short(file.IntField(0)),
			file.FloatField(1) ));
}

#define AVRG_CAPACITY	6

AvrFragLengths::~AvrFragLengths()
{
	if( _isChanged ) {
		LineFile file(_fileName, TAB);

		file.BeginWrite(2+2+CHRLEN_CAPAC);
		for(vector<AvrFragLen>::iterator it = _avrLens.begin(); it != _avrLens.end(); it++)
			file.WriteLine(it->AmplCoeff, it->Value, AVRG_CAPACITY);
		file.Write();
	}
}

bool AvrFragLengths::Get(short amplCoeff, float *sizeFactor, float *commonAverage, float *savedAverage) const
{
	vector<AvrFragLen>::const_iterator it;
	for(it = _avrLens.begin(); it != _avrLens.end(); it++)
		if( it->AmplCoeff == amplCoeff ) {
			*savedAverage = it->Value;
			*sizeFactor = _avrLens[0].Value;
			*commonAverage = _avrLens[1].Value;
			return true;
		}
	return false;
}

void AvrFragLengths::Add(short amplCoeff, float sizeFactor, float commonAverage, float savedAverage)
{
	if( _avrLens.size() == 0 ) {
		_avrLens.push_back(AvrFragLen(-1, sizeFactor));
		_avrLens.push_back(AvrFragLen(0, commonAverage));
		_avrLens.push_back(AvrFragLen(amplCoeff, savedAverage));
	}
	else {
		// insert savedAverage in ascending order
		_isChanged = false;
		vector<AvrFragLen>::iterator it;
		for(it = _avrLens.begin(); it != _avrLens.end(); it++)
			if( it->AmplCoeff > amplCoeff ) {
				_avrLens.insert(it, AvrFragLen(amplCoeff, savedAverage));
				_isChanged = true;
				break;
			}
		if( !_isChanged )	// insert maximal to the end
			_avrLens.push_back(AvrFragLen(amplCoeff, savedAverage));
		// average commonAverage
		_avrLens[1].Value = (_avrLens[1].Value + commonAverage) / 2;
	}
	_isChanged = true;
}
#ifdef DEBUG
void AvrFragLengths::Print() const
{
	//for(auto it = _avrLens.begin(); it != _avrLens.end(); it++)
	vector<AvrFragLen>::const_iterator it=_avrLens.begin();
	for(; it<_avrLens.end(); it++)
		cout << it->AmplCoeff << TAB << it->Value << endl;
}
#endif	// DEBUG
/************************ end of class AvrFragLengths ************************/

/************************ class ChromCutter ************************/

// Creates instance
//	@imitator: the owner
//	@csThread: thread contained treated chromosomes
//	@calcAverage: true if averages should be calculated
Imitator::ChromCutter::ChromCutter(
	const Imitator* imitator, ChromsThreads::ChromsThread* csThread, bool calcAverage) :
	_chrFiles(imitator->_chrFiles),
	_partoFile(&(imitator->_oFile)),	// keep pointers to main output files
	_isTerminated(false),
	_thread(*csThread)
{
	ClearCounters();
	_ampl.calcAverage = calcAverage;
	_ampl.random = &_lnDist;
	if( _thread.IsTrial() )
		_lnDist._saveAverage = true;	// to calculate samples
	else if( _thread.IsSlave() )
		_partoFile = new OutFile(imitator->_oFile, _thread.Numb);
	_partoFile->SetEmptyMode(csThread->IsTrial());
}

Imitator::ChromCutter::~ChromCutter ()
{
	if( _thread.IsSlave() ) {
		if( !_isTerminated ) {
			_partoFile->Write();
			InterlockedExchangeAdd(&(Imitator::TotalSlaveWrReadsCnt), _partoFile->Count());
		}
		delete _partoFile;
	}
}

// Outputs count and percent of writes Reads
//  gr: fore/background
//  title: title printed at first and separated by ": ", or NULL
void Imitator::ChromCutter::OutputReadCnt(Imitator::eGround gr, const char* title)
{
	static BYTE charsCnt[2] = {0, 0};
	// print count of reads
	if(title)	cout << title << SepCl;
	cout << setw(Imitator::DigitsCnt[gr]) << setfill(BLANK)
		 << (_wrReadsCnt[gr]<<OutFile::PairedEnd());
	// print percent
	if( NoAmplification ) {
		float percent = Percent(_wrReadsCnt[gr], _selReadsCnt[gr]);
		int charCnt = SFCNT(percent, 3);
		charCnt = max(charCnt, 3);
		if( charsCnt[gr] < charCnt )	charsCnt[gr] = charCnt;
		cout << sPercent(percent, 3, charsCnt[gr]);
	}
}

// Outputs chromosome's name and treatment info
//  @nts: current chromosome
//	@timer: current timer to thread-saves time output or NULL
//	@exceedLimit: true if limit is exceeded
void Imitator::ChromCutter::OutputChromInfo (const Nts& nts, Timer& timer, bool exceedLimit)
{
	if( !Verbose(vRT) )	return;
	Mutex::Lock(Mutex::OUTPUT);
	if( !RegularMode )
		OutputReadCnt(Imitator::FG, TestMode ? "fg" : NULL);
	if( TestMode )
		OutputReadCnt(Imitator::BG, "  bg");
	if(Verbose(vPAR)) {
		cout << "\tN" << SepCl << sPercent(nts.CountN(), nts.Length(), 2, 4, false);
		if( !LetN )
			cout << ", discard " << sPercent(nts.DefLength(), nts.Length(), 3, 0, false);
		//timer.Stop("\t", false, false);
		if( exceedLimit )
			cout << " ! stopped due to exceeding of reads limit of " << Read::MaxCount;
	}
	timer.Stop("\t", false, false);
	cout << endl;
	Mutex::Unlock(Mutex::OUTPUT);
}

// Sets terminate's sign and output message
void Imitator::ChromCutter::Terminate(const char*msg)
{
	_isTerminated = true;
	cerr << "thread " << int(_thread.Numb) << SepCl << msg << endl;
}

// Treats chromosomes given for current thread
//	@singleThread: true if single thread execution: just for print chrom name
void Imitator::ChromCutter::Execute(bool singleThread)
{
	chrid	cID;
	BedF::cIter	cit;	// template chrom's iterator
	ULONG	k, n, cnt;	// count of cells
	chrlen	currPos, cntFtrs;
	short	res;
	Timer	timer;

	try {
		for(vector<chrid>::iterator it=_thread.chrIDs.begin(); it!=_thread.chrIDs.end(); it++) {
			timer.Start();
			cID = *it;
			OutputChromName(cID, singleThread);			// print before cutting
			cnt = CellsCnt << _chrFiles[cID].Numeric();	// multiply twice for digits
			cntFtrs = ( Bed && (cit=Bed->GetIter(cID)) != Bed->cEnd() ) ?
				Bed->FeaturesCount(cit) : 0;
			Nts nts(_chrFiles.FileName(cID), LetN);
			const Featr defRegion = nts.DefRegion();
			ClearCounters();
			_chrName = Chrom::AbbrName(cID) + string(Read::NmDelimiter);
			_partoFile->BeginWriteChrom(cID);

			for(n = 0; n < cnt; n++) {
				res = 0;
				// random shift from the beginning
				currPos=nts.Start() + _lnDist.Range(Imitator::FragLenMax);	
				for(k=0; k < cntFtrs; k++)
					if( res = CutChrom(nts, &currPos, Bed->Feature(cit, k), true) )
						break;	
				if( res < 0 )	// achievement of limit
					break;	
				// add background after last 'end' position
				if( Imitator::TreatOutFtrs
				&& (res = CutChrom(nts, &currPos, defRegion, ControlMode)) < 0 )
					break;				// achievement of limit
			}
			OutputChromName(cID, !singleThread);			// print before cutting
			OutputChromInfo(nts, timer, res < 0);
			for(BYTE i=0; i<2; i++) {
				InterlockedExchangeAdd(&(Imitator::TotalSelReadsCnts[i]), _selReadsCnt[i]);
				InterlockedExchangeAdd(&(Imitator::TotalWrReadsCnts[i]), _wrReadsCnt[i]);
			}
			//timer.Stop("\t", false);
			if( res < 0 )	break;		// achievement of limit
		}
	}
	catch(const Err &e)			{ Terminate(e.what()); }
	catch(const exception &e)	{ Terminate(e.what()); }
	catch(...)					{ Terminate("Unregistered error in thread"); }
	if( Verbose(vDEBUG) )	{
		Mutex::Lock(Mutex::OUTPUT);
		//cout << SignDbg << "end thread " << int(_thread.Numb) << endl;
		cout << "thr " << int(_thread.Numb) << ":\tend" << endl;
		Mutex::Unlock(Mutex::OUTPUT);
	}
}

// Cuts chromosome 
//	@nts: cutted chromosome
//	@currPos: cutting start position
//	@feature: current treated feature
//	@fgInFeature: if true accept foreground keeps insinde feature
//	return: 0 if success,
//		1 if end chromosome is reached (continue treatment),
//		-1 if limit is achieved (cancel treatment)
int Imitator::ChromCutter::CutChrom	(
	const Nts& nts,
	chrlen* const currPos,
	const Featr& feature,
	bool fgInFeature)
{
	bool	selByCorrBounds;// selection by corrected bounds; always true for BG
	bool	reverse;		// reverse Read (set minus strand)
	BYTE	indGr;			// ground index: 0 - FG, 1 - BG
	fraglen fragLen,		// fragment's length
			fracLen,		// fraction's length
			fracShift,		// fraction's start position within fragment
			szselDev;		// fragment's length deviation
	chrlen	start,			// corrected feature's start position
			end;			// corrected feature's end position
	short	addRdRes;
	// variables used to imitate unexpected fragments on the opposite side of BS
	chrlen	fracCentre,		// centre of fraction
			featrBound;		// feature's bound: start or end

	readscr	score[] = { 				// current score:
		Imitator::UniformScore ? 1 : feature.Score,	// FG
		1											// BG
	};			

	for(; *currPos <= feature.End; *currPos += fragLen)	// ChIP: control right mark
	{
		fragLen = _lnDist.NextWithAccum();
		szselDev = _lnDist.NormalNext();
		if( szselDev < 0 )	szselDev = -szselDev;
		//szselDev = 0;
		if( fragLen < FragLenMin - szselDev )	continue;	// size selection: skip short fragment

		// control left mark: 
		// TestMode: foreground (indGr==0) is inside and
		// background (indGr==1) outside features for chromosomes from BED-file;
		// for other chromosomes background (indGr==0) is inside feature,
		// which is the whole chromosome.
		// ControlMode: foreground (indGr==0) is always inside feature,
		// which is the whole chromosome.
		indGr = BYTE(fgInFeature ^ (*currPos + fragLen >= feature.Start));
		if( RequestSample(indGr) ) {
			if(indGr)							// background?
				selByCorrBounds = true;
			else {								// foreground
				start = feature.Start;
				end = feature.End;
				if( Imitator::FlatLen ) {
					fraglen halfShrinkLen = Range(fragLen - FlatLen)>>1;
					//start = *currPos + Range(fragLen - FlatLen);
					start += halfShrinkLen;
					end -= halfShrinkLen;
				}
				//start = *currPos + (FlatLen ?	// smoothing ON?
				//	Range(fragLen - FlatLen):	// correct by smoothing
				//	BSLen);						// correct by BS
				selByCorrBounds = *currPos <= end && *currPos + fragLen >= start;
			}
			if( selByCorrBounds && _lnDist.RequestSample(score[indGr]) )
				for( _ampl.Reset(fragLen); fracLen = _ampl.GetFraction(&fracShift); )
					if( fracLen <= FragLenMax + szselDev	// size selection: skip long fragments
					&& RequestAdjSample() ) {				// adjusted sample?
						*currPos += fracShift;				// in case of BG fracShift is always 0
						
						if( !indGr && Imitator::StrandAdmix ) {	// FG and admix opposite strand?
	// Admix opposite strand:
	// if fragment is on the left site of BS,
	// likelihood of negative strand is linearly decreasing from 1 to 0.5 while moving fragment right
	// if fragment is on the right site of BS,
	// likelihood of negative strand is linearly decreasing from 0.5 to 0 while moving fragment right
	// if fragment's centre is inside BS,
	// likelihood of negative strand is 0.5
							fracCentre = *currPos + (fracLen>>1);
							if(fracCentre < start)			featrBound = start;
							else if(fracCentre > end)		featrBound = end;
							else							goto A;
							reverse = _lnDist.RequestSample( float(featrBound - *currPos)/fracLen );
						}
						else	// likelihood of neg strand is always 0.5
A:							reverse = _lnDist.Boolean();
						
						addRdRes = AddRead(nts, *currPos, fracLen, reverse);

						if( addRdRes < 0 )		return 1;	// end of chromosome: continue treatment
						if( addRdRes > 0 ) {
							// increment of writed Reads in thread
							// file may be NULL in case of SetSample()
							_wrReadsCnt[indGr]++;
							if( Read::IncrementCounter() )							
								return -1;	// achieved of limit: cancel treatment
						}
					}
		}
		if( fragLen <= FragLenMax + szselDev )
			_selReadsCnt[indGr]++;	// increment of selected Reads in thread
	}
	return 0;
}

// Adds read(s) to output file
//	@nts: cutted chromosome
//	@pos: current cutting position
//	@fragLen: length of current fragment
//	@isReverse: true if read has negative strand
//	return: -1 if fragment is NULL,
//		0 if limitN is exceeded,
//		1 if Read(s) is(are) added,
//		2 if output file is NULL
inline int Imitator::ChromCutter::AddRead(
	const Nts& nts, chrlen pos, short fragLen, bool isReverse)
{
	return _partoFile->AddRead(_chrName, nts, 
		// +1 since counters are not incremented yet
		Read::IsNameAsNumber() ? _wrReadsCnt[0] + _wrReadsCnt[1] + 1 : 0,
		pos, fragLen, isReverse);
}

/************************ end of class ChromCutter ************************/

/************************  class Imitator ************************/
ULONG	Imitator::TotalWrReadsCnts[GR_CNT] = {0,0};	// total counts of writed Reads: [0] - fg, [1] - bg
ULLONG	Imitator::TotalSelReadsCnts[GR_CNT] = {0,0};// total counts of selected Reads: [0] - fg, [1] - bg
ULLONG	Imitator::TotalSlaveWrReadsCnt = 0;			// total count of all writed Reads in slaved threads
ULONG	Imitator::CellsCnt;
float	Imitator::AdjSample = 1.0;
float	Imitator::Samples[GR_CNT] = {1.0,1.0};
readlen	Imitator::FlatLen = 0;
fraglen Imitator::FragLenMin;	// Minimal length of selected fragments
fraglen Imitator::FragLenMax;	// Maximal length of selected fragments
BYTE	Imitator::Verb;
BYTE	Imitator::DigitsCnt[GR_CNT] = {0,0};
bool	Imitator::TreatOutFtrs;
bool	Imitator::LetN;
bool	Imitator::UniformScore;
bool	Imitator::StrandAdmix;
bool	Imitator::All;
eMode	Imitator::Mode;			// Current task mode
Imitator	*Imitator::Imit = NULL;
const BedF	*Imitator::Bed = NULL;

// Prints chromosome's name and treatment info
//	@cID: chromosomes ID
//	@isOutput: true if chromosomes name should be printed
void Imitator::OutputChromName(chrid cID, bool isOutput)
{
	if( isOutput && Verbose(vRT) ) {
		cout << Chrom::TitleName(cID) << SepClTab;
		fflush(stdout);
	}
}

// Prints number of recorded Reads
//	@gr: fore/back ground
//	@title: title printed before and separated by ": "
void Imitator::OutputReadCnt(eGround gr, const char* title)
{
	ULONG cnt = TotalWrReadsCnts[gr]<<OutFile::PairedEnd();
	cout << title << SepCl << cnt;		// print count
	if( NoAmplification )		// print percent
		cout << sPercent(cnt, TotalSelReadsCnts[gr]<<OutFile::PairedEnd(), 3);
}

// Runs task in current mode and write result to output files
void Imitator::Execute()
{
	bool res = RegularMode ? CutRegular() : CutGenome();
	_oFile.Write();

	if( Verbose(vRES) ) {
		cout << "Total recorded reads" << SepCl << (_oFile.Count() + TotalSlaveWrReadsCnt);
		if( TestMode ) {
			OutputReadCnt(FG, ", from wich foreground");
			OutputReadCnt(BG, ", background");
		}
		cout  << endl;
	}
}

// Writes each Read started from += RGL_SHIFT positions
bool Imitator::CutRegular	()
{
	short RegShift = RGL_SHIFT();
	Nts nts(_chrFiles.FileName(0), LetN);	// chromosome 1
	ChromsThreads::ChromsThread csThread(true);
	ChromCutter chrCutter(this, &csThread, false);
	Timer timer;

	// Chromosome pnts is regulary cutted with writing to file
	for( ULONG currPos = 0;
		chrCutter.AddRead(nts, currPos, RegShift, false);
		//0);	// never generate reverse strand
		currPos+=RegShift );
	OutputChromName(1, true);
	chrCutter.OutputChromInfo(nts, timer, false);
	timer.Stop();
	return false;
}

// Curs genome into fragments and generate output
bool Imitator::CutGenome	()
{
	BYTE i;
	Array<Thread*> slaves( min(THREADS_CNT(),  _chrFiles.TreatedCount()) - 1) ;
	ChromsThreads cThreads(threadnumb(slaves.Length()+1), _chrFiles);
	
	SetSample();
	if( slaves.Length() && Verbose(vDEBUG))	cThreads.Print();	//return true;

	for(i=0; i<slaves.Length(); i++)			// run slave threads
		slaves[i] = new Thread(StatCutChrom, &cThreads[i+1]);
	retThreadValType res = CutChrom(&cThreads[0], !slaves.Length());	// run main thread
	for(i=0; i<slaves.Length(); i++) {			// wait for slave threads finishing
		slaves[i]->WaitFor();
		delete slaves[i];
	}
	return bool(res);
}

// Sets adjusted Samples and clear all counter and means.
// Samples are needed for the control of BF&FG levels by percent (given by user),
// and to prorate number of written reads for each chromosome depending on reads limit.
void Imitator::SetSample()
{
	float sizeFactor;	// ratio lenth_of_nts / size_of_file
	float commonAvrg;	// lognormal average
	float savedAvrg;	// saved average
	AvrFragLengths avrgs(_chrFiles.Path());

	// *** Get averages from file if it exists, otherwise calculate and save theirs
	if( !avrgs.Get(Amplification::Coefficient, &sizeFactor, &commonAvrg, &savedAvrg) ) {
		Timer timer;
		chrlen pos = 0;
		if( Verbose(vRT) ) {
			cout << "Samples for given distribution do not exist. Generate...";
			fflush(stdout);
			timer.Start();
		}
		Nts nts(_chrFiles.FileName(), true);
		ChromsThreads::ChromsThread cThread(false);	// 'trial' ChromsThread without writing to file
		ChromCutter chrCutter(this, &cThread, true);
		chrCutter.CutChrom(nts, &pos, Region(0, nts.Length()), true);

		sizeFactor = float(nts.Length()) / _chrFiles.FirstFileLength();
		commonAvrg = LognormDistribution::Mean();
		savedAvrg = LognormDistribution::SavedMean(chrCutter._wrReadsCnt[0]);
		avrgs.Add(Amplification::Coefficient, sizeFactor, commonAvrg, savedAvrg);
		if( Verbose(vRT) ) {
			cout << MsgDone;
			timer.Stop();
			//fflush(stdout);
		}
		if( Verbose(vDEBUG) ) {
			cout << SignDbg << "CHROMOSOME " << Chrom::Name(_chrFiles.FirstChromID()) << ":\n";
			cout << SignDbg << "Numbers::\tfrags" << SepCl << LognormDistribution::CallsCnt()
				 << "\treads" << SepCl << chrCutter._wrReadsCnt[0] << EOL;
			cout << SignDbg << "coeff of ampl" << SepCl << Amplification::SimpleMean() << EOL;
		}
	}
	Samples[0] = float(SAMPLE_FG()/100);
	if( TestMode )
		Samples[1] = float(Samples[0] * SAMPLE_BG()/100);
	if( Verbose(vDEBUG) ) {
		cout<< EOL << SignDbg << cAverage << SepCl << "frag" << SepCl << commonAvrg
			<< "\tsaved" << SepCl << savedAvrg << EOL;
		cout << SignDbg << "sizeFactor" << SepCl << sizeFactor << EOL;
	}
	// *** Determine the total possible numbers of saved reads
	ULLONG totalCnt = 0;	// total number of saved reads
	ULLONG FtrsLen;			// length of all features
	ULLONG cnt;
	float countFactor = CellsCnt / savedAvrg;	// coefficient in formula:
												// ntsCount = countFactor * Sample * ntsLen
	
	for(ChromFiles::cIter it=_chrFiles.cBegin(); it!=_chrFiles.cEnd(); it++)
		if( _chrFiles.IsTreated(it) ) {
			//cout << TAB << Chrom::AbbrName(CID(it)) << TAB;
			if( Bed && Bed->FindChrom(CID(it)) ) {
				FtrsLen = Bed->FeaturesTreatLength(CID(it), it->second.Numeric(), commonAvrg);
				// count of foreground Reads
				cnt = (ULONG)(Samples[0] * FtrsLen * countFactor);
				SetMaxDigitCnt(FG, cnt/3);	// 3 just to reduce digits number to 1
				totalCnt += cnt;
				//cout << "fg cnt" << SepCl << cnt << TAB;
			}
			else
				FtrsLen = 0;
			// count of background Reads
			cnt = ULLONG(Samples[1] *
				(_chrFiles.ChromTreatLength(it, sizeFactor) - FtrsLen) * countFactor);
			//cout << "bg cnt"<< SepCl << cnt << EOL;
			SetMaxDigitCnt(TestMode ? BG : FG, cnt);
			totalCnt += cnt;
		}

	if( Verbose(vDEBUG) )
		cout << SignDbg << "CONTROL TOTAL READS TO WRITE" << SepCl << totalCnt << EOL
			 << SignDbg << "MaxDigitCnt: fg = " << (int)DigitsCnt[0]
			 << "\tbg = " << (int)DigitsCnt[1] << EOL;
		
	// *** Estimate adjusted Sample
	if( totalCnt > Read::MaxCount ) {
		AdjSample = (float)Read::MaxCount / totalCnt;
		if( ControlMode )
			AdjSample *= NoAmplification ? 1.03f : 1.03f;	// the best in practise
		if( Verbose(vRES) )
			cout << "Added recovery sample = " << setprecision(2) << (AdjSample * 100)
			<< "% due to reads limit of " << Options::GetDVal(oREAD_LIMIT) << endl;
	}
	TreatOutFtrs = All || Samples[1];
	if( Verbose(vDEBUG) )	cout << endl;
}

// Sets maximal count of digits for given ground if val is maximal
void Imitator::SetMaxDigitCnt(eGround gr, ULLONG val)
{
	BYTE ndigits = DigitsCount(val);
	if( ndigits > DigitsCnt[gr] )
		DigitsCnt[gr] = ndigits;
}

/************************  end of class Imitator ************************/

/************************ class AlterFQ ************************/
//long	AlterFQ::_cntReplSeqs = 0;
//
//void AlterFQ::ReplaceReads(
//	FqFile* inFile, FqFile* outFile,
//	char *extReads, BYTE perc, int eachReplSeq, const BYTE quality)
//{
//	bool isLimited = perc==0;
//	BYTE percent = isLimited ? 10 : perc;
//	BYTE cntInsPercent = 0;	// inserted percent counter (in ten)
//	BYTE cntPercent = 0;		// total percent counter (in ten)
//	size_t	startExtReadPos = 0;		// current index of the beginning of external Read
//
//	char*	recPtr;
//
//	_cntReplSeqs = 0;
//	while( (recPtr = inFile->GetSequence()) != NULL )
//	{
//		if( percent == 10 || cntInsPercent < percent )	// processing percent counters
//		{
//			// replace current Read by external Read if they aren't finished
//			if( inFile->CntSeqs() % eachReplSeq == eachReplSeq-1 
//			&& extReads[startExtReadPos] )
//			{
//				//cout << "\a";
//				inFile->CopyRead(extReads + startExtReadPos);
//				if( quality )
//					inFile->SetQuality(quality);
//				startExtReadPos += Read::Len+1;
//				_cntReplSeqs++;
//			}
//			outFile->AddSequence(recPtr, inFile->RecordLength());
//		}
//		if( cntPercent >= 10)	// dumping percent counters
//		{
//			cntPercent=1;
//			cntInsPercent=0;
//		}
//		else					// increasing percent counters
//		{
//			cntPercent++; 
//			cntInsPercent++;
//		}
//	}
//}

