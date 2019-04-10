/**********************************************************
Imitator.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 10.04.2019
-------------------------
Provides chip-seq imitation functionality
***********************************************************/

#include "isChIP.h"
#include "Imitator.h"
#include <algorithm>    // std::sort
#include <cfloat>		// FLT_MAX
#include <cwchar>		// long '-'

#define PI 3.141593

/************************  class Average ************************/
UINT Imitator::Average::operator+=(UINT val)
{
	static const char* cAverage = "Average";
	if(_keep)
		if( _summator < (ULLONG_MAX-val) ) {
			_summator += val;
			_count++;
			_keep = _count < ULONG_MAX;		//  cut off treatment due to counter overflow
			if( !_keep && Imitator::Verbose(vDBG) )
				Err("block count", cAverage).Throw();
		}
		else {
			_keep = false;					//  cut off treatment due to summator overflow
			if( Imitator::Verbose(vDBG) )	
				Err("block summation", cAverage).Throw();
		}
	return val;
}
/************************  end of class Average ************************/

/************************ FragDistr ************************/

	  float Imitator::ChromCutter::FragDistr::ssFactor0;
const float Imitator::ChromCutter::FragDistr::ssFactor1 = 2.5f / (float)sqrt(2 * PI);
/*
Current release (2 vars defined lognormal distribution):
X_lognorm = exp( Normal() * Sigma + Mean );
where Mean = 5.46, Sigma = 0.4,
Normal(): mean = 0, SD = 1
*/

// Normal distribution with mean=0 and variance=1 (standard deviation = 1)
//	return:  value with Gaussian likelihood between about -5 and +5
//	(from -6 to 6 in 1000000000 cycles) 
double Imitator::ChromCutter::FragDistr::Normal()
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

	w = sqrt( log(w) * (-2./w) );
	_normal_x2 *= w;	// normal_x1 and normal_x2 are independent normally distributed variates
	_phase = 1;
	return normal_x1 * w;	// return normal distributed value
}

/************************ end of FragDistr ************************/

/************************  Amplification ************************/
a_cycle Imitator::ChromCutter::MDA::Coeff = 0;

// Recursively fills _fracts by splitted fragments (segments).
//  @firstCall: true if it is a first call
//	@fragLenMin: current minimum fragment length
void Imitator::ChromCutter::MDA::Split(bool firstCall, fraglen fragLenMin)
{
	fraglen leftFrac, rightFrac;
	_iterCnt = _currCnt;
	for(fraglen i=(_initCnt+_currCnt)/2; i<_iterCnt; i++)	// start from second half of added fractions
	{	// temporary use rightFrac as a Current Fragment Length
		if( firstCall )
			rightFrac = _fragLen;
		else if((rightFrac=_fracts[i].Length) < 0)	// skip "hole"
			continue;
		leftFrac = fraglen(_rng.Range(rightFrac));
		rightFrac -= leftFrac;						// use rightFrac as Right Fraction Length
		if(leftFrac >= fragLenMin) {				// save left fraction on the same position
			_fracts[i].Length = leftFrac;			// shift doesn't change for left fraction
			if(rightFrac >= fragLenMin) {			// save right fraction on the new position
				if(_currCnt+1 >= _memFactor*Coeff) {
					// increase memory for array if necessary
					Fraction* tmp = _fracts;
					_fracts = new Fraction[(++_memFactor) * Coeff];
					memset(_fracts, 0, FractionsSize());		// clear new block
					memcpy(_fracts, tmp, sizeof(struct Fraction)*(_memFactor-1)*Coeff);
					delete [] tmp;
					if( Imitator::Verbose(vDBG) ) {
						Mutex::Lock(Mutex::OUTPUT);
						cout << "Amplification's buffer was increased; new memFactor = "
							 << (int)_memFactor << endl;
						Mutex::Unlock(Mutex::OUTPUT);
					}
				}
				// save right fraction on the new position
				_fracts[_currCnt++].Set(rightFrac, _fracts[i].Shift + leftFrac);	
			}
		}
		else			// discard left fraction
			if( rightFrac >= fragLenMin)
				// save right fraction on the same position
				_fracts[i].Set(rightFrac, _fracts[i].Shift + leftFrac);	
			else
				// discard right fraction: position should be free (a "hole" is created)
				// mark "hole" and clear shift for case saving right fraction in the future
				_fracts[i].Set(-1, 0);	
	}
	if( !firstCall )
		_initCnt = _iterCnt;
	if( _currCnt > _initCnt )
		Split(false, fragLenMin);
}

// Prepare instance to the new generate cycle
//	@fragLen: current fragment length
//	@fragLenMin: current minimal fragment length
void Imitator::ChromCutter::MDA::Reset(fraglen fragLen, fraglen fragLenMin)
{
	_fragLen = fragLen;
	_initCnt = -(Coeff>>1);
	if(_fracts) {
		// fill the array by amplified fractions
		_currCnt = -_initCnt;
		memset(_fracts, 0, FractionsSize());	// clear shifts
		Split(true, fragLenMin);
		_currCnt = -(Coeff>>1) - 1;		// clear current index to prepare for reading
	}
	else {
		_currCnt = -2;
		_iterCnt = 0;
	}
}

// Gets length and shift of amplificated fraction (MDA)
//  shift: returned shift of fraction
//  return: length of fraction exceeded min frag length, established in Reset(),
//	or 0 if fractions are ended
fraglen Imitator::ChromCutter::MDA::GetFraction	(fraglen* shift)
{
	if( ++_currCnt < 0 ) {				// "virtual" part of array: return initial fragment
		_iterCnt++;
		*shift = 0;
		return _fragLen;
	}
	if(_currCnt == _initCnt)	return 0;		// end of records
		
	while(_fracts[_currCnt].Length < 0)			// skip "holes"
		_currCnt++;
	_iterCnt++;
	*shift = _fracts[_currCnt].Shift;
	return _fracts[_currCnt].Length;			// return current fraction
}

/************************  Amplification: end ************************/

/************************ AvrFrags ************************/

// 'AvrFrags' reads and writes average legths of fragment into plain text file.
class AvrFrags
/*
 * Class 'AvrFrags' reads and writes average legths of fragment into plain text file.
 * Each line in file:
 *   field 1:	ampl coefficient or increasing negative int for the first 3 items
 *   field 2:	recorded average frag length for MDA or general value
 *   field 3:	recorded average frag length for PCR or general value
 *   field 4:	comment
 */
{
private:
	struct AvrFrag {
		a_cycle	aCoeff;		// coefficient of amplification
		float	fLen;		// first average length
		float	sLen;		// second average length

		inline AvrFrag(a_cycle coeff, float flen, float slen) : aCoeff(coeff), fLen(flen), sLen(slen) {}

		// for sorting by ascent
		inline bool operator < (const AvrFrag& avr) const {	return (avr.aCoeff > aCoeff); }
	};
	//static const BYTE	Fl_capacity = 6;	// Number of significant float digits stored

	bool	_isChanged;
	string	_fileName;
	vector<AvrFrag> _avrLens;

public:
	static const BYTE SelMIN = 0;				// values array index of selected min frag length
	static const BYTE SelMAX = SelMIN + 1;		// values array index of selected max frag length
	static const BYTE SelAVR = SelMIN + 2;		// values array index of selected average frag length
	static const BYTE GenAVR = SelMIN + 3;		// values array index of generated average frag length
	static const BYTE RecAVR = SelMIN + 4;		// values array index of recorded average frag length
												// for given MDA, PCR ampl coeff
	static const BYTE Capacity = RecAVR + 2;	// values array capacity

	// Creates a new instance of AvrFrags from given path;
	// a file name will be constructed.
	// If file does not exist, the instance is empty.
	AvrFrags	(const string path) : _isChanged(false)
	{
		ostringstream oss;
		oss << path
			<< Options::GetFVal(oFRAG_MEAN)	<< HPH
			<< Options::GetFVal(oFRAG_SIGMA) << HPH;
		if(DistrParams::IsSS())
			oss	<< setprecision(5) << DistrParams::ssMean << HPH
				<< DistrParams::ssSigma;
	#ifdef OS_Windows	
		oss << ".txt";
	#else
		oss << ".ini";
	#endif
	
		TabFile file(_fileName = oss.str(), TxtFile::READ_ANY, 3);
		_avrLens.reserve(RecAVR/2 + 2);	// 2 items for coefficients
		while( file.GetLine() )
			_avrLens.push_back(AvrFrag(a_cycle(file.IntField(0)), 
				file.FloatField(1), file.FloatField(2)));
	}
	
	// Writes instance to a file if it's changed.
	~AvrFrags()
	{
		if(!_isChanged)		return;
		const char* comms[] = { 
			"[1]:unused; selected frag len: [2]:min, [3]:max",
			"[1]:unused; frag avrg len: [2]:selected, [3]:generated",
			"[1]:unused; [2]:recorded frag avrg len; [3]:unused",
			"[1]:ampl coeff; recorded frag avrg len for: [2]:MDA, [3]:PCR"
		};
		const BYTE commsCnt = sizeof(comms)/sizeof(char*);
		BYTE i = 0;
		ofstream file;

		file.open (_fileName.c_str(), ios_base::out);
		file << "# sampling generated by isChIP; do not change\n";
		for(vector<AvrFrag>::iterator it = _avrLens.begin(); it != _avrLens.end(); it++) {
			file << it->aCoeff << TAB << it->fLen << TAB << it->sLen;
			if(i < commsCnt)	file << TAB << comms[i++];
			file << setprecision(2) << fixed << EOL;
		}
		file.close();
	}
	
	// Gets size Factor and fragment length average values.
	//	@isPCR: true if PCR is established
	//	@coeff: user-defined coefficient of amplification
	//	@avrs: extern array to return values; first RecAVR values are always returned,
	//	last 2 only for given coefficient
	//	return: true if values for given coefficient of amplification are returned successfully
	bool Get(bool isPCR, BYTE coeff, float avrs[]) const
	{
		short i = 0;
		for(vector<AvrFrag>::const_iterator it=_avrLens.begin(); it!=_avrLens.end(); it++)
			if(i < RecAVR) {
				avrs[i++] = it->fLen;
				avrs[i++] = it->sLen;
			}
			else if(it->aCoeff == coeff) {
				avrs[i++] = it->fLen;
				avrs[i++] = it->sLen;
				return isPCR ? it->sLen : it->fLen;		// return false if coeff for given ampl is 0
			}
		return false;
	}

	// Adds size Factor and fragment length average values.
	//	@amplCoeff: coefficient pf amplification
	//	@avrs: array of values
	void Add(BYTE coeff, float avrs[])
	{
		if(_avrLens.size()) {	// existed instance
			for(vector<AvrFrag>::iterator it = _avrLens.begin() + RecAVR/2; it != _avrLens.end(); it++)
				if( _isChanged = it->aCoeff == coeff ) {
					it->fLen = avrs[RecAVR];
					it->sLen = avrs[RecAVR+1];
					break;
				}
			if( !_isChanged ) {	// add new coefficient
				_avrLens.push_back(AvrFrag(coeff, avrs[RecAVR], avrs[RecAVR+1]));
				sort(_avrLens.begin() + RecAVR/2, _avrLens.end());
			}
			// do average with new frag average length
			_avrLens[SelAVR/2].fLen = (_avrLens[SelAVR/2].fLen + avrs[SelAVR]) / 2;
		}
		else {					// new instance
			short i = 0;
			for(; i<RecAVR; i+=2)
				_avrLens.push_back(AvrFrag(-1, avrs[i], avrs[i+1]));
			_avrLens.push_back(AvrFrag(coeff, avrs[i], avrs[i+1]));
		}
		_isChanged = true;
	}

#ifdef DEBUG
	void Print() const
	{
		for(vector<AvrFrag>::const_iterator it=_avrLens.begin(); it<_avrLens.end(); it++)
			cout << it->aCoeff << TAB << it->fLen << TAB << it->sLen << EOL;
	}
#endif
};

/************************ AvrFrags: end ************************/

// Thread-safety increment sizes by RefSeq
void Imitator::GenomeSizes::IncrSizes(const RefSeq& seq)
{
	Mutex::Lock(Mutex::INCR_SUM);
	Total = seq.Length();
	Defined = seq.DefRegion().Length();
	Gaps = seq.GapLen();
	Mutex::Unlock(Mutex::INCR_SUM);
}

/************************ ChromView ************************/

/*
	difference between cout & printf precision:

	float v1 = 0.1234, v2 = 0.001234;
	int prec = 2;
	cout << setprecision(prec) << v1 << TAB << v2 << EOL;
	printf("%.*f\t%.*f\n", prec, v1, prec, v2);

	== output
	0.12    0.0012
	0.12    0.00
*/

//#define cFIL_VAL	'*'
#define cFIL_VAL BLANK
const char* Imitator::ChromView::tGapsExcl = "g_excl";
const char* Imitator::ChromView::tTime = "mm:ss";

// Prints str on given field and return field width
//	@width: field width
int PrFittedStr(const char* str, int width)
{
	cout << setw(width) << str;
	return width;
}

// Prints unzero float val or nothing on given field
//	@width: field width
void PrFittedInt(ULONG val, int width)
{
	cout << right << setw(width);
	if(val)	cout << val;
	else	cout << cFIL_VAL;
}

// Prints unzero float value or nothing on given field
//	@percent: if true then print percent
//	@width: field width
//	@padd: right padding - space to the right of the value
void PrFittedFloat(bool percent, float val, int precision, int width, int padd = 0)
{
	if(val)	{
		if(percent)	printf("%s", sPercent(val, precision, width - padd).c_str());
		else		printf("%*.*f", width - padd, precision, val);
		if(padd)	printf("%c", cFIL_VAL);
	}
	else	printf("%*c", width, cFIL_VAL);
}

// Prints chrom's header: Reads
//	@FgHeader: true for foreground header
//	return: header width
int Imitator::ChromView::PrintHeader(bool FgHeader)
{
	int lineW = 0;
	// shift title to right in case of little reads count
	const string& rTitle = FT::ItemTitle(FT::ABED);
	const bool w = CountW < rTitle.length();	

	cout << right;
	// if little reads count then extend "reads" field to 1, otherwhise extend "sample" field to 1
	lineW += PrFittedStr(rTitle.c_str(), 
		GrTitleW() + CountW - w + marg_R*ControlMode);		// "reads"
	lineW += PrFittedStr("sample", marg_S + SampleW);		// "sample"
	// now we should compensate extention to 1:
	// if MDA is set then shrink "MDAc" to 1, otherwhise shrink "unit densuty" field to 1
	if(ChromCutter::IsMDA())
		lineW += PrFittedStr("MDAc", marg_A + AmplW - 1);	// MDA coefficient 
	lineW += PrFittedStr(UnitDens, marg_D + DensW);			// unit density
	if(FgHeader)	PrintMarg(margD_), lineW += margD_;		// empty space
	
	return lineW;
}

// Prints chrom's header: gaps statistics
int Imitator::ChromView::PrintHeaderGaps()
{
	int lineW = 0;

	if(Verbose(vPAR)) {
		lineW = PrFittedStr("gaps", margD_ + GapsW);
		if(!RefSeq::LetGaps )
			lineW += PrFittedStr(tGapsExcl, margG_excl + strlen(tGapsExcl));
	}
	if(Timer::Enabled)
		lineW += PrFittedStr(tTime, marg_T + strlen(tTime));
	return lineW;
}

// Print info about recorded Reads (count, sample, ampl, density) for ground type defined in this nstance
//	@gMode: gen mode for which info is printed
//  @fragCnt: frags counters
//	@densLen: length to print density
void Imitator::ChromView::PrintReads(GM::Mode gMode, const FragCnt fragCnt, chrlen densLen)
{
	ULONG rCnt = ULONG(fragCnt.RecCnt() << Seq::Mode());	// count of reads

	if(TestMode)						// Ground title
		if(gMode == GM::Control)		// can be true only if Imitator::MakeControl is true
			PrintMarg(Gr::TitleLength + 1);					// print blanks instead of Gr title
		else
			cout << right << Gr::Title(GrType) << COLON;	// print Gr title
	PrFittedInt(rCnt, marg_R + CountW);												// count of Reads
	PrFittedFloat(true, fragCnt.Sample(), SamplePr, marg_S + SampleW, SampleWRP);	// sample
	if(ChromCutter::IsMDA())
		PrFittedFloat(false, fragCnt.RealAmplCoeff(), AmplPr, marg_A + AmplW);		// MDA coef
	PrFittedFloat(false, ReadDens(rCnt, densLen), DensPr, marg_D + DensW, DensWRP);	// density
	if(GrType == Gr::FG)		PrintMarg(margD_);
}

// Prints chroms gaps statistics
void Imitator::ChromView::PrintGaps(const GenomeSizes& s)
{
	PrFittedFloat(true, s.GapsInPers(), GapsPr, margD_ + GapsW);
	if(RefSeq::LetGaps)		return;
	PrFittedFloat(true, s.UndefInPers(), GapsPr, margG_excl + strlen(tGapsExcl));
}

// Initializes viewing data
void Imitator::ChromView::Init(ULONG maxCnt, float sample, float maxDens)
{
	//CountW = (BYTE)max( DigitsCount(maxCnt, Options::GetBVal(oLOCALE)) + 1,	
	//	FT::ItemTitle(FT::ABED).length() );		// lengh or "reads"
	CountW = DigitsCount(maxCnt, Options::GetBVal(oLOCALE)) + 1;	// +1 for total
	if(CountW < FT::ItemTitle(FT::ABED).length())
		CountW = FT::ItemTitle(FT::ABED).length();

	//if(sample <= 0.1)		SampleW = 6;	// 0.099% or <0.01%
	//else if(sample <= 1.1)	SampleW = 5;	// 0.99%
	//else if(sample < 10.1)	SampleW = 4;	// 9.9%
	//else if(sample < 100)	SampleW = 3;	// 99%
	//else					SampleW = 4;	// 100%
	SampleWRP = sample * 100 > 9.9;		// '9.9% '
	DensWRP = 0;
	if(maxDens < 0.015)		DensPr = 3;				// 0.0012
	else if(maxDens < 0.15)	DensPr = 3;				// 0.012
	else if(maxDens < 1.5)	DensPr = 2;				// 0.123
	else if(maxDens < 15)	DensPr = 1, DensWRP = 1;// 1.23
	else if(maxDens < 150)	DensPr = 0;				// 123.4
	else					DensPr = 0;				// 1234.5
}

#ifdef DEBUG
void Imitator::ChromView::Print(ULONG maxCnt)
{
	cout << Gr::Title(GrType) << SepCl;
	cout << "MaxCnt = " << setw(6) << setfill(BLANK) << maxCnt << SepCm;
	cout << "CountW = " << int(CountW) << SepCm;
	cout << "SampleW = " << int(SampleW) << SepCm;
	cout << "DensPr = " << int(DensPr) << SepCm;
	cout << "DensW = " << int(DensW) << SepCm;
}
#endif

/************************ ChromView: end ************************/

// pointer to the 'Read counter increment' method
Imitator::FragCnt::pRecIncr	Imitator::FragCnt::pRecIncrSaved;
// pointer to the thread-saved 'selected frag's number adding' method
Imitator::FragCnt::pSelAdd	Imitator::FragCnt::pSelAddSaved;

/************************ class ChromCutter ************************/

fraglen	Imitator::ChromCutter::_SsDev;			// deviation of frag size selection
a_cycle	Imitator::ChromCutter::_PCRdcycles = 1;	// PCR cycles: read doubling cycles
a_coeff	Imitator::ChromCutter::_amplCoeff = 0;	// user-stated amplification coefficient

// Sets global mode
void Imitator::ChromCutter::SetGMode(GM::Mode gmode)
{
	_gMode = gmode;
	_fragCnt.SetGMode(gmode);
	_output->SetGMode(gmode);
}

// Increments counters of local and total recorded fragments thread-safely
//	@g: ground
//	@primer: true if increment derived (amplified) frag's counter
//	return: true if Reads limit is exceeded.
bool Imitator::ChromCutter::IncrRecFragCount(Gr::Type g, bool primer)
{
	_fragCnt[g].RecIncr(primer);	// incr of local recorded Reads
	return GlobContext[_gMode].IncrRecFragCount(g, primer);
}

// Increments counter of total selected fragments thread-safely
void Imitator::ChromCutter::IncrTotalSelFragCount()
{
	for(BYTE i=0; i<Gr::Cnt; i++)
		GlobContext[_gMode].fCnts[i].SelAddSaved(_fragCnt[i].SelCnt());
}

// Set amplification
//	@isPCR: true if PCR set
//	@coeff: amplification coefficient; 0 if no amplification
void Imitator::ChromCutter::SetAmpl(bool isPCR, a_coeff coeff) {
	if(isPCR)	_PCRdcycles = 1<<coeff;
	else		MDA::Coeff = coeff<<1;
	_amplCoeff = coeff;
}

// Creates instance
//	@imitator: the owner
//	@avr: recorded frags average (trial for sampling), or NULL
//	@slave: if true then this instance it slave
Imitator::ChromCutter::ChromCutter(const Imitator* imitator, Average* avr, bool slave) :
	_cSizes(imitator->_cSizes),
	_ampl(_fragDistr),
	_fragDistr(avr),
	_slave(slave),
	_gMode(GM::Test)
{
	OutFiles::SetSeqMode(avr);
	_fragCnt.Clear();
	_fragCnt.SetGMode(GM::Test);
	_output = slave ? new OutFiles(imitator->_oFile) : &(imitator->_oFile);
}

// Prints thread-safe info about treated chroms and stops timer
//	@cID: chrom ID
//  @seq: current reference chromosome
//	@enrRegLen: length of all enriched regions
//	@timer: current timer to thread-saves time output or NULL
//	@exceedLimit: true if limit is exceeded
void Imitator::ChromCutter::PrintChrom (
	chrid cID, const RefSeq& seq, chrlen enrRegLen, Timer& timer, bool excLimit)
{
	if( !Verbose(vRT) )	return;
	const ULONG rgnLens[] = { enrRegLen, seq.DefRegion().Length() - enrRegLen }; // FG, BG region's lengths

	Mutex::Lock(Mutex::OUTPUT);

	PrintChromInfo(cID, _gMode, _fragCnt.GetFragCnts(), rgnLens, !IsSingleThread());
	if(Verbose(vPAR)) {
		const GenomeSizes s(seq);
		ChromView::PrintGaps(s);
	}
	timer.Stop(ChromView::marg_T);		// print time
	if(excLimit) {
		cout << " exceeded limit";
		if(!Verbose(vDBG))	
			cout << " of " << Seq::ReadsLimit() << BLANK << FT::ItemTitle(FT::ABED, true);
	}
	cout << endl;
	
	Mutex::Unlock(Mutex::OUTPUT);
}

ULLONG Imitator::ChromCutter::PrepareCutting(GM::Mode gm, chrid cID, Timer& timer)
{
	SetGMode(gm);
	PrintChromName(cID, gm, IsSingleThread());		// print chrom name before cutting
	timer.Start();
	return CellCnt(gm) << _cSizes.Autosome(cID);	// multiply twice for autosomes
}

// Treats chromosomes given for current thread
//	@cSubset: pointer to Subset - set of chrom IDs treated in this thread
thrRetValType Imitator::ChromCutter::Execute(const effPartition::Subset& cIDSet)
{
	chrid	cID;
	BedF::cIter	cit;				// template chrom's iterator
	ULONG	n, cellCnt;				// count of cells, length of enriched regions
	chrlen	currPos, k, fCnt;		// count of features
	chrlen*	pCurrPos = &currPos;	// pointer to current position
	chrlen	enrRegLen;				// length of enriched regions
	short	lineW = 0;
	Timer	timer(Verbose(vRT));	// print local time on Verbose 'runtime info and above'

	try {
		for(effPartition::numb_id_cit it=cIDSet.NumbIDs().begin(); it!=cIDSet.NumbIDs().end(); it++) {
			cID = *it;
			_fragCnt.Clear();
			cellCnt = PrepareCutting(GM::Test, cID, timer);
			_output->BeginWriteChrom(cID);

			if(Bed && (cit=Bed->GetIter(cID)) != Bed->cEnd()) {
				fCnt = Bed->FeatureCount(cit);
				enrRegLen = Bed->EnrRegLength(cit, 0, SelFragAvr);
			}
			else	enrRegLen = fCnt = 0;
			RefSeq seq(_cSizes.RefName(cID));
			
			for(n = 0; n < cellCnt; n++) {
				currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
				for(k=0; k < fCnt; k++)
					if(lineW = CutChrom(seq, pCurrPos, Bed->Feature(cit, k), false))
						goto A;			// achievement of Reads limit
				// add background after last 'end' position
				if((lineW = CutChrom(seq, pCurrPos, seq.DefRegion(), true)) < 0 && All) {
					_output->EndWriteChrom(cID);
					break;				// achievement of Reads limit
				}
			}
A:			PrintChrom(cID, seq, enrRegLen, timer, lineW < 0);		// timer stops and printed in here
			IncrTotalSelFragCount();
			// collect total enriched regions length to calculate total density
			IncrementTotalLength(seq, enrRegLen);
			if(lineW < 0) {
				_output->EndWriteChrom(cID);
				break;			// achievement of Reads limit
			}
			if(MakeControl) {
				cellCnt = PrepareCutting(GM::Control, cID, timer);
				for(n = 0; n < cellCnt; n++) {
					currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
					CutChrom(seq, pCurrPos, seq.DefRegion(), true);
				}
				PrintChrom(cID, seq, enrRegLen, timer, false);		// timer stops in here
				IncrTotalSelFragCount();
			}
			_output->EndWriteChrom(cID);
		}
	}
	catch(const Err &e)			{ Terminate(cIDSet.ID(), e.what()); }
	catch(const exception &e)	{ Terminate(cIDSet.ID(), e.what()); }
	catch(...)					{ Terminate(cIDSet.ID(), "Unregistered error in thread"); }
	if(!IsSingleThread() && Verbose(vDBG))	{
		Mutex::Lock(Mutex::OUTPUT);
		cout << SignDbg << sThread << int(cIDSet.ID()) << ":  end" << endl;
		Mutex::Unlock(Mutex::OUTPUT);
	}
	return thrRetValFalse;
}

// Cuts chromosome until reaching end position of current treated feature
//	@seq: cutted reference chromosome
//	@currPos: cutting start position
//	@feature: current treated feature
//	@bg: if true then generate background (swap foreground and background)
//	@fragStat: statistics of selected fragments average counter, or NULL under work mode
//	return: 0 if success,
//		1 if end chromosome is reached (continue treatment),
//		-1 if Reads limit is achieved (cancel treatment)
int Imitator::ChromCutter::CutChrom	(
	const RefSeq& seq,
	chrlen* const currPos,
	const Featr& feature,
	bool bg,
	FragLenStat* fragStat
	)
{
	const chrlen cLen = seq.End();	// chrom defined end position
	const readscr featScore[] = { UniformScore ? 1 : feature.Score, 1 };	// FG, BG feature's score
	bool	select;			// selection by corrected bounds; always true for BG
	bool	inDarkZone;
	bool	primer;			// true for the primer (not amplified) frag
	Gr::Type	g;			// ground: 0 - FG, 1 - BG
	fraglen fragLen,		// fragment's length
			fragLenMin = Read::Len,		// minimal fragment's length after size selection
			fragLenMax = FRAG_MAX,		// maximal fragment's length after size selection
			fracLen,		// fraction's length
			fracShift,		// fraction's start position within fragment
			ssLen;			// fragment's size selection random length
	fraglen darkZone;		// uncertain zone, with increasing probability of frag binding
	//fraglen	flatLen = FlatLen > 0 ? FlatLen : 0;	// positive flattering length
	chrlen	fragEnd;		// fragment's end position
	a_cycle	i;				// PCR amplification counter
	int	res;				// result of AddRead() method

	for(; *currPos<=feature.End; *currPos += fragLen) {	// ChIP: control right mark
		fragLen = _fragDistr.LognormNext();
		if(*currPos + fragLen > cLen)	return 1;		// end of chrom
		if(DistrParams::IsSS())
			_fragDistr.SizeSelLimNext(fragLenMin, fragLenMax);
		/*
		 * Since the lower limit of the fragment length after size selection remains unchanged
		 * after amplification, in order to increase the efficiency, it is cut off immediately.
		 * The upper limit is checked after amplification, since the fragment length can decrease with MDA.
		 */
		if(fragLen < fragLenMin)	continue;		// size selection 1: skip short fragment
		/* 
		 * control left mark: 
		 * TestMode: foreground (g==0) is inside and background (g==1) outside template features;
		 * for other chromosomes background (g==0) is inside feature==chrom's length.
		 * ControlMode: foreground (g==0) is always inside feature==chrom's length
		 */
		res = 0;
		fragEnd = *currPos + fragLen;
		fracLen = fragLen;
		g = Gr::Type(!bg ^ (fragEnd >= feature.Start));
		
		if(PerSample(g)) {	// ground sampling
			if(!(select = bool(g))) {		// foreground; for background select always true
				select = fragEnd >= feature.Start;
				if(select && FlatLen) {		// flattening is ON
					darkZone = feature.End - *currPos;
					inDarkZone = darkZone < FlatLen;		// frag's start is in uncertain zone
					if(!inDarkZone) {
						darkZone = fragEnd - feature.Start;
						inDarkZone = darkZone < FlatLen;	// frag's end is in uncertain zone
					}
					if(inDarkZone)	select = _fragDistr.Sample( float(darkZone) / FlatLen );
				}
			}
			if(select && _fragDistr.Sample(featScore[g]))	// selection by bounds & feature score
				// MDA amplification
				for(primer=true, _ampl.Reset(fragLen, fragLenMin); 
					fracLen=_ampl.GetFraction(&fracShift); primer=false)
					if(PerAutoSample()				// adjusted sample
					&& fracLen <= fragLenMax)		// ** size selection 2: skip long fragments
						// PCR amplification
						for(i=0; i<_PCRdcycles; i++) {	// _PCRdcycles is number of read doubling cycles 
							if(i)	primer=false;
							// while BG or without MDA fracShift is always 0
							if(!_output->AddRead(seq, *currPos+fracShift, fragLen, g)
							&& IncrRecFragCount(g, primer))
								return -1;						// Reads limit is exceeded
						}
		}
		// ** size selection 2: statistics record

		if(fragLen <= fragLenMax) {
			_fragCnt[g].SelIncr();
			if(fragStat)	fragStat->TakeFragLen(fragLen);
		}
	}
	return 0;
}

/************************ end of class ChromCutter ************************/

/************************  class Imitator ************************/

Imitator::GenomeSizes Imitator::gSizes;
Imitator::Context	Imitator::GlobContext[2];	// global generation context
Imitator::ChromView Imitator::ChrView[] = {Gr::FG,Gr::BG};
ULONG	Imitator::TreatedLen[] = {0,0};	// FG, BF genome treated length
float	Imitator::AutoSample = 1;		// adjusted FG sample to stay in limit
//readlen	Imitator::BindLen;
short	Imitator::FlatLen = 0;
BYTE	Imitator::Verb;
BYTE	Imitator::ThrCnt;			// actual number of threads
bool	Imitator::MakeControl;		// true if control file (input) should be produced
bool	Imitator::UniformScore;
bool	Imitator::All;
eMode	Imitator::TMode;			// Current task mode
float	Imitator::SelFragAvr;		// mean length of selected fragments
Imitator	*Imitator::Imit = NULL;
const BedF	*Imitator::Bed = NULL;

// Prints chromosome's name
//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
//	@gm: generation mode Test|Control
//	@print: true if chromosomes name should be printed
void Imitator::PrintChromName(chrid cID, GM::Mode gm, bool print)
{
	if( Verbose(vRT) && print) {
		if(MakeControl) 	cout << *GM::Title(gm) << BLANK;
		cout << setw(ChromView::ChromNameW()) << left << setfill(BLANK)
			 << (cID==Chrom::UnID ? Total : Chrom::AbbrName(cID, true)) + COLON;
		fflush(stdout);		// including reset to default right and setfill
	}
}

// Prints chromosome's name and full statistics
//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
//	@fCnts: array of fragment's counters, for FG and BG
//	@rgnLens: array of region's lengths to print FG|BG density
//	@prChrName: true if chromosome's name should be printed
void Imitator::PrintChromInfo(
	chrid cID, GM::Mode gMode, const FragCnt fCnts[], const ULONG rgnLens[], bool prChrName)
{
	PrintChromName(cID, gMode, prChrName);
	if(TestMode)	PrintReadInfo(Gr::FG, gMode, fCnts, rgnLens);
	PrintReadInfo(Gr::BG, gMode, fCnts, rgnLens);
}

// Prints FG and BG, or total, header
///	@header: if false, then print solid line only
void Imitator::PrintHeader(bool header)
{
	if(!Verbose(vRT))	return;
	static int w = 0;

	if(header) {
		w = 0;
		if(Verbose(vPAR))	cout << EOL;
		if(MakeControl)
			w += PrFittedStr(sBLANK, 2);	// "t " or "c "
		cout << setfill(BLANK) << left;
		w += PrFittedStr(Chrom::Short.c_str(), ChromView::ChromNameW());
		if(TestMode)
			w += ChrView[Gr::FG].PrintHeader(true);
		w += ChrView[Gr::BG].PrintHeader(false);
		w += ChromView::PrintHeaderGaps();
		cout << EOL;
	}
	PrintHorLine(w);
}

// Prints total outcome
void Imitator::PrintTotal()
{
	cout << Total << " recorded "; OutFiles::PrintItemTitle();
	if(MakeControl)
		cout << SepCl << GM::Title(GM::Test) << COLON;
	cout << BLANK << GlobContext[GM::Test].RecCnt();
	if(MakeControl)	
		cout << SepCm << GM::Title(GM::Control) << SepCl << GlobContext[GM::Control].RecCnt();
	cout << endl;
}

// Print amplification info
void Imitator::PrintAmpl(const char* ampls[])
{
	cout << SignPar << "Amplification" << SepCl;
	if( ChromCutter::AmplCoeff() )
		cout << ampls[ChromCutter::IsPCR()]
				<< SepSCl << "coefficient" << Equel << int(ChromCutter::AmplCoeff());
	else	cout << Options::BoolToStr(false);
	cout << EOL;
}

// Increments grounds total length.
void Imitator::IncrementTotalLength(const RefSeq& seq, chrlen enrRgnLen)
{
	if(enrRgnLen)	InterlockedExchangeAdd(&(TreatedLen[Gr::FG]), enrRgnLen);
	InterlockedExchangeAdd(&(TreatedLen[Gr::BG]), seq.DefRegion().Length() - enrRgnLen);
	if(Verbose(vPAR))	gSizes.IncrSizes(seq);
}

// Runs task in current mode and write result to output files
void Imitator::Execute(BedF* templ)
{
	Bed = templ;
	CutGenome();

	// print statistics
	if(Verb == vRES)		PrintTotal();
	else if( Verbose(vRT) )	{
		if(TestMode || _cSizes.TreatedCount() > 1) {	// print summary test statistics?
			PrintChromInfo(Chrom::UnID, GM::Test, GlobContext[GM::Test].fCnts, TreatedLen);
			if(Verbose(vPAR))	ChrView[Gr::BG].PrintGaps(gSizes);
			cout << EOL;
		}
		if(TestMode)	PrintTotal();
	}
}

// Curs genome into fragments and generate output
void Imitator::CutGenome	()
{
	effPartition cSets(_cSizes, thrid(ThrCnt));
	vector<Thread*> slaves;
	slaves.reserve(ThrCnt);
	
	if( ThrCnt>1 && Verbose(vDBG))	cSets.Print();
	if(FlatLen < 0)		FlatLen = -FlatLen;
	SetSample();
	PrintHeader(true);
	for(BYTE i=1; i<ThrCnt; i++)	// invoke slave threads
		slaves.push_back(new Thread(StatCutChrom, &cSets[i]));
	CutChrom(cSets[0], false);		// invoke main thread
	for(BYTE i=0; i<ThrCnt-1; i++) {
		slaves[i]->WaitFor();		// wait for slave threads finishing
		delete slaves[i];
	}
	PrintHeader(false);
}

// Gets Reads estimated count per chrom and corrects maximum count and density
//	@g: ground
//	@densLen: length on which density is defined
//	@factor: CellsCnt / recAvr
//	@numeric: 1 for diploid chromosomes, 0 for haploid ones
//	@maxCnt: FG|BG maximum counters to fill
//	@maxDens: FG|BG maximum densities to fill
//	return: estimated number of Reads
ULONG Imitator::GetReadsCnt(Gr::Type g, chrlen densLen, float factor, 
	BYTE numeric, ULONG maxCnt[], float maxDens[])
{
	// Reads count
	ULONG cnt = (ULONG)(Sample(GM::Test, g) * densLen * factor) << Seq::Mode();
	if(maxCnt[g] < cnt)		maxCnt[g] = cnt;
	// Reads density
	factor = ReadDens(cnt, densLen >> numeric);
	if(maxDens[g] < factor)		maxDens[g] = factor;
	return cnt;
}

// Sets adjusted Sample and clear all counter and means.
// Sample are needed for the control of BF&FG levels by percent (given by user),
// and to prorate number of written reads for each chromosome depending on reads limit.
void Imitator::SetSample()
{
	float avr[AvrFrags::Capacity];
	AvrFrags avrgs(_cSizes.Path());

	memset(avr, 0, sizeof(avr));
	// *** Get averages from file if it exists, otherwise calculate and save ones
	if(!avrgs.Get(ChromCutter::IsPCR(), ChromCutter::AmplCoeff(), avr)) {
		chrlen pos = 0;
		Average genFrAvr;					// generated frags average length
		FragLenStat frLenStat;
		ChromCutter cCutter(this, &genFrAvr, false);
		const RefSeq seq(_cSizes);
	
		// generate statistics based on first chrom; Reads are not recorded
		GlobContext[GM::Test].SetSample(1.0);
		for(int i=0; i<2; i++) {
			pos = 0;
			cCutter.CutChrom(seq, &pos, seq.DefRegion(), false, &frLenStat);
		}
		GlobContext[GM::Test].ClearFragCounters();
		
		avr[AvrFrags::SelMIN] = frLenStat.Min;
		avr[AvrFrags::SelMAX] = frLenStat.Max;
		avr[AvrFrags::SelAVR] = frLenStat.SelAvr.Mean();	// selected frags average length
		avr[AvrFrags::GenAVR] = genFrAvr.Mean();			// generated frags average length
		avr[AvrFrags::RecAVR + ChromCutter::IsPCR()] = float(genFrAvr.Sum()) / cCutter.RecFgFragCnt();
		avrgs.Add(ChromCutter::AmplCoeff(), avr);	// avrgs is saved in destructor

		//avrgs.Print();
		//if( Verbose(vDBG) ) {
		//	cout << SignDbg << "Generate by " << Chrom::ShortName(_cFiles.FirstChromID()) << ":\n";
		//	cout << SignDbg << "frags" << SepCl << FragDistr::CallsCnt()
		//		 //<< "  reads" << SepCl << cCutter._recFragCnt[0] << EOL;
		//		 << "  reads" << SepCl << cCutter._fragCnt[0].Rec << EOL;
		//	cout << SignDbg << "coeff of ampl" << SepCl << MDA::SimpleMean() << EOL;
		//}
	}
	SelFragAvr = avr[AvrFrags::SelAVR];
	
	// *** Set Test samples and Control number of cells and sample
	GlobContext[GM::Test].SetSample(SAMPLE_FG()/100);
	if(TestMode)
		if( !(GlobContext[GM::Test].Sample[Gr::BG] *= SAMPLE_BG()/100) )	All = false;
	if(MakeControl)	// can be true only in TEST task mode
		GlobContext[GM::Control].SetControlSample( GlobContext[GM::Test].GetExactBGCellCnt() );

	// *** Determine the total possible numbers of recorded reads
	ULLONG	totalCnt = 0;		// total number of recorded reads
	ULONG	maxCnt[] = {0,0};
	float	maxDens[] = {0,0};
	chrlen	enRgnLen;		// length of enriched regions
	// coefficient in formula: fragCnt_onLen = CellCnt * Sample * Len / recordedAvrLen
	// constant countFactor = CellCnt / recordedAvrLen
	const float	countFactor = float(GlobContext[GM::Test].CellCnt) /
		avr[AvrFrags::RecAVR+ChromCutter::IsPCR()];

	// *** estimate total number of frags
	for(ChromSizes::cIter it=_cSizes.cBegin(); it!=_cSizes.cEnd(); it++) {
		if( !_cSizes.IsTreated(it) )	continue;
		// count is estimated according to diploid (numerical) sign,
		//	but density not, because basic length is single!
		if( Bed && Bed->FindChrom(CID(it)) ) {
			enRgnLen = Bed->EnrRegLength(CID(it), 0, SelFragAvr);
			totalCnt += GetReadsCnt(Gr::FG, enRgnLen, countFactor, 0, maxCnt, maxDens);
		}
		else	enRgnLen = 0;
		totalCnt += GetReadsCnt(Gr::BG,	_cSizes.DefEffLength(it) - enRgnLen,
			countFactor, _cSizes.Autosome(CID(it)), maxCnt, maxDens);
	}
	// *** Estimate adjusted Sample
	if(totalCnt > Seq::ReadsLimit())
		AutoSample = Seq::ReadsLimit() / totalCnt;
	// *** print debug info
	if(Verbose(vPAR)) {
		cout << SignPar << "Selected fragments size" << SepCl
			 << sActual << "mean" << Equel << SelFragAvr
			 << SepSCl << "minimum ~ " << avr[AvrFrags::SelMIN]
			 << SepSCl << "maximum ~ " << avr[AvrFrags::SelMAX]
			 << EOL;
		if(MakeControl)
			cout << SignPar << "Generated " << GM::Title(GM::Control) << SepDCl	
				 << "Count of cells" << Equel << CellCnt(GM::Control) << SepSCl
				 << "sample" << Equel << sPercent(Sample(GM::Control, Gr::BG)*100, 2) << EOL;
	}
	if(Verbose(vDBG))
		cout << SignDbg << "Total recorded reads initial estimate" << SepCl << totalCnt << EOL;
	if(AutoSample < 1 && Verbose(vRES))
		cout << "Added recovery sample = " << setprecision(3) << (AutoSample * 100)
			<< "% due to reads limit of " << Seq::ReadsLimit() << endl;

	// *** set Reads statistics params
	if(TestMode)	
		InitReadsView(Gr::FG, maxCnt, maxDens);
	InitReadsView(Gr::BG, maxCnt, maxDens);
}

/************************  end of Imitator ************************/

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

//const char* HORZ = u8"─";     // HORiZontal line
//const char* VERT = u8"│";     // VERTical line
//const char* TL = u8"┌";     // Top Left
//const char* TR = u8"┐";     // Top Right
//const char* BL = u8"└";     // Bottom Left
//const char* BR = u8"┘";     // Bottom Right
//setlocale(LC_CTYPE,"");
//wchar_t z = L'\304';		// ASCII: gorizontal line
//wchar_t z = L'\315';		// ASCII: gorizontal double line
//wchar_t z = L'\u2500';	// Unicode: thick line
//wchar_t z = L'\u2015';	// Unicode: thin line
