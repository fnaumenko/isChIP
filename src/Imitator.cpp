/**********************************************************
Imitator.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 21.08.2019
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

void Imitator::ChromCutter::MDA::Split(fraglen shift, fraglen len)
{
	if(len < _minLen)		return;
	_fracs.push_back(make_pair(shift, len));
	if(!Imitator::IsMDA)	return;
	fraglen len1 = fraglen(_rng.Range(len));	// new right fraction relative position
	Split(shift + len1, len - len1);			// split right fraction
	Split(shift, len1 - 1);						// split left fraction
}

// Prepare instance to the new generating cycle
//	@fragLen: current fragment length
//	@fragLenMin: current minimal fragment length
void Imitator::ChromCutter::MDA::Reset(fraglen fragLen, fraglen fragLenMin)
{
	_fracs.clear();
	_minLen = fragLenMin;
	Split(0, fragLen);
	_it = _fracs.begin();
	
	//cout << "total " << _fracs.size() << EOL;
	//for(vector<Frac>::const_iterator it=_fracs.cbegin(); it!=_fracs.cend(); it++)
	//	cout << it->Shift << TAB << it->Len << EOL;
	//cout << EOL;
}

// Gets length and shift of amplificated fraction (MDA)
//  shift: returned shift of fraction
//  return: length of fraction exceeded min frag length, established in Reset(),
//	or 0 if fractions are ended
fraglen Imitator::ChromCutter::MDA::GetFraction	(fraglen* shift)
{
	if(_it == _fracs.end())	return 0;
	*shift = _it->first;
	return _it++->second;
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
public:
	struct AvrFrag {
		int		ReadLen;
		float	SelLen;
		float	RecLen;
		float	MdaLen;

		inline AvrFrag(int readLen=0, float selLen=0, float recLen=0, float mdaLen=0) 
			: ReadLen(readLen), SelLen(selLen), RecLen(recLen), MdaLen(mdaLen)  {}

		// for sorting by ascent
		inline bool operator < (const AvrFrag& avr) const {	return (avr.ReadLen > ReadLen); }
	};
	//static const BYTE	Fl_capacity = 6;	// Number of significant float digits stored

private:
	bool	_isChanged;
	string	_fName;
	vector<AvrFrag> _avrs;
	vector<AvrFrag>::iterator _it;	// iterator to record for current Read length

public:
	// Creates a new instance of AvrFrags from given service path;
	// a file name will be constructed.
	// If file does not exist, the instance is empty.
	AvrFrags(const string& path) : _isChanged(false)
	{
		if(!path.length())	{ _fName = strEmpty; return; }	// empty path means don't save file

		ostringstream ss;
		ss << path << DistrParams::lnMean << HPH << DistrParams::lnSigma;
		if(DistrParams::IsSS())
			ss << HPH << int(DistrParams::ssMean) << HPH << DistrParams::ssSigma;
	#ifdef OS_Windows	
		ss << ".txt";
	#else
		ss << ".ini";
	#endif
	
		TabFile file(_fName = ss.str(), TxtFile::READ_ANY, 4);
		_avrs.reserve(3);		// with a margin
		while(file.GetLine())
			_avrs.push_back(AvrFrag(
				file.IntField(0), file.FloatField(1), file.FloatField(2), file.FloatField(3)));
	}
	
	// Writes instance to a file if it's changed.
	~AvrFrags()
	{
		if(!_isChanged || !_fName.length())		return;

		const char* comms[] = { 
			"[1][2]:unused; generated: [3]:min, [4]:max",
			"[1]:read len; avr: [2]:selected, [3]:recorded; [4]:mda recorded"
		};
		const BYTE commsCnt = sizeof(comms)/sizeof(char*);
		BYTE i = 0;
		ofstream file;

		file.open (_fName.c_str(), ios_base::out);
		file << "# info for sampling generated by isChIP; do not change. Fragment length:\n";
		//for(vector<AvrFrag>::iterator it = _avrs.begin(); it != _avrs.end(); it++) {
		sort(_avrs.begin() + 1, _avrs.end());
		for(_it = _avrs.begin(); _it != _avrs.end(); _it++) {
			file << _it->ReadLen << TAB << _it->SelLen << TAB << _it->RecLen << TAB << _it->MdaLen;
			if(i < commsCnt)	file << TAB << comms[i++];
			file << setprecision(2) << fixed << EOL;
		}
		file.close();
	}

	// Returns records for current Read len
	AvrFrag& Get()
	{
		if(_avrs.size()) {
			for(_it=_avrs.begin() + 1; _it<_avrs.end(); _it++)
				if(_it->ReadLen == Read::Len)
					return *_it;
			_avrs.push_back(AvrFrag());
			_isChanged = true;
			return *(_it = _avrs.end() - 1);
		}
		else {		// empty container
			_avrs.push_back(AvrFrag());
			_avrs.push_back(AvrFrag());
			_isChanged = true;
			return *(_it = _avrs.begin() + 1);
		}
	}

	void SetMinMax(fraglen min, fraglen max)
	{
		AvrFrag& avr = *_avrs.begin();
		if(!avr.RecLen)
			avr.RecLen = min,
			avr.MdaLen = max;
	}

	void SetMDAlen(float len) {	_it->MdaLen = len; _isChanged = true; }

	inline fraglen Min() const { return  _avrs.begin()->RecLen; }

	inline fraglen Max() const { return  _avrs.begin()->MdaLen; }

#ifdef DEBUG
	void Print() const
	{
		for(vector<AvrFrag>::const_iterator it=_avrs.begin(); it<_avrs.end(); it++)
			cout << it->ReadLen << TAB << it->SelLen << TAB << it->RecLen << TAB << it->MdaLen << EOL;
	}
#endif
};

/************************ AvrFrags: end ************************/

// Thread-safety increment sizes by RefSeq
void Imitator::GenomeSizes::IncrSizes(const RefSeq& seq)
{
	Mutex::Lock(Mutex::INCR_SUM);
	Total += seq.Length();
	Defined += seq.DefRegion().Length();
	Gaps += seq.GapLen();
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
BYTE Imitator::ChromView::GapsWexcl;

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
	if(Imitator::IsMDA || AmplCoeff)
		lineW += PrFittedStr("ampl", marg_A + AmplW);	// amplification coefficient 
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
			lineW += PrFittedStr(tGapsExcl, margGexcl + GapsWexcl);
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
	if(Imitator::IsMDA || AmplCoeff)
		PrFittedFloat(false, fragCnt.RealAmplCoeff(), AmplPr, marg_A + AmplW);		// MDA coef
	PrFittedFloat(false, ReadDens(rCnt, densLen), DensPr, marg_D + DensW, DensWRP);	// density
	if(GrType == Gr::FG)		PrintMarg(margD_);
}

// Prints chroms gaps statistics
void Imitator::ChromView::PrintGaps(const GenomeSizes& s)
{
	PrFittedFloat(true, s.GapsInPers(), GapsPr, margD_ + GapsW);
	if(RefSeq::LetGaps)		return;
	PrFittedFloat(true, s.UndefInPers(), GapsPr, margGexcl + strlen(tGapsExcl));
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

	GapsWexcl = strlen(tGapsExcl);
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
	Output::SetSeqMode(avr);
	_fragCnt.Clear();
	_fragCnt.SetGMode(GM::Test);
	_output = slave ? new Output(imitator->_oFile) : &(imitator->_oFile);
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
	if(Verbose(vPAR))
		if(_gMode == GM::Test) {
			const GenomeSizes s(seq);
			ChromView::PrintGaps(s);
		}
		else if(Timer::Enabled)	ChromView::PrintGapsMarg();
	timer.Stop(ChromView::marg_T);		// print time
	if(excLimit) {
		cout << " exceeded limit";
		if(!Verbose(vPAR))	
			cout << " of " << Seq::ReadsLimit() << BLANK << FT::ItemTitle(FT::ABED, true);
	}
	cout << endl;
	
	Mutex::Unlock(Mutex::OUTPUT);
}

ULLONG Imitator::ChromCutter::PrepareCutting(GM::Mode gm, chrid cID, Timer& timer)
{
	SetGMode(gm);
	PrintChromName(cID, gm, IsSingleThread());			// print chrom name before cutting
	timer.Start();
	return CellCnt(gm) << int(_cSizes.IsAutosome(cID));	// multiply twice for autosomes
}

// Treats chromosomes given for current thread
//	@cSubset: pointer to Subset - set of chrom IDs treated in this thread
thrRetValType Imitator::ChromCutter::Execute(const effPartition::Subset& cIDSet)
{
	chrid	cID;
	Features::cIter	cit;			// template chrom's iterator
	ULONG	n, cellCnt;				// count of cells, length of enriched regions
	chrlen	currPos, k, fCnt;		// count of features
	chrlen*	pCurrPos = &currPos;	// pointer to current position
	chrlen	enrRegLen;				// length of enriched regions
	int		res = 0;				// result of cutting
	Timer	timer(Verbose(vRT));	// print local time on Verbose 'runtime info and above'

	try {
		for(effPartition::numb_id_cit it=cIDSet.NumbIDs().begin(); it!=cIDSet.NumbIDs().end(); it++) {
			_fragCnt.Clear();
			_output->BeginWriteChrom(cID = *it);
			cellCnt = PrepareCutting(GM::Test, cID, timer);

			if(Templ && (cit=Templ->GetIter(cID)) != Templ->cEnd()) {
				fCnt = Templ->Count(cit);
				enrRegLen = Templ->EnrRegLength(cit, 0, SelFragAvr);
			}
			else	enrRegLen = fCnt = 0;
			RefSeq seq(cID, _cSizes);
			
			for(n = 0; n < cellCnt; n++) {
				currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
				for(k=0; k < fCnt; k++)
					if(res = CutChrom(seq, pCurrPos, Templ->Feature(cit, k), false))
						goto A;			// achievement of Reads limit
				// add background after last 'end' position
				if((res = CutChrom(seq, pCurrPos, seq.DefRegion(), true)) < 0)
					break;				// achievement of Reads limit
			}
A:			PrintChrom(cID, seq, enrRegLen, timer, res < 0);		// timer stops and printed in here
			IncrTotalSelFragCount();
			// collect total enriched regions length to calculate total density
			IncrementTotalLength(seq, enrRegLen);
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
			if(res < 0)		break;			// achievement of Reads limit
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
//	@fragStart: fragment start position
//	@feature: current treated feature
//	@bg: if true then generate background (swap foreground and background)
//	@fragStat: statistics of selected fragments average counter, or NULL under work mode
//	return: 0 if success,
//		1 if end chromosome is reached (continue treatment),
//		-1 if Reads limit is achieved (cancel treatment)
int Imitator::ChromCutter::CutChrom	(
	const RefSeq& seq,
	chrlen* const pfragStart,
	const Featr& feature,
	bool bg,
	FragLenStat* fragStat
	)
{
	const chrlen cLen = seq.End();	// chrom defined end position
	const float featScore[] = { UniformScore ? 1 : feature.Score, 1 };	// FG, BG feature's score
	bool	select;			// selection by corrected bounds; always true for BG
	bool	primer;			// true for the primer (not amplified) frag
	Gr::Type	g;			// ground: 0 - FG, 1 - BG
	fraglen fragStart,		// fragment's start position
			fragLen,		// fragment's length
			fragLenMin = Read::Len,		// minimal fragment's length after size selection
			fragLenMax = FRAG_MAX,		// maximal fragment's length after size selection
			fracLen,		// fraction's length
			fracShift;		// fraction's start position within fragment
	short	uZone;			// ustable zone, with increasing probability of frag binding
	chrlen	fragEnd;		// fragment's end position
	a_cycle	i;				// PCR amplification counter
	//int	res;				// result of AddRead() method

	for(; *pfragStart<=feature.End; *pfragStart += fragLen) {	// ChIP: control right mark
		fragLen = _fragDistr.LognormNext();
		fragStart = *pfragStart;
		fragEnd = fragStart + fragLen;
		if(fragEnd > cLen)	
			if(fragStart >= cLen - Read::Len)	return 1;	// end of chrom
			else	fragLen = (fragEnd = cLen) - fragStart;	// cut last fragment
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
		//res = 0;
		fracLen = fragLen;
		g = Gr::Type(!bg ^ (fragEnd >= feature.Start));
		
		if(PerSample(g)) {	// ground sampling
			if(!(select = bool(g))					// foreground; for background select always true
			&& (select = fragEnd >= feature.Start)	// fragment captures feature?
			&& FlatLen) {							// flattening is ON?
				uZone = short(fragEnd - feature.Start);		// is frags end in unstable zone?
				if(uZone > FlatLen)
					uZone = short(feature.End - fragStart);	// is frags start in unstable zone?
				if(uZone <= FlatLen)						// is frags start or end in unstable zone?
					select = _fragDistr.Sample(float(uZone) / FlatLen);	// new unstable select
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
							if(!_output->AddRead(seq, fragStart+fracShift, fracLen, g)
							&& IncrRecFragCount(g, primer))		// recorded reads
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
short	Imitator::FlatLen = 0;		// BS edge flattening length
a_coeff	Imitator::AmplCoeff = 0;	// user-stated amplification coefficient
BYTE	Imitator::Verb;
BYTE	Imitator::ThrCnt;			// actual number of threads
bool	Imitator::IsMDA;
bool	Imitator::MakeControl;		// true if control file (input) should be produced
bool	Imitator::UniformScore;
bool	Imitator::All;
eMode	Imitator::TMode;			// Current task mode
float	Imitator::SelFragAvr;		// mean length of selected fragments
Imitator	*Imitator::Imit = NULL;
const Features	*Imitator::Templ = NULL;

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
//	@gMode: Test|Control
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
	cout << Total << " recorded "; Output::PrintItemTitle();
	if(MakeControl)
		cout << SepCl << GM::Title(GM::Test) << COLON;
	cout << BLANK << GlobContext[GM::Test].RecCnt();
	if(MakeControl)	
		cout << SepCm << GM::Title(GM::Control) << SepCl << GlobContext[GM::Control].RecCnt();
	cout << endl;
}

// Print amplification info
void Imitator::PrintAmpl()
{
	cout << SignPar << "Amplification" << SepCl;
	if(IsMDA)	cout << "MDA";
	if(AmplCoeff) {
		if(IsMDA)	cout << SepSCl;
		cout << "PCR cycles" << Equel << int(AmplCoeff);
	}
	else if(!IsMDA)	cout << Options::BoolToStr(false);
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
void Imitator::Execute(Features* templ)
{
	Templ = templ;
	CutGenome();

	// print statistics
	if(Verb == vRES)		PrintTotal();
	else if( Verbose(vRT) )	{
		if(_cSizes.TreatedCount() > 1) {	// print summary test statistics?
			PrintChromInfo(Chrom::UnID, GM::Test, GlobContext[GM::Test].fCnts, TreatedLen);
			if(Verbose(vPAR))	ChrView[Gr::BG].PrintGaps(gSizes);	// ground doesn't matter
			cout << EOL;
		}
		if(TestMode)	PrintTotal();
	}
}

// Cuts genome into fragments and generate output
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

// Imitates cutting chromosome to reach statistics
ULLONG Imitator::CutForSample(Average& genFrAvr, FragLenStat* frLenStat)
{
	chrlen		pos = 0;
	ChromCutter cCutter(this, &genFrAvr, false);
	const RefSeq seq(_cSizes[0]);
	
	// generate statistics based on first chrom; Reads are not recorded
	GlobContext[GM::Test].SetSample(1.0);
	for(int i=0; i<4; pos=0, i++)
		cCutter.CutChrom(seq, &pos, seq.DefRegion(), false, frLenStat);
	GlobContext[GM::Test].ClearFragCounters();
	return cCutter.RecFgFragCnt();
}

// Sets adjusted Sample and clear all counter and means.
// Sample are needed for the control of BF&FG levels by percent (given by user),
// and to prorate number of written reads for each chromosome depending on reads limit.
void Imitator::SetSample()
{
	// *** Get averages from file if it exists, otherwise calculate and save ones
	AvrFrags avrs(_cSizes.ServPath());
	AvrFrags::AvrFrag& avr = avrs.Get();

	if(!avr.ReadLen) {
		Average		genFrAvr;					// generated frags average length
		FragLenStat frLenStat;
		bool isMDA = IsMDA;
		IsMDA = false;
		ULLONG recCnt = CutForSample(genFrAvr, &frLenStat);	// cut without MDA

		avrs.SetMinMax(frLenStat.Min, frLenStat.Max);
		avr.ReadLen = Read::Len;
		avr.SelLen	= frLenStat.SelAvr.Mean();	// selected frags average length
		avr.RecLen	= float(genFrAvr.Sum()) / recCnt;
		IsMDA = isMDA;							// restore MDA
	}
	// *** MDA accounting
	if(IsMDA && !avr.MdaLen) {
		Average		genFrAvr;
		ULLONG recCnt = CutForSample(genFrAvr, NULL);

		avrs.SetMDAlen(float(genFrAvr.Sum()) / recCnt);
	}
	//avrs.Print();
	SelFragAvr = avr.SelLen;
	
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
	const float	countFactor = float(GlobContext[GM::Test].CellCnt) / (IsMDA ? avr.MdaLen : avr.RecLen);

	// *** estimate total number of frags
	for(ChromSizes::cIter it=_cSizes.cBegin(); it!=_cSizes.cEnd(); it++) {
		if( !_cSizes.IsTreated(it) )	continue;
		// count is estimated according to diploid (numerical) sign,
		//	but density not, because basic length is single!
		if( Templ && Templ->FindChrom(CID(it)) ) {
			enRgnLen = Templ->EnrRegLength(CID(it), 0, SelFragAvr);
			totalCnt += GetReadsCnt(Gr::FG, enRgnLen, countFactor, 0, maxCnt, maxDens);
		}
		else	enRgnLen = 0;
		totalCnt += GetReadsCnt(Gr::BG,	_cSizes.DefEffLength(it) - enRgnLen,
			countFactor, _cSizes.IsAutosome(CID(it)), maxCnt, maxDens);
	}
	//if(IsMDA)	totalCnt *= 1.1;	// empirical coefficient
	if(AmplCoeff) {
		ChromCutter::SetAmpl();
		totalCnt *= pow(2.f, int(AmplCoeff));
	}
	// *** Estimate adjusted Sample
	if(totalCnt > Seq::ReadsLimit())
		AutoSample = Seq::ReadsLimit() / totalCnt;
	// *** print debug info
	if(Verbose(vPAR)) {
		cout << SignPar << "Selected fragments size" << SepCl
			 << sActual << "Mean" << Equel << SelFragAvr
			 << SepSCl << "minimum ~ " << avrs.Min()
			 << SepSCl << "maximum ~ " << avrs.Max()
			 << EOL;
		if(MakeControl)
			cout << SignPar << "Generated " << GM::Title(GM::Control) << SepDCl	
				 << "Count of cells" << Equel << CellCnt(GM::Control) << SepSCl
				 << "sample" << Equel << sPercent(Sample(GM::Control, Gr::BG)*100, 2) << EOL;
	}
	if(Verbose(vDBG))
		cout << SignDbg << "Total recorded reads number estimate" << SepCl << totalCnt << EOL;
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
