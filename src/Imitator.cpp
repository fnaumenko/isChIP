/**********************************************************
Imitator.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 25.03.2021
-------------------------
Provides chip-seq imitation functionality
***********************************************************/

#include "isChIP.h"
#include "Imitator.h"
#include <algorithm>    // std::sort
#include <cfloat>		// FLT_MAX
#include <cwchar>		// long '-'
#include <random>		// std::exponential_distribution

//#define PI 3.141593

/************************  class Average ************************/
UINT Imitator::Average::operator+=(UINT val)
{
	static const char* cAverage = "Average";
	if(_keep)
		if( _summator < (ULLONG_MAX-val) ) {
			_summator += val;
			_count++;
			_keep = _count < ULONG_MAX;		//  cut off treatment due to counter overflow
			if( !_keep && Imitator::Verbose(eVerb::DBG) )
				Err("block count", cAverage).Throw();
		}
		else {
			_keep = false;					//  cut off treatment due to summator overflow
			if( Imitator::Verbose(eVerb::DBG) )
				Err("block summation", cAverage).Throw();
		}
	return val;
}
/************************  end of class Average ************************/

/************************ FragDistr ************************/

/*
Current release (2 vars defined lognormal distribution):
X_lognorm = exp( Normal() * Sigma + Mean );
where Mean = 5.46, Sigma = 0.4,
Normal(): mean = 0, SD = 1
*/

/************************ end of FragDistr ************************/

/************************  Amplification ************************/

void Imitator::ChromCutter::MDA::Split(fraglen shift, fraglen len, fraglen minLen)
{
	if(len < minLen)		return;
	emplace_back(shift, len);
	fraglen len1 = fraglen(_rng.Range(len));	// new right fraction relative position
	Split(shift, len1 - 1, minLen);				// split left fraction
	Split(shift + len1, len - len1, minLen);	// split right fraction
}

// Generats amplified fragment collection
//	@fLen: current fragment length
//	@fLenMin: current minimal fragment length
void Imitator::ChromCutter::MDA::Generate(fraglen fLen, fraglen fLenMin)
{
	clear();
	if (fLen >= fLenMin)
		if (IsMDA)	Split(0, fLen, fLenMin);
		else 		emplace_back(0, fLen);		// original fragment
	
	//cout << "MDA " << size() << LF;
	//for(Fraction f : *this)	cout << f.first << TAB << f.second << LF;
}

/************************  Amplification: end ************************/

/************************ AvrFrags ************************/

// 'AvrFrags' reads and writes average legths of fragment into plain text file.
class AvrFrags
/*
 *	Class 'AvrFrags' reads and writes average legths of fragment into plain text file.
 *	Each line in file:
 *   field 1:	const Read length
 *   field 2:	selected average frag length
 *   field 3:	recorded average frag length without MDA
 *   field 4:	recorded average frag length with MDA, or 0 if MDA hahas never been applied
 *	Last value (field 4) depends on Read length, therefore file keeps separate line for each different Read length
 *	
 *	Data for variable Read length  is not stored.
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
		ss << FT::Ext(FT::eType::INI);
	
		TabFile file(_fName = ss.str(), FT::eType::INI, TxtFile::eAction::READ_ANY);
		_avrs.reserve(3);		// with a margin
		while(file.GetNextLine())
			_avrs.emplace_back(
				file.IntField(0), file.FloatField(1), file.FloatField(2), file.FloatField(3));
	}
	
	// Writes instance to a file if it's changed.
	~AvrFrags()
	{
		if(!_isChanged || !_fName.length())		return;

		const char* comms[] = { 
			"[1][2]:unused; generated: [3]:min, [4]:max",
			"[1]:read len; avr: [2]:selected, [3]:recorded; [4]:mda recorded for given read len"
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
			file << setprecision(2) << fixed << LF;
		}
		file.close();
	}

	// Returns records for current Read len
	AvrFrag& Get()
	{
		if(_avrs.size()) {
			for(_it=_avrs.begin() + 1; _it<_avrs.end(); _it++)
				if(_it->ReadLen == Read::FixedLen)
					return *_it;
			//_avrs.push_back(AvrFrag());
			_avrs.emplace_back();
			_isChanged = true;
			return *(_it = _avrs.end() - 1);
		}
		else
		{		// empty container
			_avrs.emplace_back();
			_avrs.emplace_back();
			_isChanged = true;
			return *(_it = _avrs.begin() + 1);
		}
	}

	// Returns min frag len from the first line if it's not empty, otherwise 0
	inline chrlen GetMin() const { return _avrs.size() ? _avrs.begin()->RecLen : 0; }

	// Returns max frag len from the first line if it's not empty, otherwise 0
	inline chrlen GetMax() const { return _avrs.size() ? _avrs.begin()->MdaLen : 0; }

	// Writes min, max values to the second line
	void SetMinMax(fraglen min, fraglen max)
	{
		AvrFrag& avr = *_avrs.begin();
		if(!avr.RecLen)
			avr.RecLen = float(min),
			avr.MdaLen = float(max);
	}

	void SetMDAlen(float len) {	_it->MdaLen = len; _isChanged = true; }

#ifdef DEBUG
	void Print() const
	{
		for(vector<AvrFrag>::const_iterator it=_avrs.begin(); it<_avrs.end(); it++)
			cout << it->ReadLen << TAB << it->SelLen << TAB << it->RecLen << TAB << it->MdaLen << LF;
	}
#endif
};

/************************ AvrFrags: end ************************/

// Thread-safety increment sizes by RefSeq
void Imitator::GenomeSizes::IncrSizes(const RefSeq& seq)
{
	Mutex::Lock(Mutex::eType::INCR_SUM);
	Total += seq.Length();
	Defined += seq.DefRegion().Length();
	Gaps += seq.GapLen();
	Mutex::Unlock(Mutex::eType::INCR_SUM);
}

/************************ ChromView ************************/

/*
	difference between cout & printf precision:

	float v1 = 0.1234, v2 = 0.001234;
	int prec = 2;
	cout << setprecision(prec) << v1 << TAB << v2 << LF;
	printf("%.*f\t%.*f\n", prec, v1, prec, v2);

	== output
	0.12    0.0012
	0.12    0.00
*/

//#define cFIL_VAL	'*'
#define cFIL_VAL SPACE
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
	const string& rTitle = FT::ItemTitle(FT::eType::ABED);
	const bool w = CountW < rTitle.length();	

	cout << right;
	// if little reads count then extend "reads" field to 1, otherwhise extend "sample" field to 1
	lineW += PrFittedStr(rTitle.c_str(), 
		GrTitleW() + CountW - w + marg_R*ControlMode);		// "reads"
	lineW += PrFittedStr("sample", marg_S + SampleW);		// "sample"
	// now we should compensate extention to 1:
	// if MDA is set then shrink "MDAc" to 1, otherwhise shrink "unit densuty" field to 1
	if(Imitator::IsMDA || PCRCoeff)
		lineW += PrFittedStr("ampl", marg_A + AmplW);	// amplification coefficient 
	lineW += PrFittedStr(sUnitDens, marg_D + DensW);			// unit density
	if(FgHeader)	PrintMarg(margD_), lineW += margD_;		// empty space
	
	return lineW;
}

// Prints chrom's header: gaps statistics
int Imitator::ChromView::PrintHeaderGaps()
{
	int lineW = 0;

	if(Verbose(eVerb::PAR)) {
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
void Imitator::ChromView::PrintReads(GM::eMode gMode, const FragCnt fragCnt, chrlen densLen)
{
	ULONG rCnt = ULONG(fragCnt.RecCnt() << Seq::Mode());	// count of reads

	if(TestMode)						// Ground title
		if(gMode == GM::eMode::Control)		// can be true only if Imitator::MakeControl is true
			PrintMarg(Gr::TitleLength + 1);					// print blanks instead of Gr title
		else
			cout << right << Gr::Title(GrType) << COLON;	// print Gr title
	PrFittedInt(rCnt, marg_R + CountW);												// count of Reads
	PrFittedFloat(true, fragCnt.Sample(), SamplePr, marg_S + SampleW, SampleWRP);	// sample
	if(Imitator::IsMDA || PCRCoeff)
		PrFittedFloat(false, fragCnt.RealAmplCoeff(), AmplPr, marg_A + AmplW);			// MDA coef
	PrFittedFloat(false, LinearDens(rCnt, densLen), DensPr, marg_D + DensW, DensWRP);	// density
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
	if(CountW < FT::ItemTitle(FT::eType::ABED).length())
		CountW = BYTE(FT::ItemTitle(FT::eType::ABED).length());

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

	GapsWexcl = BYTE(strlen(tGapsExcl));
}

#ifdef DEBUG
void Imitator::ChromView::Print(ULONG maxCnt)
{
	cout << Gr::Title(GrType) << SepCl;
	cout << "MaxCnt = " << setw(6) << setfill(SPACE) << maxCnt << SepCm;
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
void Imitator::ChromCutter::SetGMode(GM::eMode gmode)
{
	_gMode = gmode;
	_fragCnt.SetGMode(gmode);
	_output->SetGMode(gmode);
}

// Increments counters of local and total recorded fragments thread-safely
//	@g: ground
//	@primer: true if increment derived (amplified) frag's counter
//	return: true if Reads limit is exceeded.
bool Imitator::ChromCutter::IncrRecFragCount(Gr::eType g, bool primer)
{
	_fragCnt[g].RecIncr(primer);	// incr of local recorded Reads
	return GlobContext[int(_gMode)].IncrRecFragCount(g, primer);
}

// Increments counter of total selected fragments thread-safely
void Imitator::ChromCutter::IncrTotalSelFragCount()
{
	for(BYTE i=0; i<Gr::Cnt; i++)
		GlobContext[int(_gMode)].fCnts[i].SelAddSaved(_fragCnt[i].SelCnt());
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
	_gMode(GM::eMode::Test)
{
	Output::SetSeqMode(avr);
	_fragCnt.Clear();
	_fragCnt.SetGMode(GM::eMode::Test);
	_output = slave ? new Output(imitator->_oFile) : &(imitator->_oFile);
}

// Prints thread-safe info about treated chroms and stops timer
//  @seq: current reference chromosome
//	@enrRegLen: length of all enriched regions
//	@timer: current timer to thread-saves time output or NULL
//	@exceedLimit: true if limit is exceeded
void Imitator::ChromCutter::PrintChrom (
	const RefSeq& seq, chrlen enrRegLen, Timer& timer, bool excLimit)
{
	if( !Verbose(eVerb::RT) )	return;
	const ULONG rgnLens[] = { enrRegLen, seq.DefRegion().Length() - enrRegLen }; // FG, BG region's lengths

	Mutex::Lock(Mutex::eType::OUTPUT);

	PrintChromInfo(seq.ID(), _gMode, _fragCnt.GetFragCnts(), rgnLens, !IsSingleThread());
	if(Verbose(eVerb::PAR))
		if(_gMode == GM::eMode::Test) {
			const GenomeSizes s(seq);
			ChromView::PrintGaps(s);
		}
		else if(Timer::Enabled)	ChromView::PrintGapsMarg();
	timer.Stop(ChromView::marg_T);		// print time
	if(excLimit) {
		cout << " exceeded limit";
		if(!Verbose(eVerb::PAR))
			cout << " of " << Seq::ReadsLimit() << SPACE << FT::ItemTitle(FT::eType::ABED, true);
	}
	cout << endl;
	
	Mutex::Unlock(Mutex::eType::OUTPUT);
}

ULLONG Imitator::ChromCutter::PrepareCutting(GM::eMode gm, chrid cID, Timer& timer)
{
	SetGMode(gm);
	PrintChromName(cID, gm, IsSingleThread());			// print chrom name before cutting
	timer.Start();
	return ULLONG(CellCnt(gm) << int(_cSizes.IsAutosome(cID)));	// multiply twice for autosomes
}

// Treats chromosomes given for current thread
//	@cSubset: pointer to Subset - set of chrom IDs treated in this thread
thrRetValType Imitator::ChromCutter::Execute(const effPartition::Subset& cIDSet)
{
	Features::cIter	cit;			// template chrom's iterator
	ULONG	n, cellCnt;				// count of cells, length of enriched regions
	chrlen	currPos, k, fCnt;		// count of features
	chrlen	enrRegLen;				// length of enriched regions
	int		res = 0;				// result of cutting
	Timer	timer(Verbose(eVerb::RT));	// print local time on Verbose 'runtime info and above'

	try {
		for(effPartition::numb_id_cit it=cIDSet.NumbIDs().begin(); it!=cIDSet.NumbIDs().end(); it++) {
			_fragCnt.Clear();
			cellCnt = ULONG(PrepareCutting(GM::eMode::Test, *it, timer));

			if(Templ && (cit=Templ->GetIter(*it)) != Templ->cEnd()) {
				fCnt = Templ->ItemsCount(cit);
				enrRegLen = Templ->EnrRegLength(cit, 0, SelFragAvr);
			}
			else	enrRegLen = fCnt = 0;
			const RefSeq seq(*it, _cSizes);
			const chrlen cLen = seq.End();		// chrom 'end' position
			_output->BeginWriteChrom(seq);

			for(n = 0; n < cellCnt; n++) {
				currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
				for(k=0; k < fCnt; k++)
					if(res = CutChrom(cLen, currPos, Templ->Feature(cit, k), false))
						goto A;			// achievement of Reads limit
				// add background after last 'end' position
				if((res = CutChrom(cLen, currPos, seq.DefRegion(), true)) < 0)
					break;				// achievement of Reads limit
			}
A:			PrintChrom(seq, enrRegLen, timer, res < 0);		// timer stops and printed in here
			IncrTotalSelFragCount();
			// collect total enriched regions length to calculate total density
			IncrementTotalLength(seq, enrRegLen);
			if(MakeControl) {
				cellCnt = ULONG(PrepareCutting(GM::eMode::Control, seq.ID(), timer));
				for(n = 0; n < cellCnt; n++) {
					currPos = seq.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
					CutChrom(cLen, currPos, seq.DefRegion(), true);
				}
				PrintChrom(seq, enrRegLen, timer, false);		// timer stops in here
				IncrTotalSelFragCount();
			}
			_output->EndWriteChrom();
			if(res < 0)		break;			// achievement of Reads limit
		}
	}
	catch(const Err &e)			{ Terminate(cIDSet.ID(), e.what()); }
	catch(const exception &e)	{ Terminate(cIDSet.ID(), e.what()); }
	catch(...)					{ Terminate(cIDSet.ID(), "Unregistered error in thread"); }
	if(!IsSingleThread() && Verbose(eVerb::DBG))	{
		Mutex::Lock(Mutex::eType::OUTPUT);
		cout << SignDbg << sThread << int(cIDSet.ID()) << ":  end" << endl;
		Mutex::Unlock(Mutex::eType::OUTPUT);
	}
	return thrRetValFalse;
}

typedef pair<fraglen, fraglen>	frag;

// Returns sample of Flattening of binding site suburb
//	@fStart: tested frag's start
//	@fEnd: tested frag's end
//	@ft: binding site
//	@sample: corrected sample
void Imitator::ChromCutter::GetFlattSample(chrlen fStart, chrlen fEnd, const Featr& ft, bool& sample)
{
	short uZone = short(fEnd - ft.Start);	// ustable zone, 
											// with increasing probability of frag binding
	if (uZone > FlatLen)					// is frags end in unstable zone?
		uZone = short(ft.End - fStart);		// is frags start in unstable zone?
	if (uZone <= FlatLen)					// is frags start or end in unstable zone?
		sample = _fragDistr.Sample(float(uZone) / FlatLen);	// new unstable select
}

// Cuts chromosome until reaching end position of current treated feature
//	@cLen: chromosome's 'end' position
//	@fStart: fragment start position
//	@ft: current treated feature
//	@bg: if true then generate background (swap foreground and background)
//	@fStat: statistics of selected fragments average counter, or NULL under work mode
//	return: 0 if success,
//		1 if end chromosome is reached (continue treatment),
//		-1 if Reads limit is achieved (cancel treatment)
int Imitator::ChromCutter::CutChrom	(
	chrlen cLen,
	chrlen& fStart,
	const Featr& ft,
	bool bg,
	FragLenStat* fStat
	)
{
	const float ftScore[] = { UniScore ? 1 : ft.Value, 1 };	// FG, BG feature's score
	fraglen	fLenMin = Read::FixedLen;	// minimal fragment's length after size selection
	fraglen	fLenMax = FRAG_MAX;			// maximal fragment's length after size selection
	Region frags[4];		// [0],[1] - sheared fragments of forward and backward strand DNA,
							// [2],[3] - additional fragments generated by EXO

	for(fraglen	fLen = 0; fStart<=ft.End; fStart += fLen + 1) {	// ChIP: control right mark
		fLen = _fragDistr.LognormNext();
		chrlen fEnd = fStart + fLen;						// fragment's end position
		if(fEnd > cLen)	
			if(fStart >= cLen - Read::FixedLen)	return 1;	// end of chrom
			else	fLen = (fEnd = cLen) - fStart;			// cut last fragment
		
		if(DistrParams::IsSS())								// get next size seletion limits?
			_fragDistr.SizeSelLimits(fLenMin, fLenMax);
		/*
		 * Since the lower limit of the fragment length after size selection remains unchanged
		 * after amplification, in order to increase the efficiency, it is cut off immediately.
		 * The upper limit is checked after amplification, since the fragment length can decrease with MDA.
		 */

		 //== size selection check 1: skip short fragment
		if(fLen < fLenMin)	continue;
		/* 
		 * control left mark: 
		 * TestMode: foreground (g==0) is inside and background (g==1) outside template features;
		 * for other chromosomes background (g==0) is inside feature==chrom's length.
		 * ControlMode: foreground (g==0) is always inside feature==chrom's length
		 */
		Gr::eType g = Gr::eType(!bg ^ (fEnd >= ft.Start));

		//== EXO processing
		frags[0].Start = frags[1].Start = fStart;
		frags[0].End = frags[1].End = fEnd;
		if (IsExo) {
			int diff = ft.Start - fStart - _fragDistr.Expo();	// left difference
			if (diff > 0)	frags[0].Start += diff;
			diff = fEnd - ft.End - _fragDistr.Expo();			// right difference
			if (diff > 0)	frags[1].End -= diff;
			frags[2] = frags[1];
			frags[3] = frags[0];
		}

		//== sequencing
		bool select = bool(g);				// selection by corrected bounds
		if (!select							// foreground; always true for BG
		&& (select = fEnd >= ft.Start)		// fragment captures feature?
		&& FlatLen)							// flattening is ON?
			GetFlattSample(fStart, fEnd, ft, select);
			
		if (select && _fragDistr.Sample(ftScore[g]))	// selection by bounds & feature score
			for (int x = (2 << IsExo) - 1; x >= 0; x--)	// loop through the fragments from frags[]: 2 or 4
				if (PerSample(g) && frags[x].Length() >= Read::FixedLen) {	// FG/BG loss && not short fragment after EXO
					// ** MDA amplification
					bool primer = true;
					_ampl.Generate(frags[x].Length(), fLenMin);
					for (const Fraction& frac : _ampl)
						if (PerAutoSample()				// adjusted sample
						&& frac.second <= fLenMax)		// ** size selection 2: skip long fragments
							// ** PCR amplification: _PCRdcycles is number of read doubling cycles
							for (a_cycle i = 0; i < _PCRdcycles; primer = false, i++)
								// while BG or without MDA frac.first is always 0
								if (!_output->AddRead(frags[x].Start + frac.first, frac.second, x % 2)
								&& IncrRecFragCount(g, primer))		// recorded reads
									return -1;						// Reads limit is exceeded
				}
		//== size selection check 2: statistics record
		if(fLen <= fLenMax) {
			_fragCnt[g].SelIncr();
			if(fStat)	fStat->TakeFragLen(fLen);
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
a_coeff	Imitator::PCRCoeff = 0;		// user-stated amplification coefficient
eVerb	Imitator::Verb;
BYTE	Imitator::ThrCnt;			// actual number of threads
bool	Imitator::IsExo;
bool	Imitator::IsMDA;
bool	Imitator::MakeControl;		// true if control file (input) should be produced
bool	Imitator::UniScore;
bool	Imitator::All;
eMode	Imitator::TMode;			// Current task mode
float	Imitator::SelFragAvr;		// mean length of selected fragments
Imitator *Imitator::Imit = NULL;
const Features *Imitator::Templ = NULL;

// Prints chromosome's name
//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
//	@gm: generation mode Test|Control
//	@print: true if chromosomes name should be printed
void Imitator::PrintChromName(chrid cID, GM::eMode gm, bool print)
{
	if( Verbose(eVerb::RT) && print) {
		if(MakeControl) 	cout << *GM::Title(gm) << SPACE;
		cout << setw(ChromView::ChromNameW()) << left << setfill(SPACE)
			 << (cID==Chrom::UnID ? sTotal : Chrom::AbbrName(cID, true)) + COLON;
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
	chrid cID, GM::eMode gMode, const FragCnt fCnts[], const ULONG rgnLens[], bool prChrName)
{
	PrintChromName(cID, gMode, prChrName);
	if(TestMode)	PrintReadInfo(Gr::FG, gMode, fCnts, rgnLens);
	PrintReadInfo(Gr::BG, gMode, fCnts, rgnLens);
}

// Prints FG and BG, or total, header
///	@header: if false, then print solid line only
void Imitator::PrintHeader(bool header)
{
	if(!Verbose(eVerb::RT))	return;
	static int w = 0;	// width

	if(header) {
		w = 0;
		if(Verbose(eVerb::PAR))	cout << LF;
		if(MakeControl)
			w += PrFittedStr(sBLANK, 2);	// "t " or "c "
		cout << setfill(SPACE) << left;
		w += PrFittedStr(Chrom::Short.c_str(), ChromView::ChromNameW());
		if(TestMode)
			w += ChrView[Gr::FG].PrintHeader(true);
		w += ChrView[Gr::BG].PrintHeader(false);
		w += ChromView::PrintHeaderGaps();
		cout << LF;
	}
	PrintHorLine(w);
}

// Prints total outcome
void Imitator::PrintTotal()
{
	cout << sTotal << " recorded "; Output::PrintItemTitle();
	if(MakeControl)		// add "test:"
		cout << SepCl << GM::Title(GM::eMode::Test) << COLON;
	cout << SPACE;
	Output::PrintItemCount( GlobContext[int(GM::eMode::Test)].RecCnt() );
	if (MakeControl) {	// add "control:"
		cout << SepCm << GM::Title(GM::eMode::Control) << SepCl;
		Output::PrintItemCount( GlobContext[int(GM::eMode::Control)].RecCnt() );
	}
	cout << endl;		// flash cout buffer
}

// Print amplification info
//	@signOut: output marker
void Imitator::PrintAmpl(const char* signOut)
{
	cout << signOut << "Amplification" << SepCl;
	if(IsMDA)	cout << "MDA";
	if(PCRCoeff) {
		if(IsMDA)	cout << SepSCl;
		cout << "PCR cycles" << Equel << int(PCRCoeff);
	}
	else if(!IsMDA)	cout << Options::BoolToStr(false);
	cout << LF;
}

// Increments grounds total length.
void Imitator::IncrementTotalLength(const RefSeq& seq, chrlen enrRgnLen)
{
	if(enrRgnLen)	InterlockedExchangeAdd(&(TreatedLen[Gr::FG]), enrRgnLen);
	InterlockedExchangeAdd(&(TreatedLen[Gr::BG]), seq.DefRegion().Length() - enrRgnLen);
	if(Verbose(eVerb::PAR))	gSizes.IncrSizes(seq);
}

// Initializes static values
void	Imitator::Init(
	eMode	tmode,				// task mode
	bool	input,				// true if control should be generated as well
	ULONG	cellsCnt,			// count of cells
	bool	isExo,				// true if EXO mode
	bool	isReadLenAssigned,	// true if Read length is assigned by user
	bool	isMDA,				// true if MDA is assigned
	a_coeff	amplCoeff,			// coefficient of PCR
	UINT	verb,				// verbosity level
	bool	allBg,				// true if all background mode is assigned
	bool	uniScore,			// true if uniform template score is assigned
	//readlen bindLen,
	//const pairVal& flattens
	UINT	unstBindLen			// unstable binding length
) {
	TMode = tmode;
	MakeControl = TestMode ? input : false;
	GlobContext[int(GM::eMode::Test)].CellCnt = cellsCnt;
	IsExo = isExo;
	IsMDA = isMDA;
	PCRCoeff = amplCoeff;	// the actual ChromCutter ampl coeff will be set in Sample()
	ChromCutter::FragDistr::Init();
	Verb = eVerb(verb);
	All = (tmode == CONTROL) || allBg;
	UniScore = uniScore;
	//BindLen = bindLen;
	//FlatLen = flattens.first + flattens.second;
	FlatLen = unstBindLen;
	SelFragAvr = DistrParams::LnMean();	// for the Sample(), before set actual ln mean
	if (DistrParams::IsRVL())			// Read variable length mode is set
		if (!isReadLenAssigned)
			Read::FixedLen = Read::VarMinLen;
		else if (Read::FixedLen > Read::VarMinLen && Verb >= eVerb::CRIT)
			Err(string("assigned read length exceeds the default value of ")
				+ to_string(Read::VarMinLen) + " when read variable mode is activated!").Warning();
}

// Runs task in current mode and write result to output files
void Imitator::Execute(Features* templ)
{
	Templ = templ;
	CutGenome();

	// print statistics
	if(Verb == eVerb::RES)		PrintTotal();
	else if( Verbose(eVerb::RT) )	{
		if(_cSizes.TreatedCount() > 1) {	// print summary test statistics?
			PrintChromInfo(Chrom::UnID, GM::eMode::Test, GlobContext[int(GM::eMode::Test)].fCnts, TreatedLen);
			if(Verbose(eVerb::PAR))	ChrView[Gr::BG].PrintGaps(gSizes);	// ground doesn't matter
			cout << LF;
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
	
	if( ThrCnt>1 && Verbose(eVerb::DBG))	cSets.Print();
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
ULONG Imitator::GetReadsCnt(Gr::eType g, chrlen densLen, float factor, 
	BYTE numeric, ULONG maxCnt[], float maxDens[])
{
	// Reads count
	ULONG cnt = (ULONG)(Sample(GM::eMode::Test, g) * densLen * factor) << Seq::Mode();
	if(maxCnt[g] < cnt)		maxCnt[g] = cnt;
	// Reads density
	factor = LinearDens(cnt, densLen >> numeric);
	if(maxDens[g] < factor)		maxDens[g] = factor;
	return cnt;
}

// Imitates cutting chromosome to reach statistics
ULLONG Imitator::CutForSample(Average& genFrAvr, FragLenStat* fLenStat)
{
	chrlen		pos = 0;
	ChromCutter cCutter(this, &genFrAvr, false);
	const RefSeq seq(_cSizes[0]);
	const chrlen cLen = seq.End();		// chrom defined 'end' position
	
	// generate statistics based on first chrom; Reads are not recorded
	GlobContext[int(GM::eMode::Test)].SetSample(1.0);
	for(int i=0; i<2; pos=0, i++)
		cCutter.CutChrom(cLen, pos, seq.DefRegion(), false, fLenStat);
	GlobContext[int(GM::eMode::Test)].ClearFragCounters();
	return cCutter.RecFgFragCnt();
}

// Sets adjusted Sample and clear all counter and means.
// Sample are needed for the control of BF&FG levels by percent (given by user),
// and to prorate number of written reads for each chromosome depending on reads limit.
void Imitator::SetSample()
{
	/***
	If Read constant length is set,
	get averages from file if it exists, otherwise calculate and save ones.
	If Read variable length is set, don't use file of averages at all.
	***/
	bool isRFL = !DistrParams::IsRVL();		// is Read fixed length set
	AvrFrags::AvrFrag avrLocal;				// local temporary averages
	unique_ptr<AvrFrags> avrs;				// averages from file
	if (isRFL)	avrs.reset(new AvrFrags(_cSizes.ServPath()));

	AvrFrags::AvrFrag& avr = isRFL ? avrs->Get() : avrLocal;
	FragLenStat fLenStat;

	if(!avr.ReadLen) {						// new averages file or local averages
		Average		genFrAvr;				// generated frags average length
		bool isMDA = IsMDA;
		IsMDA = false;
		// trial sharing without MDA to find 'clear' average frag length for given distribution
		ULLONG recCnt = CutForSample(genFrAvr, &fLenStat);

		if (isRFL)	avrs->SetMinMax(fLenStat.Min, fLenStat.Max);
		avr.ReadLen = Read::FixedLen;
		avr.SelLen	= fLenStat.SelAvr.Mean();	// selected frags average length
		avr.RecLen	= float(genFrAvr.Sum()) / recCnt;
		IsMDA = isMDA;							// restore MDA
	}
	else {
		fLenStat.Min = avrs->GetMin();
		fLenStat.Max = avrs->GetMax();
	}
	// *** MDA accounting
	if (IsMDA && !avr.MdaLen) {
		// trial sharing to find average frag length for given distribution with MDA
		Average	genFrAvr;
		ULLONG recCnt = CutForSample(genFrAvr, NULL);
		float mdaLen = float(genFrAvr.Sum()) / recCnt;

		if (isRFL)	avrs->SetMDAlen(mdaLen);
		else avr.MdaLen = mdaLen;
	}
	//avrs.Print();
	SelFragAvr = avr.SelLen;

	// *** Set Test samples and Control number of cells and sample
	GlobContext[int(GM::eMode::Test)].SetSample(SAMPLE_FG()/100);
	if(TestMode)
		if( !(GlobContext[int(GM::eMode::Test)].Sample[Gr::BG] *= SAMPLE_BG()/100) )	All = false;
	if(MakeControl)	// can be true only in TEST task mode
		GlobContext[int(GM::eMode::Control)].SetControlSample( GlobContext[int(GM::eMode::Test)].GetExactBGCellCnt() );

	// *** Determine the total possible numbers of recorded reads
	ULLONG	totalCnt = 0;		// total number of recorded reads
	ULONG	maxCnt[] = {0,0};
	float	maxDens[] = {0,0};
	chrlen	enRgnLen;		// length of enriched regions
	// coefficient in formula: fragCnt_onLen = CellCnt * Sample * Len / recordedAvrLen
	// constant countFactor = CellCnt / recordedAvrLen
	const float	countFactor = float(GlobContext[int(GM::eMode::Test)].CellCnt) / 
		(IsMDA ? avr.MdaLen : avr.RecLen);

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
	//if (IsMDA)	totalCnt += totalCnt/5;	// empirical coefficient 1.2: right for small read cnt, but failed for big one
	if(PCRCoeff) {
		ChromCutter::SetAmpl();
		totalCnt *= ULLONG(pow(2.f, int(PCRCoeff)));
	}
	// *** Estimate adjusted Sample
	if(totalCnt > Seq::ReadsLimit())
		AutoSample = Seq::ReadsLimit() / totalCnt;
	// *** print debug info
	if(Verbose(eVerb::PAR)) {
		cout << SignPar << "Actual fragments size" << SepCl
			 << "Mean" << Equel << SelFragAvr << SepSCl
			<< "minimum ~ " << fLenStat.Min << SepSCl
			<< "maximum ~ " << fLenStat.Max << endl;
	}
	if (Verbose(eVerb::DBG)) {
		if (MakeControl)
			cout << SignDbg << "Generated " << GM::Title(GM::eMode::Control) << SepDCl
			<< "Count of cells" << Equel << CellCnt(GM::eMode::Control) << SepSCl
			<< "sample" << Equel << sPercent(Sample(GM::eMode::Control, Gr::BG) * 100, 2) << LF;
		cout << SignDbg << "Total recorded reads number estimate" << SepCl << totalCnt << endl;
	}
	if(AutoSample < 1 && Verbose(eVerb::RES))
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
//				startExtReadPos += Read::FixedLen+1;
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
