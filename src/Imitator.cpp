#include "isChIP.h"
#include "Imitator.h"
#include <algorithm>    // std::sort
#include <cfloat>		// FLT_MAX
#include <cwchar>		// long '-'

static const char* cAverage = "Average";
#define PI 3.141593

/************************  class Average ************************/
UINT Average::operator+=(UINT val)
{
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

/************************  FragDistr ************************/

//Imitator::ChromCutter::FragDistr::tpLognormNext Imitator::ChromCutter::FragDistr::pLognormNext;

//float Imitator::ChromCutter::FragDistr::ssDVar;	// Doubled variance 2*sigma^2 in the size sel norm distr
//float Imitator::ChromCutter::FragDistr::ssFactor;	// Factor sigma*sqr(2PI) in the size sel norm distr
//float Imitator::ChromCutter::FragDistr::ssAlRatio;// Size sel normal distribution aligning up ratio
float Imitator::ChromCutter::FragDistr::ssFactor0;	// factor sigma*sqrt(2) in the size sel norm distr
float Imitator::ChromCutter::FragDistr::ssFactor1;	// factor 2.5/sqrt(2PI) in the size sel norm distr

/*
First release (4 vars defined lognormal distribution):
X_lognorm = exp( (Normal()*Sigma + Mean) / LnFactor + LnTerm );
where LnFactor = 500, LnTerm = 5.1,
Normal(): mean = 200, SD = 200

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
	double normal_x1;		// first random coordinate
	double w;				// radius
	
	if (_phase) {			// we have a valid result from last call
		_phase = 0;
		return _normal_x2;
	}    
	do {					// make two normally distributed variates by Box-Muller transformation
		normal_x1 = 2. * DRand() - 1.;
		_normal_x2 = 2. * DRand() - 1.;
		w = normal_x1 * normal_x1 + _normal_x2 * _normal_x2;
	} while (w >= 1. || w < 1E-30);

	w = sqrt( log(w) * (-2./w) );
	_phase = 1;
	_normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
	
	return normal_x1 * w;		// return normal distributed value
}

// Initializes size sel
void Imitator::ChromCutter::FragDistr::Init() 
{
	if(DistrParams::IsSS()) {
		//ssDVar = 2 * ssSigma * ssSigma;
		//ssFactor = ssSigma * sqrt(2 * PI);
		//ssAlRatio = 2.5 * ssSigma;	// // *2.5 because normilized (mean=0, sigma=1) Gauss() = 0.4
		ssFactor0 = (float)(DistrParams::ssSigma * sqrt(2.f));
		ssFactor1 = 2.5f / (float)sqrt(2 * PI);
	}
}

/*
float FragDistr::LognormalFactorCalc()
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

float FragDistr::LognormalMeanCalc()
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

/************************  FragDistr: end ************************/

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

/************************ AvrFragLengths ************************/

#define	SZ_FCTR	0	// index of file size factor
#define	SEL_MIN		2	// index of selected min frag length
#define	SEL_MAX		3	// index of selected max frag length
#define	SEL_AVR		4	// index of selected average frag length
#define	GEN_AVR		5	// index of generated average frag length
#define	REC_AVR		6	// index of recorded average frag length for given MDA, PCR ampl coeff

// 'AvrFragLengths' reads and writes average legths of fragment into plain text file.
class AvrFragLengths
/*
 * Class 'AvrFragLengths' reads and writes average legths of fragment into plain text file.
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
	//static const BYTE	AmplInd = 6;		// Index of amplification coefficient in param (averages array)

	// Creates a new instance of AvrFragLengths from given path;
	// a file name will be constructed.
	// If file does not exist, the instance is empty.
	AvrFragLengths	(const string path) : _isChanged(false)
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
		_avrLens.reserve(REC_AVR/2 + 2);	// 2 items for coefficients
		while( file.GetLine() )
			_avrLens.push_back(AvrFrag(a_cycle(file.IntField(0)), 
				file.FloatField(1), file.FloatField(2)));
	}
	
	// Writes instance to a file if it's changed.
	~AvrFragLengths()
	{
		if(!_isChanged)		return;
		const char* title = "# sampling generated by isChIP; do not change";
		const char* comms[] = { 
			"[1]:unused; [2]:file size factor; [3]:unused",
			"[1]:unused; selected frag len: [2]:min, [3]:max",
			"[1]:unused; frag avrg len: [2]:selected, [3]:generated",
			"[1]:ampl coeff; recorded frag avrg len: [2]:MDA, [3]:PCR"
		};
		const BYTE commsCnt = sizeof(comms)/sizeof(char*);
		TxtOutFile file(_fileName);

		BYTE i = 0, capacity = 5;	// size factor float capacity = 5, others = 6
		file.SetLineBuff(2*(2+CHRLEN_CAPAC) + strlen(comms[3]));	// comms[3] is the longest string
		file.WriteLine(title);
		for(vector<AvrFrag>::iterator it = _avrLens.begin(); it != _avrLens.end(); it++) {
			file.WriteLine(it->aCoeff, it->fLen, it->sLen, capacity, i<commsCnt ? comms[i++] : "");
			capacity = 6;
		}
	}
	
	// Gets size Factor and fragment length average values.
	//	@isPCR: true if PCR is established
	//	@coeff: user-defined coefficient of amplification
	//	@avrs: extern array to return values; first REC_AVR values are always returned,
	//	last 2 only for given coefficient
	//	return: true if values for given coefficient of amplification are returned successfully
	bool Get(bool isPCR, BYTE coeff, float avrs[]) const
	{
		short i = 0;
		for(vector<AvrFrag>::const_iterator it=_avrLens.begin(); it!=_avrLens.end(); it++)
			if(i < REC_AVR) {
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
			for(vector<AvrFrag>::iterator it = _avrLens.begin() + REC_AVR/2; it != _avrLens.end(); it++)
				if( _isChanged = it->aCoeff == coeff ) {
					it->fLen = avrs[REC_AVR];
					it->sLen = avrs[REC_AVR+1];
					break;
				}
			if( !_isChanged ) {	// add new coefficient
				_avrLens.push_back(AvrFrag(coeff, avrs[REC_AVR], avrs[REC_AVR+1]));
				sort(_avrLens.begin() + REC_AVR/2, _avrLens.end());
			}
			// do average with new frag average length
			_avrLens[SEL_AVR/2].fLen = (_avrLens[SEL_AVR/2].fLen + avrs[SEL_AVR]) / 2;
		}
		else {					// new instance
			short i = SZ_FCTR;
			for(; i<REC_AVR; i+=2)
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

/************************ AvrFragLengths: end ************************/

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

const char* Imitator::ChromView::sMargTime = "    ";	// margin string before time output

// Prints str on given field and return field width
//	@width: field width
int PrFittedStr(const char* str, int width)
{
	cout << setw(width) << str;
	return width;
}

// Prints str on given field and return field width
//	@preWidth: width of field before str
//	@width: str field width
void PrFittedStr(const char* str, int preWidth, int width)
{
	printf("%*c%*s", preWidth, BLANK, width, str);

	//cout << setw(preWidth) << BLANK;
	//cout << setw(width) << str;
	//return preWidth + width;
}

// Prints percent value on given field or nothing if percent is 0
//	@val: value
//	@preWidth: width of field before value
//	@width: value field width
void PrFittedPerc(float val, int precision, int preWidth, int width)
{
	cout << setw(preWidth) << BLANK;
	cout << setw(width);
	if(val)	cout << sPercent(val, precision, width);
	else	cout << BLANK;
}

// Prints unzero float val or nothing on given field
//	@width: field width
void PrFittedFloat(float val, int preWidth, int width, int precision)
{
	printf("%*c", preWidth, BLANK);
	if(val)	printf("%*.*f", width, precision, val);
	else	printf("%*c", width, BLANK);

	//cout << setw(preWidth) << BLANK;
	//cout << setw(width);
	//if(val)	cout << setprecision(precision) << val;
	//else	cout << BLANK;
}

// Prints unzero float val or nothing on given field
//	@width: field width
void PrFittedInt(ULONG val, int preWidth, int width)
{
	cout << setw(preWidth) << BLANK;
	cout << right << setw(width);
	if(val)	cout << val;
	else	cout << BLANK;
}

// Prints chrom's header: Reads
//	@FgHeader: true for foreground header
//	return: header width
int Imitator::ChromView::PrintHeader(bool FgHeader)
{
	int lineW = 0;
	// shift title to right in case of little reads count
	bool w = CountW+1 < FT::ItemTitle(FT::ABED, true).length();	

	cout << right;
	// if little reads count then extend "reads" field to 1, otherwhise extend "sample" field to 1
	lineW += PrFittedStr(FT::ItemTitle(FT::ABED, true).c_str(), 
		GrTitleW() + CountW + marg_R*ControlMode + w);		// "reads"
	lineW += PrFittedStr("sample", margR_S + SampleW + !w);						// "sample"
	// now we should compensate extention to 1:
	// if MDA is set then shrink "MDAc" to 1, otherwhise shrink "unit densuty" field to 1
	if(w = ChromCutter::IsMDA())
		lineW += PrFittedStr("MDAc", margS_ + AmplW - 1);	// MDA coefficient 
	lineW += PrFittedStr(UnitDens, marg_D + DensW - !w);	// unit density
	if(FgHeader)	PrintMarg(margD_), lineW += margD_;		// empty space
	
	return lineW;
}

// Prints chrom's header: 'N' statistics
int Imitator::ChromView::PrintHeaderN()
{
	int lineW = 0, w = 3;	// 3 for right shift 

	if(Nts::StatN) {
		lineW = PrFittedStr("N", margD_ + NW - 1);
		if(!Nts::LetN )
			lineW += PrFittedStr("N_excl", w + margN_Nex + NW);	// reduce shift
	}
	if(Timer::Enabled)
		lineW += PrFittedStr("mm:ss", (w-=2) + strlen(sMargTime) + 2);
	return lineW;
}

// Print info about recorded Reads (count, sample, ampl, density) for ground type defined in this nstance
//  @fragCnt: frags counters
//	@gMode: gen mode for which info is printed
//	@densLen: length to print density
void Imitator::ChromView::PrintReads(const FragCnt fragCnt, GM::Mode gMode, chrlen densLen)
{
	ULONG rCnt = ULONG(fragCnt.RecCnt() << Seq::Mode());	// count of reads

	// == Gr title
	if(TestMode)
		if(gMode == GM::Control)		// can be true only if Imitator::MakeControl is true
			PrintMarg(Gr::TitleLength + 1);				// print blanks instead of Ground title
		else
			cout << right << Gr::Title(GrType) << COLON;// print of Ground title
	// == count of Reads
	PrFittedInt(rCnt, marg_R, CountW);
	// == sample
	if(MakeControl && gMode==GM::Control)
		PrintMarg(margR_S + SampleW);	// don't print sample for accompanying input
	else
		PrFittedPerc(fragCnt.Sample(), SamplePr, margR_S, SampleW);
	// == print MDA coef
	if(ChromCutter::IsMDA())
		PrFittedFloat(fragCnt.RealAmplCoeff(), margS_, AmplW, AmplPr);
	// == print density
	PrFittedFloat(ReadDens(rCnt, densLen), marg_D, DensW, DensPr);
	//float val = ReadDens(rCnt, densLen);
	//PrFittedFloat(val, marg_D, DensW, /*val < DensMax/2 ? DensPr-1 :*/ DensPr);//, val && val < 10);

	if(GrType == Gr::FG)		PrintMarg(margD_);
}

// Prints chroms 'N' statistics
void Imitator::ChromView::PrintN(const Nts& nts)
{
	if(!Nts::StatN)		return;
	float p = Percent(nts.CountN(), nts.Length());
	PrFittedStr(p? sPercent(p, 2, NW).c_str() : sBLANK, margD_, NW);
	if(Nts::LetN)		return;
	PrFittedStr(p ? sPercent(nts.UndefLengthInPerc(), 2, NW).c_str() : sBLANK, margN_Nex, NW);
}

// Initializes viewing data
void Imitator::ChromView::Init(ULONG maxCnt, float sample, float maxDens)
{
	 // +1 not really good; need to estimate total by ground
	CountW = DigitsCount(maxCnt, Options::GetBVal(oLOCALE)) + 1;	// +1 for total
	//if(!CountW)		CountW = 1;		// if BG level is 0
	sample *= 100;
	if(sample <= 0.01)		SampleW = 7;		// 0.0012%
	else if(sample <= 0.1)	SampleW = 6;		// 0.012%
	else if(sample <= 1)	SampleW = 5;		// 0.12%
	else if(sample < 10)	SampleW = 4;		// 1.2%
	else if(sample < 100)	SampleW = 3;		// 12%
	else					SampleW = 4;		// 100%
	DensMax = maxDens;
	if(maxDens < 0.015)		DensW = 6, DensPr = 3;	// 0.0012
	else if(maxDens < 0.15)	DensW = 5, DensPr = 3;	// 0.012
	else if(maxDens < 1.5)	DensW = 5, DensPr = 2;	// 0.123
	else if(maxDens < 15)	DensW = 4, DensPr = 1;	// 1.23
	else if(maxDens < 150)	DensW = 5, DensPr = 0;	// 123.4
	else					DensW = 6, DensPr = 0;	// 1234.5
	// this variant for cout (std::setprecision() works different with printf(%*.*f",width,precision,val)
	//DensPr = 2;
	//if(maxDens < 0.015)		DensW = 6;				// 0.0012
	//else if(maxDens < 0.15)	DensW = 5;				// 0.012
	//else if(maxDens < 1.5)	DensW = 5, DensPr = 3;	// 0.123
	//else if(maxDens < 15)	DensW = 4, DensPr = 3;	// 1.23
	//else if(maxDens < 150)	DensW = 5, DensPr = 4;	// 123.4
	//else					DensW = 6, DensPr = 5;	// 1234.5
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
	cout << "DensMax = " << DensMax << EOL;
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
ULONG Imitator::ChromCutter::SetGMode(GM::Mode gmode)
{
	_gMode = gmode;
	_fragCnt.SetGMode(gmode);
	_output->SetGMode(gmode);
	return CellCnt(_gMode);
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
	_cFiles(imitator->_cFiles),
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
//  @nts: current chromosome
//	@enrRegLen: length of all enriched regions
//	@timer: current timer to thread-saves time output or NULL
//	@exceedLimit: true if limit is exceeded
void Imitator::ChromCutter::PrintChrom (
	chrid cID, const Nts& nts, chrlen enrRegLen, Timer& timer, bool excLimit)
{
	if( !Verbose(vRT) )	return;
	const ULONG regLens[] = { enrRegLen, nts.Length() - enrRegLen }; // FG, BG region's lengths

	Mutex::Lock(Mutex::OUTPUT);

	PrintChromInfo(cID, _gMode, _fragCnt.GetFragCnts(), regLens, !IsSingleThread());
	ChromView::PrintN(nts);							// print 'N' statistics
	timer.Stop(ChromView::sMargTime, false, false);	// print time
	if(excLimit)	cout << " exceeded limit of " << Seq::ReadsLimit() << " reads";
	cout << endl;
	
	Mutex::Unlock(Mutex::OUTPUT);
}

// Treats chromosomes given for current thread
//	@cSubset: pointer to Subset - set of chrom IDs treated in this thread
thrRetValType Imitator::ChromCutter::Execute(const effPartition::Subset& cIDSet)
{
	chrid	cID;
	BedF::cIter	cit;			// template chrom's iterator
	ULONG	n, cellCnt;			// count of cells, length of enriched regions
	chrlen	currPos, k, fCnt;	// count of features
	chrlen*	pCurrPos = &currPos;// pointer to current position
	chrlen	enrRegLen;			// length of enriched regions
	short	lineW = 0;
	Timer	timer(Verbose(vRT));// print local time on Verbose 'runtime info and above'

	try {
		for(effPartition::numb_id_cit it=cIDSet.NumbIDs().begin(); it!=cIDSet.NumbIDs().end(); it++) {
			timer.Start();
			cID = *it;
			_fragCnt.Clear();
			// multiply twice for digits (diploid) chrom
			cellCnt = SetGMode(GM::Test) << _cFiles[cID].Numeric();
			// print chrom name before cutting
			PrintChromName(cID, _gMode, IsSingleThread());
			_output->BeginWriteChrom(cID);

			if(Bed && (cit=Bed->GetIter(cID)) != Bed->cEnd()) {
				fCnt = Bed->FeatureCount(cit);
				enrRegLen = Bed->EnrRegLength(cit, 0, SelFragAvr);
			}
			else	enrRegLen = fCnt = 0;
			Nts nts(_cFiles.FileName(cID));

			for(n = 0; n < cellCnt; n++) {
				currPos = nts.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
				for(k=0; k < fCnt; k++)
					if(lineW = CutChrom(nts, pCurrPos, Bed->Feature(cit, k), false))
						goto A;			// achievement of Reads limit
				// add background after last 'end' position
				if((lineW = CutChrom(nts, pCurrPos, nts.DefRegion(), true)) < 0 && All) {
					_output->EndWriteChrom(cID);
					break;				// achievement of Reads limit
				}
			}
A:			PrintChrom(cID, nts, enrRegLen, timer, lineW < 0);		// timer stops and printed in here
			IncrTotalSelFragCount();
			// correct total enriched regions length to calculate total density
			IncrementTotalLength(nts.DefRegion().Length(), enrRegLen);
			if(lineW < 0) {
				_output->EndWriteChrom(cID);
				break;			// achievement of Reads limit
			}
			if(MakeControl) {
				timer.Start();
				cellCnt = SetGMode(GM::Control) << _cFiles[cID].Numeric();
				PrintChromName(cID, _gMode, IsSingleThread());	// print chrom name before cutting
				for(n = 0; n < cellCnt; n++) {
					currPos = nts.Start() + _fragDistr.RandFragLen();	// random shift from the beginning
					CutChrom(nts, pCurrPos, nts.DefRegion(), true);
				}
				PrintChrom(cID, nts, enrRegLen, timer, false);		// timer stops in here
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
//	@nts: cutted chromosome
//	@currPos: cutting start position
//	@feature: current treated feature
//	@bg: if true then generate background (swap foreground and background)
//	@fragStat: statistics of selected fragments average counter, or NULL under work mode
//	return: 0 if success,
//		1 if end chromosome is reached (continue treatment),
//		-1 if Reads limit is achieved (cancel treatment)
int Imitator::ChromCutter::CutChrom	(
	const Nts& nts,
	chrlen* const currPos,
	const Featr& feature,
	bool bg,
	FragLenStat* fragStat
	)
{
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
	const chrlen cLen = nts.Length();	// chrom end position
	a_cycle	i;				// PCR amplification counter
	int	res;				// result of AddRead() method

	readscr	featScore[] = { 	// current feature's score:
		UniformScore ? 1 : feature.Score,	// FG
		1									// BG
	};

	for(; *currPos<=feature.End; *currPos += fragLen) {	// ChIP: control right mark
		fragLen = _fragDistr.LognormNext();
		if(*currPos + fragLen > cLen)	return 1;		// end of chrom
		_fragDistr.SizeSelLimNext(fragLenMin, fragLenMax);
		/*
		 * Since the lower limit of the fragment length after size selection remains unchanged
		 * after amplification, in order to increase the efficiency, it is cut off immediately.
		 * The upper limit is checked after amplification, since the fragment length can decrease with MDA.
		 */
		if(fragLen < fragLenMin)	
			continue;		// size selection 1: skip short fragment
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
							if(!_output->AddRead(nts, *currPos+fracShift, fragLen, g)
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

Imitator::Context	Imitator::GlobContext[2];	// global generation context
Imitator::ChromView Imitator::ChrView[] = {Gr::FG,Gr::BG};	// 
ULONG	Imitator::TotalLen[] = {0,0};	// FG, BF total length; needed to define total density
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
		if(MakeControl) 
			cout << *GM::Title(gm) << BLANK;
		cout << setw(ChromView::ChromNameW()) << left << setfill(BLANK)
				<< (cID==CHRID_UNDEF ? Total + COLON : Chrom::AbbrName(cID, true) + COLON);
		fflush(stdout);		// including reset to default right and setfill
	}
}

// Prints chromosome's name and full statistics
//	@cID: chromosomes ID, or CHRID_UNDEF to print "total" instead chrom name
//	@fCnts: array of fragment's counters, for FG and BG
//	@regLens: array of region's lengths to print FG|BG density
//	@printChrName: true if chromosome's name should be printed
void Imitator::PrintChromInfo(
	chrid cID, GM::Mode gMode, const FragCnt* fCnts, const ULONG regLens[], bool printChrName)
{
	PrintChromName(cID, gMode, printChrName);
	if(TestMode)	
		ChrView[Gr::FG].PrintReads(fCnts[Gr::FG], gMode, regLens[Gr::FG]);
	ChrView[Gr::BG].PrintReads(fCnts[Gr::BG], gMode, regLens[Gr::BG]);
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
		w += ChromView::PrintHeaderN();
		cout << EOL;
	}
	PrintHorLine(w);
}

// Increments grounds total length.
void Imitator::IncrementTotalLength(chrlen defRegLen, chrlen enrRegLen)
{
	if(enrRegLen)	InterlockedExchangeAdd(&(TotalLen[Gr::FG]), enrRegLen);
	InterlockedExchangeAdd(&(TotalLen[Gr::BG]), defRegLen - enrRegLen);
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

// Runs task in current mode and write result to output files
void Imitator::Execute(BedF* templ)
{
	Bed = templ;
	CutGenome();

	// print statistics
	if( Verbose(vRES) ) {
		if(_cFiles.ChromCount() > 1)	// print summary test statistics?
			PrintChromInfo(CHRID_UNDEF, GM::Test, GlobContext[GM::Test].fCnts, TotalLen), cout << EOL;
		if(TestMode) {
			cout << Total << " recorded "; OutFiles::PrintItemTitle();
			if(MakeControl)
				cout << SepCl << GM::Title(GM::Test) << COLON;
			cout << BLANK << GlobContext[GM::Test].RecCnt();
			if(MakeControl)	
				cout << SepCm << GM::Title(GM::Control) << SepCl << GlobContext[GM::Control].RecCnt();
			cout << endl;
		}
	}
}

// Curs genome into fragments and generate output
void Imitator::CutGenome	()
{
	Array<Thread*> slaves(ThrCnt-1);
	effPartition cSets(_cFiles, thrid(slaves.Length()+1));
	
	if( slaves.Length() && Verbose(vDBG))	cSets.Print();
	
	if(FlatLen < 0)		FlatLen = -FlatLen;
	SetSample();
	PrintHeader(true);
	for(BYTE i=0; i<slaves.Length(); i++)	// invoke slave threads
		slaves[i] = new Thread(StatCutChrom, &cSets[i+1]);
	CutChrom(cSets[0], false);				// invoke main thread
	// wait for slave threads finishing
	for(BYTE i=0; i<slaves.Length(); i++) {
		slaves[i]->WaitFor();
		delete slaves[i];
	}
	PrintHeader(false);
}

// Sets adjusted Sample and clear all counter and means.
// Sample are needed for the control of BF&FG levels by percent (given by user),
// and to prorate number of written reads for each chromosome depending on reads limit.
void Imitator::SetSample()
{
	float avr[REC_AVR+2];
	FragLenStat frLenStat;
	AvrFragLengths avrgs(_cFiles.Path());
	const chrlen maxcLen = 400000000;	// conditional maximum chrom length; ~ 2*mm9:chr1

	memset(avr, 0, sizeof(avr));
	// *** Get averages from file if it exists, otherwise calculate and save ones
	if(!avrgs.Get(ChromCutter::IsPCR(), ChromCutter::AmplCoeff(), avr)) {
		chrlen pos = 0;
		Average genFrAvr;	// generated frags average length
		const ChromSizes cSizes(_cFiles, Imitator::Verbose(vRT));
		ChromCutter cCutter(this, &genFrAvr, false);
		GlobContext[GM::Test].SetSample(1.0);
		
		// generate statistics based on formal chrom length; Reads are not recorded
		cCutter.CutChrom(Nts(), &pos, Region(0, maxcLen), false,
			avr[0] ? NULL : &frLenStat);	// if szFactor = 0, file doesn't exist
		GlobContext[GM::Test].ClearFragCounters();
		
		if(!avr[SZ_FCTR]) {		// file doesn't exist, get szFactor and SelFragAvr at first time
			avr[SZ_FCTR] = float(cSizes[_cFiles.FirstTreatedChromID()]) / 
				_cFiles.FirstTreatedFileLength();
			avr[SEL_MIN] = frLenStat.Min;
			avr[SEL_MAX] = frLenStat.Max;
			avr[SEL_AVR] = frLenStat.SelAvr.Mean();
		}
		avr[GEN_AVR] = genFrAvr.Mean();				// generated frags average length
		avr[REC_AVR + ChromCutter::IsPCR()] = float(genFrAvr.Sum()) / cCutter.RecFgFragCnt();
		avrgs.Add(ChromCutter::AmplCoeff(), avr);	// avrgs is saved in destructor

		//if( Verbose(vDBG) ) {
		//	cout << SignDbg << "Generate by " << Chrom::ShortName(_cFiles.FirstChromID()) << ":\n";
		//	cout << SignDbg << "frags" << SepCl << FragDistr::CallsCnt()
		//		 //<< "  reads" << SepCl << cCutter._recFragCnt[0] << EOL;
		//		 << "  reads" << SepCl << cCutter._fragCnt[0].Rec << EOL;
		//	cout << SignDbg << "coeff of ampl" << SepCl << MDA::SimpleMean() << EOL;
		//}
	}
	SelFragAvr = avr[SEL_AVR];
	//avrgs.Print();
	
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
	chrlen	enRegLen;		// length of enriched regions
	// coefficient in formula: fragCnt_onLen = CellCnt * Sample * Len / recordedAvrLen
	// constant countFactor = CellCnt / recordedAvrLen
	float	countFactor = float(GlobContext[GM::Test].CellCnt) / avr[REC_AVR+ChromCutter::IsPCR()];

	// *** estimate total number of frags
	for(ChromFiles::cIter it=_cFiles.cBegin(); it!=_cFiles.cEnd(); it++) {
		if( !_cFiles.IsTreated(it) )	continue;
		// count is estimated according to diploid (numerical) sign,
		//	but density not, because basic length is single!
		if( Bed && Bed->FindChrom(CID(it)) ) {
			enRegLen = Bed->EnrRegLength(CID(it), 0, SelFragAvr);
			totalCnt += GetReadsCnt(Gr::FG, enRegLen, countFactor, 0, maxCnt, maxDens);
		}
		else	enRegLen = 0;
		totalCnt += GetReadsCnt(Gr::BG,	
			_cFiles.TreatedLength(it, avr[SZ_FCTR]) - enRegLen,
			countFactor, it->second.Numeric(), maxCnt, maxDens);
	}
	// *** Estimate adjusted Sample
	if( totalCnt > Seq::ReadsLimit() ) {
		AutoSample = Seq::ReadsLimit() / totalCnt;
		if( Verbose(vRES) )
			cout << "Added recovery sample = " << setprecision(2) << (AutoSample * 100)
			<< "% due to reads limit of " << Seq::ReadsLimit() << setprecision(0) << endl;
	}
	// *** print debug info
	if(Verbose(vPAR)) {
		cout << SignPar << "Selected fragments size" << SepCl
			 << sActual << "mean" << Equel << SelFragAvr
			 << SepSCl << "minimum ~ " << avr[SEL_MIN]
			 << SepSCl << "maximum ~ " << avr[SEL_MAX]
			 << EOL;
		if(MakeControl)
			cout << SignPar << "Generated " << GM::Title(GM::Control) << SepDCl	
				 << "Count of cells" << Equel << CellCnt(GM::Control) << SepSCl
				 << "sample" << Equel << sPercent(Sample(GM::Control, Gr::BG)*100, 2) << EOL;
	}
	if(Verbose(vDBG))
		cout << SignDbg << "Total recorded reads initial estimate" << SepCl << totalCnt << EOL;
	// *** set Reads statistics params
	if(TestMode)	
		InitReadsView(Gr::FG, maxCnt, maxDens);
	InitReadsView(Gr::BG, maxCnt, maxDens);
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
