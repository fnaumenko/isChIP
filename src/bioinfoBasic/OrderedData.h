/**********************************************************
OrderedData.h (c) 2022 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 11/12/2023
-------------------------
Provides chromosomally sorted data files functionality
***********************************************************/
#pragma once
#include "Data.h"
#include <assert.h>

enum eStrand { TOTAL = 0, POS, NEG, CNT };

using coval = chrlen;				// coverage value
using covmap = map<chrlen, coval>;	// coverage map

static const string	sStrandEXT[] = { strEmpty, ".pos", ".neg" };
static const char*	sStrandTITLES[] = { "positive", "negative" };

static const BYTE	ColorLEN = 11;

// 'AccumCover' represents cumulative chrom's fragment coverage data
//	and implements a single method for gradual filling (incrementing) coverage
class AccumCover : public covmap
{
protected:
	// Calls functor for each point that put the chrom coverage
	template<typename Functor>
	void DoWithItem(Functor f) const { for (const value_type& item : *this)	f(item); }

	// Calls functor for each point that put the chrom coverage
	template<typename Functor>
	void DoWith2Items(Functor f) const {
		auto it0 = cbegin(), it = it0;		// previous, current entry

		for (++it; it != end(); it0 = it++)
			if (it0->second)	f(it0, it);
	}

public:
	// Default constructor
	//AccumCover() = default;

	//// Copy constructor
	//AccumCover(const covmap& cv) : covmap(cv) {}

	//AccumCover(covmap::const_iterator first, covmap::const_iterator last) : covmap(first, last) {}
	//AccumCover(const AccumCover&) = default;
	//AccumCover(AccumCover&) = default;
	//AccumCover(AccumCover&&) = default;
	//~AccumCover() = default;

	// Adds fragment to accumulate the coverage
	void AddRegion0(const Region& frag);
	void AddRegion(const Region& frag);

#ifdef _DEBUG
	void WigPrint() const
	{
		cout << "pos\tval\n";
		DoWithItem([](const auto& item) { cout << item.first << TAB << item.second << LF; });
	}

	// Prints output in BedGraph format
	void BgPrint() const
	{
		cout << "start\tend\tval\n";
		DoWith2Items([](const auto& it0, const auto& it1)
			{ cout << it0->first << TAB << it1->first << TAB << it0->second << LF; }
		);
	}
#endif	// _DEBUG
};

// 'Freq' represents cumulative position frequency
class Freq : public covmap
{
public:
	// Adds item's position
	void AddPos(chrlen pos) { (*this)[pos]++; }

	// Adds fragment's position to accumulate the density
	//	@param frag: added frag
	void AddFragPos(const Region& frag) { AddPos(frag.Centre()); }

	// Adds Read's position to accumulate the density
	//	@param tag: added Read
	//	@param reverse: if true then add complemented read
	void AddReadPos(const Read& tag, bool reverse) { AddPos(reverse ? tag.End() : tag.Start()); }
};

//=====  WRITERS

// BED and WIGGLE track fileds
struct TrackFields
{
	const string Name;
	const char* Descr;
	const string* CommLine = nullptr;
	bool ItemRgb = false;
	bool UseScore = false;
	const char* Color = NULL;

	// Basic constructor
	//	@param name: track name
	//	@param descr: track descriptor
	//	@param commLine: command line saved as a comment
	TrackFields(const string& name, const char* descr, const string* commLine)
		: Name(name), Descr(descr), CommLine(commLine) {}

	// Copy constructor width overridden fields
	TrackFields(const TrackFields& params, bool itemRgb, bool useScore, const char* color = NULL)
		: Name(params.Name), Descr(params.Descr), CommLine(params.CommLine), ItemRgb(itemRgb), UseScore(useScore), Color(color) {}

	// Copy constructor width extended name and overridden fields
	TrackFields(const TrackFields& params, const string& addName, const char* descr, bool itemRgb, bool useScore, const char* color = NULL)
		: Name(params.Name + addName), Descr(descr), CommLine(params.CommLine), ItemRgb(itemRgb), UseScore(useScore), Color(color) {}
};

class RegionWriter : public TxtOutFile
{
protected:
	//static const char* sTRACK;
	static const char* sGRAY;
	static const BYTE sClrLEN = 11;

	RegionWriter(FT::eType ftype, eStrand strand, const TrackFields& fields);

	rowlen AddChromToLine(chrid cID);

	virtual void WriteChromData(chrid cID, const covmap& cover) {};

//public:
	//const string& FileName() const { return TxtFile::FileName(); }
};

// 'WigWriter' is a base class  writing in wiggle formats
class WigWriter : public RegionWriter
{
	static string ChromMarker(chrid cID) { return " chrom=chr" + Chrom::Mark(cID); }

	void WriteFixStepDeclLine(chrid cID, chrlen pos);

public:
	// Creates new instance for writing.
	//	@param ftype: BGRAP or WIG_VAR
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	WigWriter(FT::eType ftype, eStrand strand, const TrackFields& fields)
		: RegionWriter(ftype, strand, fields) {}

	// Fill IO buffer by <position>-<value> lines
	void WriteChromVarStepData(chrid cID, const covmap& cover);

	// Adds to IO buffer declaration line and lines, each containing one value
	//	@param cID: chrom ID
	//	@param pos: range start position
	//	@param vals: range values
	//	@param closure: if true then adds zero value to 'close' the curve for the IGV view
	void WriteFixStepRange(chrid cID, chrlen pos, const vector<float>& vals, bool closure = true);

	// Writes oblique line
	//	@param cID: chrom's ID
	//	@param pos: start (left) position
	//	@param ptCnt: number of points in oblique line
	//	@param shift: shift of value at each point of the oblique line: if < 0 then direct (positive) line, otherwise reversed (negative) one
	void WriteFixStepLine(chrid cID, chrlen pos, chrlen ptCnt, float shift);
};

// 'VarWigWriter' implements methods for writing in wiggle_0 (variable step) format
class VarWigWriter : public WigWriter
{
public:
	// Creates new instance for writing cover to wiggle_0 file.
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	VarWigWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::WIG_VAR, strand, fields) {}

	// Fill IO buffer by chrom cover
	void WriteChromData(chrid cID, const covmap& cover) override { WriteChromVarStepData(cID, cover); }
};

// 'BedGrWriter' implements methods for writing cover in BedGraph format
class BedGrWriter : public WigWriter
{
public:
	// Creates new strand-separated BedGraph instance for writing
	//	@param strand: strand
	//	@param fields: BED/WIG track fields
	BedGrWriter(eStrand strand, const TrackFields& fields)
		: WigWriter(FT::eType::BGRAPH, strand, fields) {}

	// Fill IO buffer by chrom cover
	void WriteChromData(chrid cID, const covmap& cover) override;
};

inline int StrandShift(BYTE dim) { return dim != 2; }

// 'Writers' keeps the set of writers of the same type
template <typename WRITER>
class Writers
{
protected:
public:
	vector<WRITER*> _files;

public:
	// Ñreates writers according to dimension
	//	@param dim: number (dimension) of data; should be 1 (total only), 2 (strands only) or 3 (total and strands)
	//	@param fields: BED/WIG track fields
	Writers(BYTE dim, const TrackFields& fields)
	{
		assert(dim);
		_files.resize(dim, nullptr);
		BYTE shift = StrandShift(dim);
		if (shift) _files[TOTAL] = new WRITER(TOTAL, fields);
		if (dim > 1)
			_files[shift] = new WRITER(POS, TrackFields(fields, sStrandEXT[POS], fields.Descr, false, false)),
			_files[++shift] = new WRITER(NEG, TrackFields(fields, sStrandEXT[NEG], fields.Descr, false, false));
	}

	~Writers() { Do([](WRITER* f) { delete f; }); }

	// Applies function fn to each writer
	void Do(function<void(WRITER*)> fn) { for (auto& f : _files) fn(f); }
	void Do(function<void(WRITER*)> fn) const { for (auto& f : _files) fn(f); }

	// Prints output file names separated by comma
	void PrintNames() const
	{
		const char* sep = strEmpty.c_str();
		Do([&sep](auto f) { cout << sep << f->FileName(); sep = SepCm; });
	}
};

//===== ORDERED DATA

// 'DataSet' keeps the set of data of the same type - total, or strands, or total and strands
template <typename DATA>
class DataSet
{
	vector<DATA> _data;
	BYTE _strandShift;

public:
	bool Closed = true;		// true if data generation is completed or data is empty
	bool Unsaved = true;	// true if data is still unsaved

	// Constructor
	//	@param dim: number (dimension) of data; should be 1 (total only), 2 (strands only) or 3 (total and strands)
	DataSet(BYTE dim = 1) : _strandShift(StrandShift(dim)) { assert(dim); _data.resize(dim); }

	// Returns a direct pointer to the DATA array
	DATA* Data() { return _data.data(); }
	const DATA* Data() const { return _data.data(); }

	// Returnes data from dataset by index, or common data by default
	DATA& DataByInd(BYTE ind = 0) { return _data[ind]; }
	const DATA& DataByInd(BYTE ind = 0) const { return _data[ind]; }

	// Returnes strand data by strand
	DATA& StrandData(eStrand strand) { return _data[strand - !_strandShift]; }
	const DATA& StrandData(eStrand strand) const { return _data[strand - !_strandShift]; }

	// Returnes strand data by index: 0 - POS, 1 - NEG
	DATA& StrandDataByInd(BYTE ind) { return _data[ind + _strandShift]; }
	const DATA& StrandDataByInd(BYTE ind) const { return _data[ind + _strandShift]; }


	void Reinit() { Closed = false; Unsaved = true; }

	void Clear() { for (DATA& d : _data) d.clear(); }

	bool Empty() const {
		for (const DATA& d : _data) if (!d.empty()) return false;
		return true;
	}

	// Returnes true if strands are defined
	bool Strands() const { return _data.size() > 1; }
};

// 'OrderedData' keeps the chromosome datasets and optionally the set of writers these datasets to file.
// Witing to the file is done in order by chromosomes.
template <typename DATA, typename WRITER>
class OrderedData
{
	template <typename DATA>
	struct ChromDataSet : public Chroms<DataSet<DATA>>
	{
		// Initializing constructor
		//	@param cSizes: chrom sizes instance
		//	@param dim: number (dimension) of data of the same type
		ChromDataSet(const ChromSizes& cSizes, BYTE dim) {
			for (const auto& cs : cSizes)
				if (cs.second.Treated)
					Chroms<DataSet<DATA>>::AddVal(cs.first, move(DataSet<DATA>(dim)));
		}
	};

	mutable mutex _mutex;						// only primer mutex is used
	unique_ptr<ChromDataSet<DATA>> _chromsData;	// common chroms data collection (datasets)
	const OrderedData& _primer;					// primer instance to share mutex and collections

protected:
	DataSet<DATA>* _data{};						// current accumulated chromosome data; used 
	unique_ptr <Writers<WRITER>> _writers;		// writers set

	// Primer constructor without writers
	//	@param cSizes: chrom sizes
	//	@param dim: chrom's dataset dimension
	OrderedData(const ChromSizes& cSizes, BYTE dim) : _primer(*this)
	{
		assert(dim);
		_chromsData.reset(new ChromDataSet<DATA>(cSizes, dim));
	}

private:
	// Writes chrom's data to writer
	//	@param cID: chrom ID
	//	@param data: chrom's data
	//	@param clearData: if true then the chrom's data will be clear
	void WriteChromData(const chrid cID, DataSet<DATA>& data, bool clearData)
	{
		if (data.Empty())	return;		// no chrom data in input file
		data.Unsaved = false;
		BYTE i = 0;
		_primer._writers->Do([&](auto f) { f->WriteChromData(cID, data.DataByInd(i++)); });
		if (clearData)
			data.Clear();
	}

	// Primer constructor
//	@param cSizes: chrom sizes
//	@param dim: number (dimension) of data
//	@param write: if true then data sould be save to output file
//	@param fields: BED/WIG track fields
	//OrderedData(const ChromSizes& cSizes, BYTE dim, bool write, const TrackFields& fields, const string* commLine = nullptr)
	//	: OrderedData(cSizes, dim)
	//{
	//	if (write)	_writers.reset(new Writers<WRITER>(dim, fields, commLine));
	//}

public:

	// Primer constructor
	//	@param cSizes: chrom sizes
	//	@param dim: number (dimension) of data
	//	@param write: if true then data sould be save to output file
	//	@param fname: name of output file
	//	@param descr: track decsription in declaration line
	//	@param commLine: command line saved as a comment
	OrderedData(const ChromSizes& cSizes, BYTE dim, bool write, const string& fname, const char* descr, const string* commLine = nullptr)
		: OrderedData(cSizes, dim)
	{
		if (write)	_writers.reset(new Writers<WRITER>(dim, TrackFields(fname, descr, commLine)));
	}

	// Clone constructor for multithreading
	//	@param data: primer data
	OrderedData(const OrderedData& data) : _primer(data) {}

	// Returns chrom's data
	DataSet<DATA>& ChromData(chrid cID) { return _chromsData->Data(cID); }
	const DataSet<DATA>& ChromData(chrid cID) const { return _chromsData->Data(cID); }

	void Clear() { _data->Clear(); }

	// Sets chrom as current (for accumulating only)
	void SetChrom(chrid cID) { (_data = &(_primer._chromsData->Data(cID)))->Reinit(); }

	// Save chrom's data by defined writers in chromosome order; of no writes are set, does nothing
	//	@param cID: chrom ID
	//	@param clearData: true if chrom's data should be cleaned after all
	void WriteChrom(const chrid cID, bool clearData = true)
	{
		// Chroms data are filled in different threads independently.
		// To save them in chrom sorted order the total pool is examined each time the next data is completed.
		// If previous data are filled without gups in chrom order, they are recorded and optionally removed from the pool.

		if (!_primer._writers)	return;
		if (Mutex::isOn())	
			_primer._mutex.lock();
		DataSet<DATA>& data = _primer._chromsData->Data(cID);
		bool save = data.Closed = true;					// mark the recorded chromosome data as closed
		for (auto it = _primer._chromsData->cBegin(); CID(it) < cID; it++)
			if (it->second.Data.Unsaved)
				if (bool(save = it->second.Data.Closed))	// chrom data is closed, save it
					WriteChromData(CID(it), _primer._chromsData->Data(CID(it)), clearData);
				else
					break;								// some of chrom data aren't saved
		if (save)
			WriteChromData(cID, data, clearData);		// write the recorded chromosome data
		if (Mutex::isOn())	_primer._mutex.unlock();
	}

	// For current chrom adds SE fragment to total coverage, and to strands coverage if strands are defined
	//	@param frag: added fragment
	//	@param reverse: true if read is reversed (neg strand)
	void AddFrag(const Region& frag, bool reverse) {
		_data->DataByInd().AddRegion(frag);
		if (_data->Strands())
			_data->StrandDataByInd(reverse).AddRegion(frag);
	}

	// For current chrom adds fragment to total density
	//	@param frag: added fragment
	void AddFragDens(const Region& frag) { _data->DataByInd().AddFragPos(frag); }

	// For current chrom adds SE read to total coverage
	//	@param tag: added read
	//	@param reverse: true if read is reversed (neg strand)
	void AddReadDens(const Read& tag, bool reverse) { _data->DataByInd().AddReadPos(tag, reverse); }

	// Prints output file names separated by comma
	void PrintWritersName() const { if (_primer._writers)	_primer._writers->PrintNames(); }
};