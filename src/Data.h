/**********************************************************
Data.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides common data functionality
***********************************************************/
#pragma once

#ifndef _DATA_H
#define _DATA_H

#include "DataInFile.h"
#include <array>
#include <algorithm>    // std::sort

#define	CID(it)	(it)->first

//static const string range_out_msg = "ChromMap[]: invalid key ";

// 'ChromMap'
template <typename T> class ChromMap
{
//protected:
public:
	typedef map<chrid, T> chrMap;

private:
	chrMap _cMap;

protected:
	// Returns count of elements.
	inline size_t Count() const { return _cMap.size(); }

	// Adds class type to the collection without checking cID.
	// Avoids unnecessery copy constructor call
	//	return: class type collection reference
	inline T& AddElem(chrid cID, const T& val) { return _cMap[cID] = val; }

	// Adds empty class type to the collection without checking cID
	//	return: class type collection reference
	inline T& AddEmptyElem(chrid cID) { return AddElem(cID, T()); }

	inline const chrMap& Container() const { return _cMap; }

	inline chrMap& Container() { return _cMap; }

public:
	typedef typename chrMap::iterator Iter;			// iterator
	typedef typename chrMap::const_iterator cIter;	// constant iterator

	// Returns a random-access constant iterator to the first element in the container
	inline cIter cBegin() const { return _cMap.begin(); }
	// Returns the past-the-end constant iterator.
	inline cIter cEnd()	  const { return _cMap.end(); }
	// Returns a random-access iterator to the first element in the container
	inline Iter Begin()	{ return _cMap.begin(); }
	// Returns the past-the-end iterator.
	inline Iter End()	{ return _cMap.end(); }

	// Copies entry
	inline void Assign(const ChromMap& map) {	_cMap = map._cMap; }
	
	// Returns constant reference to the item at its cID
	inline const T& At(chrid cID) const { return _cMap.at(cID); }

	// Returns reference to the item at its cID
	inline T& At(chrid cID) { return _cMap.at(cID); }

	inline const T& operator[] (chrid cID) const { return _cMap(cID); }

	inline T& operator[] (chrid cID) { return _cMap[cID]; }

	// Searches the container for a cID and returns an iterator to the element if found,
	// otherwise it returns an iterator to end (the element past the end of the container)
	inline Iter GetIter(chrid cID) { return _cMap.find(cID); }

	// Searches the container for a cID and returns a constant iterator to the element if found,
	// otherwise it returns an iterator to cEnd (the element past the end of the container)
	inline const cIter GetIter(chrid cID) const { return _cMap.find(cID);	}

	// Adds value type to the collection without checking cID
	inline void AddVal(chrid cID, const T & val) { _cMap[cID] = val; }

	// Removes from the container an element cID
	inline void Erase(chrid cID) { _cMap.erase(cID); }

	// Clear content
	inline void Clear() { _cMap.clear(); }

	// Returns true if element with cID exists in the container, and false otherwise.
	inline bool FindItem (chrid cID) const { return _cMap.count(cID) > 0;	}
};

// 'ChromData' implements sign whether chrom is involved in processing, and chrom's data itself
template <typename T> struct ChromData
{
	bool Treated = true;	// true if chrom is involved in processing
	T	 Data;				// chrom's data

	inline ChromData() : Data(T()) {}
	inline ChromData(const T& data) : Data(data) {}
};

// Basic class for chromosomes collection; keyword 'abstract' doesn't compiled in gcc
template <typename T> 
class Chroms : public ChromMap<ChromData<T> >
{
public:
	// Returns true if chromosome by iterator should be treated
	//inline bool IsTreated(typename ChromMap<ChromData<T> >& data) const { return data.second.Treated; }
	
	// Returns true if chromosome by iterator should be treated
	inline bool IsTreated(typename ChromMap<ChromData<T> >::cIter it) const { return it->second.Treated; }

	// Returns true if chromosome by ID should be treated
	//inline bool IsTreated(chrid cID) const { return IsTreated(this->GetIter(cID));	}

	inline const T& Data(typename ChromMap<ChromData<T> >::cIter it) const { return it->second.Data; }
	
	inline T& Data(typename ChromMap<ChromData<T> >::Iter it) { return it->second.Data; }

	// Returns count of chromosomes.
	inline chrid ChromCount() const { return chrid(this->Count()); }	// 'this' required in gcc

	// Gets count of treated chroms
	//chrid TreatedCount() const 
	//{
	//	chrid res = 0;
	//
	//	for(auto it=this->cBegin(); it!=this->cEnd(); res += IsTreated(it++));
	//	return res;
	//}

	// Returns true if chromosome cID exists in the container, and false otherwise.
	inline bool FindChrom(chrid cID) const { return this->FindItem(cID); }

	// Adds value type to the collection without checking cID
	inline void AddValue(chrid cid, const T& val) { this->AddVal(cid, ChromData<T>(val)); }

#ifdef	_READDENS
	// Sets common chromosomes as 'Treated' in both of this instance and given Chroms.
	// Objects in both collections have to have second field as Treated.
	//	@obj: compared Chroms object
	//	@printWarn: if true print warning - uncommon chromosomes
	//	@throwExcept: if true then throw exception if no common chroms finded, otherwise print warning
	//	return: count of common chromosomes
	chrid	SetCommonChroms(Chroms<T>& obj, bool printWarn, bool throwExcept)
	{
		typename ChromMap<ChromData<T> >::Iter it;
		chrid commCnt = 0;

		// set treated chroms in this instance
		for (it = this->Begin(); it != this->End(); it++)	// 'this' required in gcc
			if (it->second.Treated = obj.FindChrom(CID(it)))
				commCnt++;
			else if (printWarn)
				Err(Chrom::Absent(CID(it), "second file")).Warning();
		// set false treated chroms in a compared object
		for (it = obj.Begin(); it != obj.End(); it++)
			if (!FindChrom(CID(it))) {
				it->second.Treated = false;
				if (printWarn)
					Err(Chrom::Absent(CID(it), "first file")).Warning();
			}
		if( !commCnt )	Err("no common " + Chrom::Title(true)).Throw(throwExcept);
		return commCnt;
	}
#endif	// _READDENS
};

// 'ItemIndexes' representes a range of chromosome's features/reads indexes,
struct ItemIndexes
{
	chrlen	FirstInd;	// first index in items container
	chrlen	LastInd;	// last index in items container

	inline ItemIndexes(chrlen first=0, chrlen last=1) : FirstInd(first), LastInd(last-1) {}
		
	// Returns count of items
	inline size_t ItemsCount() const { return LastInd - FirstInd + 1; }
};

template <typename I>
class Items : public Chroms<ItemIndexes>
{
public:
	typedef typename vector<I>::const_iterator cItemsIter;

protected:
	typedef typename vector<I>::iterator ItemsIter;

	vector<I> _items;	// vector of bed items

	// Applies function fn to each of the item for chrom defined by cit
	void ForChrItems(cIter cit, function<void(cItemsIter)> fn) const {
		const auto itEnd = ItemsEnd(Data(cit));
		for (auto it = ItemsBegin(Data(cit)); it != itEnd; it++)
			fn(it);
	}

	// Applies function fn to each of the item for chrom defined by cit
	//void ForChrItems(const value_type& c, function<void(cItemsIter)> fn) const {
	//	//const ItemIndexes& data = 
	//	const auto itEnd = ItemsEnd(c.Data.
	//		//Data(cit));
	//	for (auto it = ItemsBegin(Data(cit)); it != itEnd; it++)
	//		fn(it);
	//}

	// Applies function fn to each of the item for all chroms
	void ForAllItems(function<void(cItemsIter)> fn) const {
		for (const auto& c : Container()) {
			const auto itEnd = ItemsEnd(c.second.Data);
			for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
				fn(it);
		}
	}

	// Initializes size of positions container.
	inline void ReserveItems(size_t size) { _items.reserve(size); }

	// Gets item.
	//	@it: chromosome's iterator
	//	@iInd: index of item
	const I& Item(cIter it, chrlen iInd) const { return _items[Data(it).FirstInd + iInd]; }

	// Returns item.
	//	@cInd: index of chromosome
	//	@iInd: index of item
	//inline const I& Item(chrid cID, chrlen iInd) const { return Item(GetIter(cID), iInd); }

#ifdef _DEBUG
	inline void PrintItem(chrlen itemInd) const { _items[itemInd].Print(); }
#endif

public:
	// Gets total count of items
	inline size_t ItemsCount() const { return _items.size(); }

	// Gets count of items for chrom
	//	@data: item indexes
	inline size_t ItemsCount(const ItemIndexes& data) const { return data.ItemsCount(); }

	// Gets count of items for chrom
	//	@cit: chrom's iterator
	size_t ItemsCount(cIter cit) const { return ItemsCount(Data(cit)); }

	// Gets count of items for chrom; used in isChIP
	//	@cit: chrom's ID
	inline size_t ItemsCount(chrid cID) const { return ItemsCount(GetIter(cID)); }

	void PrintEst(ULONG estCnt) const {	cout << " est/fact: " << double(estCnt) / ItemsCount() << LF; }

#ifdef _ISCHIP
	// Prints items name and count, adding chrom name if the instance holds only one chrom
	//	@ftype: file type
	//	@prLF: if true then print line feed
	void PrintItemCount(FT::eType ftype, bool prLF = true) const
	{
		size_t iCnt = ItemsCount();
		dout << iCnt << SPACE << FT::ItemTitle(ftype, iCnt > 1);
		if (ChromCount() == 1)		dout << " per " << Chrom::TitleName(CID(cBegin()));
		if (prLF)	dout << LF;
	}
#endif
	// Returns a constan iterator referring to the first item of specified chrom
	//	@data: item indexes
	cItemsIter ItemsBegin(const ItemIndexes& data) const { return _items.begin() + data.FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@data: item indexes
	cItemsIter ItemsEnd(const ItemIndexes& data) const { return _items.begin() + data.LastInd + 1; }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	inline cItemsIter ItemsBegin(cIter cit) const { return ItemsBegin(Data(cit)); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cID: chromosome's ID
	inline cItemsIter ItemsBegin(chrid cID) const { return ItemsBegin(GetIter(cID)); }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome'sconstant  iterator
	inline cItemsIter ItemsEnd(cIter cit) const { return ItemsEnd(Data(cit)); }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cID: chromosome's ID
	inline cItemsIter ItemsEnd(chrid cID) const { return ItemsEnd(GetIter(cID)); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@data: item indexes
	ItemsIter ItemsBegin(ItemIndexes& data) { return _items.begin() + data.FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@data: item indexes
	ItemsIter ItemsEnd(ItemIndexes& data) { return _items.begin() + data.LastInd + 1; }

	// Returns an iterator referring to the first item of specified chrom
	//	@cit: chromosome's iterator
	inline ItemsIter ItemsBegin(Iter cit) { return ItemsBegin(Data(cit)); }

	// Returns an iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome's iterator
	inline  ItemsIter ItemsEnd(Iter cit) { return ItemsEnd(Data(cit)); }

#ifdef _DEBUG
	// Prints collection
	//	@title: item title
	//	@prICnt: number of printed items per each chrom, or 0 if all
	void Print(const char* title, size_t prICnt = 0) const
	{
		cout << LF << title << ": ";
		if (prICnt)	cout << "first " << prICnt << " per each chrom\n";
		else		cout << ItemsCount() << LF;
		for (const auto& inds : Container()) {
			const string& chr = Chrom::AbbrName(inds.first);
			const auto& data = inds.second.Data;
			const size_t lim = prICnt ? prICnt + data.FirstInd : UINT_MAX;
			for (chrlen i = data.FirstInd; i <= data.LastInd; i++) {
				if (i >= lim)	break;
				cout << chr << TAB;
				PrintItem(i);
			}
		}
	}
#endif	// _DEBUG
};

#ifdef _FEATURES

struct Featr : public Region
{
	float	Value;			// features's score

	inline Featr(const Region& rgn, float val = 0) : Value(val), Region(rgn) {}

	//inline Region& operator = (const Featr& f) { *this = f; }
#ifdef _DEBUG
	inline void Print() const { cout << Start << TAB << End << TAB << Value << LF; }
#endif	// _DEBUG
};

// 'Features' represents a collection of crhoms features
class Features : public Items<Featr>
{
	FBedInFile* _file = nullptr;		// valid only in constructor!
#ifdef _ISCHIP
	readlen	_minFtrLen;			// minimal length of feature
	float	_maxScore = 0;		// maximal feature score after reading
	bool	_uniScore = false;	// true if score is undefined in input data and set as 1
#elif defined _BIOCC
	// this is needed to get warning by calling Reads instead of Features (without -a option)
	bool	_narrowLenDistr = false;	// true if features length distribution is degenerate
#endif

	// Gets item's title
	// Obj abstract method implementation.
	//	@pl: true if plural form
	inline const string& ItemTitle(bool pl = false) const { return FT::ItemTitle(FT::eType::BED, pl); }

	// Gets a copy of Region by item's iterator.
	// Abstract BaseItems<> method implementation.
	inline Region const Regn(cItemsIter it) const { return Region(*it); }

	// Adds chrom to the instance
	//	@cID: chrom
	//	@cnt: count of chrom items
	void AddChrom(chrid cID, size_t cnt);

#ifdef _ISCHIP
	// Scales defined score through all features to the part of 1.
	void ScaleScores();
#endif

public:
#ifdef _ISCHIP
	// Creates new instance by bed-file name
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@joinOvrl: if true then join overlapping features, otherwise omit
	//	@scoreInd: index of 'score' field
	//	@bsLen: length of binding site: shorter features would be omitted
	//	@prfName: true if file name should be printed unconditionally
	Features(const char* fName, ChromSizes& cSizes, bool joinOvrl,
		BYTE scoreInd, readlen bsLen, bool prfName)
		: _minFtrLen(bsLen), _uniScore(!scoreInd)
	{
		FBedInFile file(fName, &cSizes, scoreInd, 
			joinOvrl ? UniBedInFile::eAction::JOIN : UniBedInFile::eAction::OMIT, 
			eOInfo::LAC, prfName, true);
#else
	// Creates new instance by bed-file name
	//	@fName: name of bed-file
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@joinOvrl: if true then join overlapping features, otherwise omit
	//	@prfName: true if file name should be printed unconditionally
	//	@abortInvalid: true if invalid instance should abort excecution
	Features(const char* fName, ChromSizes & cSizes, bool joinOvrl,
		eOInfo oinfo, bool prfName, bool abortInvalid = true)
	{
		FBedInFile file(fName, &cSizes, 5,
			joinOvrl ? UniBedInFile::eAction::JOIN : UniBedInFile::eAction::OMIT,
			oinfo, prfName, abortInvalid);
#endif
		size_t estItemCnt = file.EstItemCount();
		if (estItemCnt) {
			ReserveItems(estItemCnt);
			_file = &file;
			file.Pass(*this);
			_file = nullptr;
		}
#ifdef _BIOCC
		_narrowLenDistr = file.NarrowLenDistr();
#endif
		//PrintEst(estItemCnt);
	}


	// treats current item
	//	return: true if item is accepted
	bool operator()();

	// Closes current chrom, open next one
	//	@cID: current chrom ID
	//	@cLen: chrom length
	//	@cnt: current chrom items count
	//	@nextcID: next chrom ID
	inline void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) { AddChrom(cID, cnt); }

	// Closes last chrom
	//	@cID: last chrom ID
	//	@cLen: chrom length
	//	@cnt: last chrom items count
	//	@tCnt: total items count
	inline void operator()(chrid cID, chrlen cLen, size_t cnt, ULONG tCnt) { AddChrom(cID, cnt); }

	// Gets chromosome's feature by ID
	//	@cID: chromosome's ID
	//	@fInd: feature's index, or first feature by default
	//inline const Featr& Feature(chrid cID, chrlen fInd=0) const { return Item(cID, fInd); }

	// Gets chromosome's feature by iterator
	//	@it: chromosome's iterator
	//	@fInd: feature's index, or first feature by default
	inline const Featr& Feature(cIter it, chrlen fInd = 0) const { return Item(it, fInd); }

	// Gets chromosome's feature by iterator
	//	@it: chromosome's iterator
	//	@fInd: feature's index, or first feature by default
	inline const Region& Regn(cIter it, chrlen fInd = 0) const { return (const Region&)Item(it, fInd); }

	// Gets the sum length of all chromosome's features
	//	@it: chromosome's iterator
	chrlen FeaturesLength(cIter it) const;

	// Gets chromosome's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@it: chromosome's iterator
	//	@multiplier: 1 for numerics, 0 for letters
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	chrlen EnrRegnLength(cIter it, BYTE multiplier, float fLen) const {
		return (FeaturesLength(it) + chrlen(2 * fLen) * Data(it).ItemsCount()) << multiplier;
	}

	// Gets chrom's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@cID: chromosome's ID
	//	@multiplier: 1 for numerics, 0 for nameds
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	//	return: chrom's total enriched regions length, or 0 if chrom is absent
	chrlen EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const;

	// Return min feature length
	chrlen GetMinFeatureLength() const;

	// Return min distance between features boundaries
	chrlen GetMinDistance() const;

	// Expands all features positions on the fixed length in both directions.
	// If extended feature starts from negative, or ends after chrom length, it is fitted.
	//	@extLen: distance on which Start should be decreased, End should be increased,
	//	or inside out if it os negative
	//	@cSizes: chrom sizes
	//	@action: action for overlapping features
	//	return: true if positions have been changed
	bool Extend(chrlen extLen, const ChromSizes& cSizes, UniBedInFile::eAction action);

	// Checks whether all features length exceed given length, throws exception otherwise.
	//	@len: given control length
	//	@lenDefinition: control length definition to print in exception message
	//	@sender: exception sender to print in exception message
	void CheckFeaturesLength(chrlen len, const string& lenDefinition, const char* sender) const;

#ifdef _ISCHIP
	inline bool IsUniScore() const { return _uniScore; }
#else
	// Copies features coordinates to external DefRegions.
	void FillRegions(chrid cID, Regions& regn) const;
#endif	// _ISCHIP
#ifdef _BIOCC
	// Returns true if features length distribution is degenerate
	inline bool NarrowLenDistr() const { return _narrowLenDistr; }

	friend class JointedBeds;	// to access GetIter(chrid)
#endif
#ifdef _DEBUG
	void Print(size_t cnt = 0) const { Items::Print("features", cnt); }
#endif	// _DEBUG
};

#endif	// _FEATURES

// 'ChromSize' represents real and defined effective chrom lengths
struct ChromSize
{
	chrlen Real;				// real (actual) chrom length
#ifdef _ISCHIP
	mutable chrlen Defined = 0;	// defined effective chrom length;
								// 'effective' means double length for autosomes, single one for somatic
	
	// Sets chrom's effective (treated) real length as defined
	inline chrlen SetEffDefined(bool autosome) const { return bool(Defined = (Real << int(autosome))); }
#endif

	inline ChromSize(chrlen size = 0) : Real(size) {}
};

// 'ChromSizes' represented chrom sizes with additional file system binding attributes
// Holds path to reference genome and to service files
class ChromSizes : public Chroms<ChromSize>
{
	string	_ext;			// FA files real extention; if empty then instance is initialized by service dir
	string	_gPath;			// ref genome path
	string	_sPath;			// service path
	mutable genlen _gsize;	// size of whole genome


	// Returns length of common prefix before abbr chrom name of all file names
	//	@fName: full file name
	//	@extLen: length of file name's extention
	//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
	static int	CommonPrefixLength(const string& fName, BYTE extLen);

	// Initializes chrom sizes from file
	void Read(const string& fName);

	// Saves chrom sizes to file
	//	@fName: full file name
	void Write(const string& fName) const;

	// Fills external vector by chrom IDs relevant to file's names found in given directory.
	//	@cIDs: filling vector of chrom's IDs
	//	@gName: path to reference genome
	//	return: count of filled chrom's IDs
	chrid GetChromIDs(vector<chrid>& cIDs, const string& gName);

	// Initializes the paths
	//	@gPath: reference genome directory
	//	@sPath: service directory
	//	@prMsg: true if print message about service fodler and chrom.sizes generation
	void SetPath(const string& gPath, const char* sPath, bool prMsg);

	// returns true if service path is defined
	inline bool IsServAvail() const { return _sPath.size(); }

protected:
	inline chrlen Length(cIter it) const { return Data(it).Real; }

public:
	inline const string& RefExt() const { return _ext; }
	
	// Creates and initializes an instance.
	//	@gName: reference genome directory or chrom.sizes file
	//	@customChrOpt: id of 'custom chrom' option
	//	@prMsg: true if print message about service fodler and chrom.sizes generation
	//	@sPath: service directory
	//	checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
	ChromSizes(const char* gName, BYTE customChrOpt, bool prMsg, const char* sPath = NULL, bool checkGRef = false);

	inline ChromSizes() { _ext = _gPath = _sPath = strEmpty; }

	// Initializes ChromSizes by SAM header
	void Init(const string& headerSAM);

	inline bool IsFilled() const { return Count(); }

	// Return true if chrom.sizes are defined explicitly, by user
	//inline bool IsExplicit() const { return !_gPath.length(); }		// false if path == strEmpty
	
	// Returns reference directory
	inline const string& RefPath() const { return _gPath; }

	// Returns service directory
	inline const string& ServPath() const { return _sPath; }

	// Returns service directory
	inline const bool IsServAsRef() const { return _gPath == _sPath; }

	// Returns full ref chrom name by chrom ID 
	inline const string RefName(chrid cid) const { return _gPath + Chrom::AbbrName(cid); }

	// Returns full service chrom name by chrom ID 
	inline const string ServName(chrid cid) const { return _sPath + Chrom::AbbrName(cid); }

	inline cIter begin() const { return cBegin(); }
	inline cIter end() const { return cEnd(); }

	//inline chrlen Size(const Iter::value_type& sz) const { return sz.second.Data.Real; }

	inline chrlen operator[] (chrid cID) const { return At(cID).Data.Real; }


	// Gets total size of genome.
	genlen GenSize() const;

#ifdef _DEBUG
	void Print() const;
#endif
};

#ifdef _ISCHIP

// 'ChromSizesExt' provides additional functionality for the ChromSizes
class ChromSizesExt : public ChromSizes
{
	mutable chrid	_treatedCnt;	// number of treated chromosomes

	// Gets chrom's effective (treated) real length: a double length for autosomes, a single somatic
	//	@it: ChromSizes iterator
	inline chrlen SetEffLength(cIter it) const { return Data(it).SetEffDefined(Chrom::IsAutosome(CID(it))); }

public:
	// Creates and initializes an instance.
	//	@gName: reference genome directory
	//	@customChrOpt: id of 'custom chrom' option
	//	@printMsg: true if print message about chrom.sizes generation (in case of reference genome)
	inline ChromSizesExt(const char* gName, BYTE customChrOpt, bool printMsg, const char* sPath)
		: _treatedCnt(0), ChromSizes(gName, customChrOpt, printMsg, sPath, true) {}
	
	// Gets chrom's defined effective (treated) length
	//	@it: ChromSizes iterator
	chrlen DefEffLength(cIter it) const;

	// Gets count of treated chromosomes.
	inline chrid TreatedCount() const	{ return _treatedCnt; }

	// Sets actually treated chromosomes according template and custom chrom
	//	@templ: template bed or NULL
	//	return: number of treated chromosomes
	chrid	SetTreated	(bool statedAll, const Features* const templ);
	
	// Prints threated chroms short names
	void	PrintTreatedChroms() const;
};

#endif	// _ISCHIP

// 'RefSeq' represented reference chromosome as an array of nucleotides
class RefSeq
{
private:
	chrid	_ID;			// chrom ID
	char*	_seq = NULL;	// the nucleotides buffer
	chrlen	_len,			// length of chromosome
			_gapLen;		// total length of gaps
	Region	_effDefRgn;		// effective defined region (except 'N' at the begining and at the end)

	// Initializes instance and/or chrom's defined regions
	//	@fName: file name
	//	@rgns: chrom's defined regions: ripe or new
	//	@fill: if true fill sequence and def regions, otherwise def regions only
	//	return: true if chrom def regions are stated
	bool Init(const string& fName, ChromDefRegions& rgns, bool fill);

public:
	static bool	LetGaps;	// if true then include gaps at the edges of the ref chrom while reading
	static bool	StatGaps;	// if true count sum gaps for statistic output

	inline ~RefSeq()	{ delete [] _seq; }

	// Gets chrom legth
	inline chrlen Length()	const { return _len; }

	// Gets Read on position or NULL if exceeding the chrom length
	//const char* Read(chrlen pos) const { return pos > _len ? NULL : _seq + pos;	}

	//const char* Read(chrlen pos, readlen len) const { return pos + len > Length() ? NULL : _seq + pos; }

	// Gets subsequence without exceeding checking 
	inline const char* Seq(chrlen pos) const { return _seq + pos; }

#if defined _ISCHIP || defined _VALIGN

	// Creates a stub instance (for sampling cutting)
	//	@len: chrom length
	inline RefSeq(chrlen len) : _ID(Chrom::UnID), _seq(NULL), _len(len), _gapLen(0)
	{ _effDefRgn.Set(0, len); }

	// Creates and fills new instance
	RefSeq(chrid cID, const ChromSizes& cSizes);

#endif
#ifdef _ISCHIP
	// Returns chrom ID
	inline chrid ID() const { return _ID; }

	// Gets count of nucleotides outside of defined region
	inline chrlen UndefLength() const { return Length() - _effDefRgn.Length(); }

	//inline float UndefLengthInPerc() const { return 100.f * (Length() - _effDefRgn.Length()) / Length(); }

	// Gets feature with defined region (without N at the begining and at the end) and score of 1
	inline const Featr DefRegion() const { return Featr(_effDefRgn); }
	
	// Gets start position of first defined nucleotide
	inline chrlen Start()	const { return _effDefRgn.Start; }
	
	// Gets end position of last defined nucleotide
	inline chrlen End()		const { return _effDefRgn.End; }
	
	// Gets total length of gaps
	inline chrlen GapLen()	const { return _gapLen; }

#elif defined _READDENS || defined _BIOCC

	// Creates an empty instance and fills chrom's defined regions
	//	@fName: FA file name with extension
	//	@rgns: new chrom's defined regions
	//	@minGapLen: minimal length which defines gap as a real gap
	inline RefSeq(const string& fName, ChromDefRegions& rgns, short minGapLen);

#endif	// _ISCHIP

//#if defined _FILE_WRITE && defined DEBUG 
//	// Saves instance to file by fname
//	void Write(const string & fname, const char *chrName) const;
//#endif
};

#if defined _READDENS || defined _BIOCC

// 'DefRegions' represents chrom's defined regions
// initialized from ChromSizes (statically, at once)
// or from .fa files (dynamically, by request).
class DefRegions : public Chroms<Regions>
{
	ChromSizes&		_cSizes;
	const chrlen	_minGapLen;	// minimal allowed length of gap
#ifdef _BIOCC
	const bool		_singleRgn = true;	// true if this instance has single Region for each chromosome
#endif


public:
	// Creates an instance by genome name, from chrom sizes file or genome.
	//	@cSizes: chrom sizes
	//	@minGapLen: minimal length which defines gap as a real gap
	inline DefRegions(ChromSizes& cSizes, chrlen minGapLen)
		: _cSizes(cSizes), _minGapLen(minGapLen)
	{ Init(); }

	void Init();

	// Returns true if regions are not initialized
	inline bool IsEmpty() const {	return !Count(); }

	inline ChromSizes& ChrSizes() { return _cSizes; }

	// Gets chrom's size by chrom's iterator
	inline chrlen Size(cIter it)	const { return Data(it).LastEnd(); }

	// Gets chrom's size by chrom's ID
	inline chrlen Size(chrid cID)	const { return At(cID).Data.LastEnd(); }

	// Gets chrom regions by chrom ID; lazy for real chrom regions
	const Regions& operator[] (chrid cID);

#ifdef _BIOCC
	// Copying constructor: creates empty copy!
	inline DefRegions(const DefRegions& gRgns) :
		_cSizes(gRgns._cSizes), _minGapLen(gRgns._minGapLen), _singleRgn(gRgns._singleRgn) {}

	// Gets miminal size of chromosome: for represented chromosomes only
	chrlen MinSize() const;

	// Gets total genome's size: for represented chromosomes only
	genlen GenSize() const;

	// Adds chromosomes and regions without check up
	inline void AddChrom (chrid cID, const Regions& rgns) {	AddVal(cID, rgns); }
#endif	// _BIOCC

#ifdef _DEBUG
	void Print() const;
#endif
};
#endif	// _READDENS || _BIOCC

#if defined _ISCHIP || defined _CALLDIST

typedef pair<float, float> fpair;

// 'LenFreq' represents a fragment's/read's length frequency statistics
class LenFreq : map<fraglen,ULONG>
{
public:
	// combined type of distribution
	enum /*class*/ eCType {		// not class to have a cast to integer by default
		NORM = 1 << 0,
		LNORM = 1 << 1,
		GAMMA = 1 << 2,
		CNT = 3,
	};

private:
	using dtype = int;	// consecutive distribution type: just to designate dist type, used as an index
	using spoint = pair<unsigned, ULONG>;	// initial sequence point
	using point = pair<unsigned, float>;	// distribution point 

	// Returns combined distribution type by consecutive distribution type
	inline static eCType GetCType(dtype type) { return eCType(1 << type); }

	// Returns consecutive distribution type by combined distribution type
	const static dtype GetDType(eCType ctype) { return RightOnePos(int(ctype)); }

	enum class eSpec {	// distribution specification
		CLEAR,		// normal quality;	exclusive
		SMOOTH,		// complementary
		MODUL,		// modulated;	complementary
		EVEN,		// exclusive
		CROP,		// cropped;	exclusive
		HCROP,		// heavily cropped;		exclusive
		SDEFECT,	// slightly defective;	exclusive
		DEFECT		// defective; exclusive
	};

	static const char* sTitle[];
	static const string sSpec[];
	static const string sParams;
	static const string sInaccurate;
	const fraglen smoothBase = 1;	// splining base for the smooth distribution
	static const int hRatio = 2;	// ratio of the summit height to height of the measuring point
	static const float lghRatio;	// log of ratio of the summit height to height of the measuring point

	// Moving window (Sliding subset)
	class MW : protected vector<ULONG>
	{
	protected:
		// Constructor by half of moving window length ("base")
		MW(unsigned base) { insert(begin(), Size(base), 0); }

		// Add last value and pop the first one (QUEUE functionality)
		void PushVal(ULONG x);

	public:
		// Returns length of moving window
		//	@base: moving window half-length
		static fraglen Size(fraglen base) { return 2 * base + 1; }
	};

	// Simple Moving Average splicer
	// https://en.wikipedia.org/wiki/Moving_average
	class MA : public MW
	{
		ULONG	_sum = 0;		// sum of adding values

	public:
		// Constructor by half of moving window length ("base")
		inline MA(fraglen base) : MW(base) {}

		// Add value and return average
		float Push(ULONG x);
	};

	// Simple Moving Median splicer
	class MM : public MW
	{
		vector<ULONG> _ss;		// sorted moving window (sliding subset)
		ULONG(MM::* _push)(ULONG) = &MM::GetMedian;	// pointer to GetMedian function: real or empty (stub)

		//	Return stub median
		inline ULONG GetMedianStub(ULONG x) { return x; }

		//	Add value and return median
		ULONG GetMedian(ULONG x);

	public:
		// Constructor
		//	@base: moving window half-length;
		//	if 0 then empty instance (initialized by empty Push() method)
		MM(fraglen base);

		// Add value and return median
		inline ULONG Push(ULONG x) { return (*this.*_push)(x); }
	};

	// Keeps distribution params: PCC, mean(alpha), sigma(beta)
	struct DParams
	{
	private:	
		static const float UndefPCC;
	public:
		float	PCC = 0;			// Pearson correlation coefficient
		fpair	Params{ 0,0 };		// mean(alpha), sigma(beta)

		bool operator >(const DParams& dp) const { return PCC > dp.PCC; }

		inline bool IsUndefPcc() const { return PCC == UndefPCC; };

		inline void SetUndefPcc() { PCC = UndefPCC; };
	};

	// 'AllDParams' represents a collection of restored distribution params for all type of distribution
	class AllDParams
	{
		// 'QualDParams' keeps restored distribution parameters
		struct QualDParams
		{
			eCType	Type;			//  combined distribution type
			const char* Title;		// distribution type title
			DParams dParams;		// PCC, mean(alpha), sigma(beta)

			// Sets combined type and title by consecutive type
			void SetTitle(dtype type) { Type = GetCType(type); Title = sTitle[type]; }

			// Returns true if distrib parameters set
			inline bool IsSet() const { return dParams.PCC != 0; }

			// Prints restored distr parameters
			//	@s: print stream
			//	@maxPCC: masimum PCC to print relative PCC percentage
			void Print(dostream& s, float maxPCC) const;
		};

		array<QualDParams, eCType::CNT>	_allParams;
		bool _sorted = false;

		// Returns true if distribution parameters set in sorted instance
		bool IsSetInSorted(eCType ctype) const;

		// Returns number of distribution parameters set in sorted instance
		int SetCntInSorted() const;

		// Returns DParams by combined distribution type
		inline DParams& Params(eCType ctype) { return _allParams[GetDType(ctype)].dParams; }

		// Sorts in PCC descending order
		void Sort();

	public:
		// Default constructor
		AllDParams();

		// Set distribution parameters by type
		//	@type: consecutive distribution type
		//	@dp: PCC, mean(alpha) & sigma(beta)
		inline void SetParams(dtype type, const DParams& dp) { _allParams[type].dParams = dp; }

		// Clear normal distribution if its PCC is less then lognorm PCC by the threshold
		void ClearNormDistBelowThreshold(float thresh) {
			if (Params(eCType::LNORM).PCC / Params(eCType::NORM).PCC > thresh)
				Params(eCType::NORM).PCC = 0;
		}

		// Returns parameters of distribution with the highest (best) PCC
		//	@dParams: returned PCC, mean(alpha) & sigma(beta)
		//	return: consecutive distribution type with the highest (best) PCC
		dtype GetBestParams(DParams& dParams);

		// Prints sorted distibutions params
		//	@s: print stream
		void Print(dostream& s);
	};
	
	// Returns specification string by specification type
	inline static const string Spec(eSpec s) { return "Distribution " + sSpec[int(s)]; }
	
	// Returns true if consecutive type is represented in combo cType
	inline static bool IsType(eCType cType, dtype type) { return cType & (1 << type); }

	// Returns true if combo type is represented in combo cType
	inline static bool IsType(eCType cType, eCType type) { return cType & type; }

	// Calls distribution parameters by consecutive distribution type
	//	@keypts: key points: X-coord of highest point, X-coord of right middle hight point
	//	@params: returned mean(alpha) & sigma(beta)
	//	Defined as a member of the class only to use the private short name lghRatio
	static void (*SetParams[])(const fpair& keypts, fpair& params);

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	AllDParams	_allParams;			// distributions parameters
#ifdef _DEBUG
	mutable vector<point> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Returns estimated moving window half-length ("base")
	//	return: estimated base, or 0 in case of degenerate distribution
	fraglen GetBase();

	// Defines key points
	//	@base: moving window half-length
	//	@summit: returned X,Y coordinates of spliced (smoothed) summit
	//	return: key points: X-coord of highest point, X-coord of right middle hight point
	fpair GetKeyPoints(fraglen base, point& summit) const;

	// Compares this sequence with calculated one by given mean & sigma, and returns PCC
	//	@type: consecutive distribution type
	//	@dParams: returned PCC, input mean(alpha) & sigma(beta)
	//	@Mode: X-coordinate of summit
	//	@full: if true then correlate from the beginning, otherwiase from summit
	//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
	void CalcPCC(dtype type, DParams& dParams, fraglen Mode, bool full = true) const;

	// Calculates distribution parameters
	//	@type: consecutive distribution type
	//	@keyPts: key points: X-coord of highest point, X-coord of right middle hight point
	//	@dParams: returned PCC, mean(alpha) & sigma(beta)
	//	@Mode: X-coordinate of summit
	void SetPCC(dtype type, const fpair& keypts, DParams& dParams, fraglen Mode) const;

	// Calculates and print called distribution parameters
	//	@type: consecutive distribution type
	//	@base: moving window half-length
	//	@summit: returned X,Y coordinates of best spliced (smoothed) summit
	void CallParams(dtype type, fraglen base, point& summit);

	// Prints original distribution features
	//	@s: print stream
	//	@base: moving window half-length
	//	@summit: X,Y coordinates of spliced (smoothed) summit
	void PrintTraits(dostream& s, fraglen base, const point& summit);

	// Prints original distribution as a set of <frequency>-<size> pairs
	//	@s: print stream
	void PrintSeq(dostream& s) const;

public:
	// Default constructor
	LenFreq() {}

	// Constructor by pre-prepared frequency distribution file
	//	@fname: name of pre-prepared frequency distribution file
	LenFreq(const char* fname);

	// Returns true if distribution has not enough size
	//inline bool IsDegenerate() const { return size() < 5; }

	// Adds fragment/read to statistics
	//	@len: frag's length
	inline void AddLen(fraglen len) { (*this)[len]++; }

	// Calculate and print dist fpair
	//	@s: print stream
	//	@type: combined type of distribution
	//	@prDistr: if true then print distribution additionally
	void Print(dostream& s, eCType type, bool prDistr = true);
};

#endif	// _ISCHIP
#if defined _ISCHIP || defined _BSDEC

using coval = UINT;					// coverage value
using covmap = map<chrlen, coval>;	// coverage map

// 'AccumCover' represents cumulative chrom's fragment coverage data
//	and implements a single method for gradual filling (incrementing) coverage
class AccumCover : public covmap
{
	bool _unsaved = true;	// true if data is still unsaved

public:
	bool Closed = false;

	// Default constructor
	inline AccumCover() {}

	// Copy constructor
	inline AccumCover(const AccumCover& cv) : covmap(cv) { }

	// Returns true if data is unsaved
	inline bool Unsaved() const { return _unsaved; }

	// Clears data and mark it as saved
	void Clear() { _unsaved = false; clear(); }

	// Adds fragment to accumulate the coverage
	void AddRegion(const Region& frag);

	// Calls functor for each point that put the chrom coverage
	template<typename Functor>
	void DoWithItem(Functor f) const {
		for (const value_type& item : *this)	f(item);
	}

	// Calls functor for each point that put the chrom coverage
	template<typename Functor>
	void DoWith2Items(Functor f) const {
		auto it0 = cbegin(), it = it0;		// previous, current entry

		for (++it; it != end(); it0 = it++)
			if (it0->second)	f(it0, it);
	}

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

#endif	// _ISCHIP || _BSDEC
#endif	// _DATA_H