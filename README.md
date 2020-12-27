# isChIP
***I**n **S**ilico* **ChIP**-seq is a fast realistic ChIP-seq simulator.

The modelling of the chromatin immunoprecipitaion followed by next generation sequencing process is based primarily on Illumina protocol. 
In addition, extra flexibility is implemented for the isChIP’s options to be straightforwardly re-formulated in other formats, 
such as Ion Torrent, Roche454, etc.<br>
More details of the model are provided in the section [Model: brief description](#model-brief-description).<br>
Suitable for [single cell simulation](#example-of-single-cell-simulation).

### Performance
On 2.5 GHz RAID HPC by default values of ground samples, in one-thread mode, within 1 minute **isChIP** records:<br>
in *test* mode up to 2 million reads (33 Kread/sec) with size selection, 
up to 10 million reads (167 Kread/sec) without size selection (both values by default background level);<br>
in *control* mode up to 66 million reads (1.1 Mread/sec).<br>
Applied amplification increases these values from 1.5 to 3 times.<br>
The first run is a bit slower because of creating service files.<br>
The required memory is linearly proportional to the number of threads. For one thread, it does not exceed 300 Mb.

### Navigation:
  [Installation](#installation)<br>
  [Usage](#usage)<br>
  [Details](#details)<br>
  [Model: brief description](#model-brief-description)<br>
  [Fragment distribution and size selection](#fragments-distribution-and-size-selection)<br>
  [Example of single cell simulation](#example-of-single-cell-simulation)

## Installation
### Executable file

[Download Linux version](https://github.com/fnaumenko/isChIP/releases/download/v1.0/isChIP-Linux-x64.gz)<br>
[Download Windows version](https://github.com/fnaumenko/isChIP/releases/download/v1.0/isChIP-Windows-x64.zip)

Alternative download in Linux:<br>
`wget -O isChIP.gz https://github.com/fnaumenko/isChIP/releases/download/v1.0/isChIP-Linux-x64.gz`

Then type<br>
```
gzip -d isChIP.gz
chmod +x isChIP
```

### Compiling in Linux
Required libraries:<br>
g++<br>
pthread<br>
zlib (optionally)

To compile from Git, type:
```
git clone https://github.com/fnaumenko/isChIP
cd isChIP
make
```
Alternative:
```
wget -O isChIP.zip https://github.com/fnaumenko/isChIP/archive/v1.0.zip
unzip isChIP.zip
cd isChIP-1.0
make
```
If **zlib** is not installed on your system, the program will be compiled without the ability to read/write compressed files.

### Prepare reference genome
Download the required reference genome from UCSC: *ftp://hgdownload.soe.ucsc.edu/goldenPath/*<br>
For example, to download human genome library **hg19**:<br>
**in Linux**:<br>
```
mkdir hg19
cd hg19
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/ ./
```
Alternative:<br>
```
wget -r ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/
mv hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes hg19
rm -r hgdownload.soe.ucsc.edu
```
Unplaced (chr\*\_random) and unlocalized (chrUn_\*\) sequences are not involved in modelling, 
so you can delete them by typing `rm hg19/*_*`<br><br>
**in Windows**<br>
copy and paste the string *ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/* into Windows browser address bar, 
then copy *.fa.gz files to your local directory.<br>
Alternatively use FTP client, e.g. [FileZilla](https://filezilla-project.org/).

To download **mm9** library, in the above example replace '**hg19**' with '**mm9**'.

## Usage
```
isChIP [options] -g|--gen <name> [<template>]
	<template> - bed file whose features specify binding sites
```

### Synopsis
`isChIP -g ref_genome –n 5`<br>
generates 'input' sequences in FastQ with read length of 50 and average read density of about 9 read/kbs, 
comparable with what is experimentally observed.

`isChIP -g mm9_dir –tz –n 300 –r 36 –f fq,sam templ.bed`<br>
generates test sequences in zipped FastQ and SAM, with read length of 36, timing, and average foreground/background read density 
comparable with what is experimentally observed in [Series GSE56098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56098).

### Help
```
Processing:
  -g|--gen <name>       reference genome library. Required
  -n|--cells <int>      number of nominal cells [1]
  -G|--ground <[float]:[float]>
                        fore- and background levels:
                        number of selected fragments inside/outside binding sites,
                        in percent of all/foreground.
                        In control mode background is ignored [50:1]
  -E|--exo [<int>]      apply ChIP-exo protocol modification with specified exonuclease 'headroom' length.
                        If the option is not specified, ChIP-exo is not applied [6]
  -D|--mda              apply MDA technique
  -a|--pcr <int>        number of PCR cycles [0]
  -c|--chr <name>       generate output for the specified chromosome only
  --bg-all <OFF|ON>     turn on/off generation background for all chromosomes.
                        For the test mode only [ON]
  -m|--smode <SE|PE>    sequencing mode: SE - single end, PE - paired end [SE]
  -s|--bscore <int>     index of the template field used to score each binding event.
                        Value '0' means ignore template scores. For the test mode only [5]
  --edge-len <int>      unstable binding length (BS edge effect). For the test mode only [0]
  -N|--full-gen         process the entire reference chromosomes (including marginal gaps)
  -P|--threads <int>    number of threads [1]
  --serv <name>         folder to store service files [-g|--gen]
  --seed <int>          fix random emission with given seed, or 0 if don't fix [0]
Fragment size distribution:
  -L|--ln <[float]:[float]>
                        mean and stand dev of fragment lognormal distribution [5.26:0.3]
  -S|--ss [<[int]:[int]>]
                        apply size selection normal distribution with specified mean and stand dev.
                        If the option is not specified, size selection is disabled [auto:30]
Reads:
  -r|--rd-len <int>     fixed length of output read, or minimum length of variable reads [50]
  -R|--rd-dist [<[int]:[int]>]
                        mean and stand dev of variable read normal distribution,
			according to Ion Torrent/Roche454 protocol.
                        If the option is not specified, the read length is fixed [200:20]
  --rd-pos              add read position to its name
  --rd-Nlim <int>       maximum permitted number of ambiguous code N in read [OFF]
  --rd-lim <long>       maximum permitted number of total recorded reads [2e+08]
  --rd-ql <char>        uniform read quality score [~]
  --rd-ql-patt <name>   read quality scores pattern
  --rd-ql-map <int>     read mapping quality in SAM and BED output [255]
Output:
  -f|--format <FQ,BED,SAM,BG,FDENS,RDENS,FDIST,RDIST>
                        format of output data, in any order [FQ]
  -C|--control          generate control simultaneously with test
  -x|--strand           generate two additional wig files, each one per strand
  -T|--sep              use 1000 separator in output
  -o|--output <name>    location of output files or existing folder
                        [TEST mode: mTest.*, CONTROL mode: mInput.*]
  -z|--gzip             compress the output
Other:
  -V|--verbose <SL|RES|RT|PAR|DBG>
                        set verbose level:
                        SL -    silent mode (show critical messages only)
                        RES -   show result summary
                        RT  -   show run-time information
                        PAR -   show actual parameters
                        DBG -   show debug messages [PAR]
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

## Details

### Mode

**isChIP** generates output in one of the two modes:<br>
*test* – simulation of site of interest sequencing; the output is test sequences/alignment<br>
*control* – simulation of control production; the output is 'input' sequences/alignment.<br>
These modes are distinguished only by involvement or elimination of sites of interest in a process. 
The sites are represented by the set of features stated in a single optional parameter – a BED file called *template*.

### Template
*Template* is a file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format whose features correspond to binding events. 
Each data line contains 3 (minimum required) or more fields. First 3 fields define the binding site. 
If there are no more fields in the line, all generating significant enriched regions have the same density. 
If there are additional fields, isChIP interprets the field indicated by the option `-s|--bscore` as the relative binding efficiency. 
Within the framework of the model, binding efficiency is defined as the probability of associating a fragment with a binding site. 
Just for simplicity, the term 'binding score' is used as a synonym for binding efficiency, although in a strict sense it is debatable. 
'Relative' means that the units specified in the field are not taken into account, only value ratios are implemented to be meaningful.<br>
BED format recommends using the 4th field as a 'name' and 5th field as a 'score' in the common sense. 
Many developers of peak callers follow these guidelines and save in the 4th field the peak name, and in the 5th field the peak p-value or local FDR. 
Both of the values are directly proportional to the peak amplitude, which makes it possible to treat them as binding efficiency. 
When using another field to store the corresponding values, it can be specified by the `-s|--bscore` option. 
That allows the use of peak caller output as a **isChIP** *template*, 
with the necessary replacement of the coordinates of the peaks to the coordinates of the proposed binding sites.<br>
Alternatively, the scores can be entered manually.<br>
At the same time different stated scores in *template* can be ignored by setting the `-s|--bscore` option to 0.

*Template* does not have to be sorted, but all the features belonging to a chromosome should be clustered together. 
The simplest way to meet this requirement is to use pre-sorted data.

### Options description

*Notes:*<br>
For options that take a pair of values, one or both of them can be omitted. 
In this case, it retains the default value. 
For example, `-G :5` means foreground level of 50% and background level of 5%; 
`–G :` does not change the settings.<br>
Enumerable option values are case insensitive.<br>
Compressed input files in gzip format are acceptable.

`-g|--gen <name>`<br>
Reference genome library. It is a directory contained nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
First **isChIP** searches .fa files in it. If there are not such files, 
or the file corresponded to chromosome specified by option `–c|--chr` is absent, **isChIP** searches .fa.gz files.<br>
In *test* mode the target references are determined by *template*. 
See [Template](#template) and `--bg-all` and `–c|--chr` options.<br>
**isChIP** omits  'random' contigs and haplotype sequences.
One can obtain a genome library in UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, 
e.g. unmasked (‘dna'), since program does not recognise mask’s types.<br>
The program also generates and saves chromosome sizes files in this directory on first launch (see `--serv <name>`).<br>
This option is required.

`-n|--cells <int>`<br>
specifies the number of nominal cells, it is usually required. 
Each nominal cell corresponds to one passage through the reference chromosome(dual for the autosomes).<br>
The convention is that this number assigns a greatest possible practical output by default 
(50% value of the foreground stated by `-G|--ground` option means 0% loss, which is unrealistic). 
If you want to deal with a real number of cells, you should assign a real level of overall loss. 
The important implication is that this proportionally increases the runtime. 
For example, assigning 100 nominal cells with zero losses leads to the generation of data equivalent to 100,000 ‘real’ cells with a ‘real’ total loss of 99.9%, 
but in a runtime about 1000 times less (e.g. 10 seconds instead of almost 3 hours per chromosome).

Since the total losses are difficult to estimate in advance, the inverse task, namely, 
determining the overall loss in a real experiment post-factum, looks more practical. 
To accomplish this, it is sufficient to select a model output with a peak density 
corresponding to the peak density in the actually observed data through several fast iterations 
with a different number of nominal cells and zero foreground level. 
Then the overall loss is calculated by the formula (1-n/N)*100%, where n is the number of nominal cells, N is the number of real cells in the experiment.

Besides the number of cells, the program’s run time heavily depends on the parameters of fragment size distribution and size selection 
(see [Fragment distribution and size selection](#fragments-distribution-and-size-selection)). 
In general, with the default fragment size distribution, a number of nominal cells from 1 to 5 
provides an output read’s mean density comparable with real ‘inputs’, 
and a value from 50 to 300 leads to generated simulated data, which is comparable with real tests in terms of density.<br>
Range: 1-2e6<br>
Default: 1

`-G|--ground <[float]:[float]>`<br>
specifies foreground and background levels.<br>
**In *test* mode:**<br>
foreground is defined as the number of selected fragments, which are intersected with the *template* binding sites, 
as a percentage of the total number of generated fragments intersecting with the *template* binding sites. 
Background is defined as the number of selected fragments, which are not intersected 
with the *template* binding sites, as a percentage of the foreground.<br>
In practice a background level of 1-3% corresponds with a good experimental data set, 
while a level of more than 6% leads to generate data with a rather low signal-to-noise ratio.<br>
**In *control* mode:**<br>
foreground is defined as the number of selected fragments, as a percentage of the total number of generated fragments. 
It should be considered as the level of overall loss.<br>
Background value is ignored.<br>
Range: 0-100 for both values<br>
Default: 50:1

`-E|--exo [<int>]`<br>
applies ChIP-exo protocol modification with specified *exonuclease 'headroom'*. 
This is the distance between the exonuclease cleavage site and the edge of the template binding site, 
at which exonucleases end their digestion with a probability exponentially decreasing from 1 to 0 
in the direction away from the binding site.<br>
This term is borrowed from the article 
[Comprehensive Genome-wide Protein-DNA Interactions Detected at Single-Nucleotide Resolution](https://www.cell.com/fulltext/S0092-8674(11)01351-1), 
as well as the default value for this parameter.<br>
Comparison of the model data with the data from the article is shown in the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/Exo-real-sim.png).<br>
The peaks detection should be done with specialized ChIP-exo peak-callers (see, for example, 
[Comparative analysis of ChIP-exo peak-callers: impact of data quality, read duplication and binding subtypes]( https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3403-3).<br>
Default (if applied): 6

`-D|--mda`<br>
applies [MDA](https://en.wikipedia.org/wiki/Multiple_displacement_amplification) technique.<br>
See [Model: brief description](#model-brief-description) for more details.<br>
Enabling this option leads to some increase in the run time.<br>
The process applies to all fragments (fore- and background).

`-a|--pcr <int>`<br>
specifies the number of [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction) amplification cycles.<br>
See [Model: brief description](#model-brief-description) for more details.<br>
If MDA is applied, PCR performs after it. The process applies to all fragments.<br>
Value 0 means the absence of amplification.<br>
Range: 0-100<br>
Default: 0

`-c|--chr <name>`<br>
Generate output for the specified chromosome only. 
The value `name` is the chromosome identifier; it is a number or character, for example, `10`, `X`.<br>
This is a strong option, which forces to ignore all other chromosomes from reference genome and *template*, 
and abolishes the impact of option `--bg-all`. 
If the specified chromosome is absent in *template*, the program has nothing to simulate.

`--bg-all <OFF|ON>`<br>
turns off/on generation of background for all chromosomes in *test* mode, whether or not they are presented in *template*.<br>
Mapping one chromosome to the whole reference genome leads to the appearance of short local discontinuities of alignment, 
‘gaps’, corresponding to low mappability regions. This observation holds true for all aligner. 
Typically, generation of background for all chromosomes eliminates this issue.<br>
Notably, this process is time consuming and is justified only in the case of subsequent use of the alignment, 
for example, for comparison with experimental FastQ or to validate aligners. 
When validating peak detectors on one or some chromosomes, it is more advantageous to use a direct SAM output format (see `-F|--format`).<br>
The established option `-c|--chr` makes this option negligible.<br>
In *control* mode is ignored.<br>
Default: `ON`

`-m|--smode <SE|PE>`<br>
generates reads according to stated sequencing mode: `SE` – single end, `PE` – paired end .<br>
Default: `SE`

`-s|--bscore <int>`<br>
specifies the 1-based index of the field in *template* which is used to score each binding event 
(in terms of binding efficiency, see [Template](#template)).<br>
Setting index to 0 leads to ignoring the specified scores: all binding events will be simulated with maximum efficiency.<br>
*Note:* specifying a field containing a string (name) or a missing field will result in zero efficiency without any warning.<br>
In *control* mode is ignored.<br>
Range: 4-12 (except for the zero value)<br>
Default: 5

`--edge-len <int>`<br>
specifies unstable binding length (binding site edge effect. 
By default fragment is considered linked to the binding site in the presence of at least one intersecting nucleotide. 
For various reasons, this connection may be unstable. 
This option specifies the length of the intersection of the fragment with the site, 
at which the probability of binding increases from 0 to 1 in direct proportion to the number of intersected nucleotides.<br>
When the option value is more than half the length of the minimum site, it decreases to this half.<br>
In *control* mode is ignored.<br>
Range: 0-10<br>
Default: 0

`-N|--full-size`<br>
forces to scan the entire reference chromosome, including gaps at each end (‘undefined telomeres’).<br>
By default these gaps are excluded from processing.<br>
This permission makes no difference in data after alignment, but slightly increases the output volume and the runtime. 
For instance, when simulation on the basis of mm9 genome, this option increases both values by about 2.4%. 
For the hg19 genome the difference is about 0.7%.

`-P|--threads <int>`<br>
specifies the number of threads. The workflow is separated between chromosomes, so the actual number of threads 
can be reduced (if the number of actual treated chromosomes is less then assigned value). 
The actual threads number is displayed in `PAR` and `DBG` verbose mode.<br>
Range: 1-20<br>
Default: 0

`--serv <name>`<br>
specifies the service directory – a place for keeping service files *chr\<x\>.region*, chromosome sizes file and sample files. 
The program generates these files on first launch, and then reuses them. 
By default, they are stored in the reference genome folder. 
But if this folder is closed for writing, or you want to store these files separately for your own reasons, that is the place.<br>
Default: reference genome directory

`--seed <int>`<br>
fixes random numbers emission to get repetitive results. 
The actual seed equals the option value increased by a certain factor to provides a noticeable difference in the of random number generation option values that differ by 1.<br>
Value 0 means non-recurring random generation.<br>
Range: 0-1000<br>
Default: 0

`-L|--ln <[float]:[float]>`<br>
specifies the mean and standard deviation of fragments sizes lognormal distribution.<br>
Default parameters give a lognormal mode (expected value) of 176 and a mean (average) of 200.<br>
When `-R|--rd-dist` option is enabled, a wide distribution is set by default with mean of 5.92 and a sigma of 0.4, 
which gives a lognormal mode of 317 and an average of 403.<br>
See [Fragment distribution and size selection](#fragments-distribution-and-size-selection).<br>
Range: 3:0.3-9:1<br>
Default: 5.26:0.3

`-S|--ss [<[int]:[int]>]`<br>
specifies the mean and standard deviation of fragment size selection normal distribution.<br>
By default size selection is deactivated. Specifying an option leads to activation of size selection with default values, 
if no option value is set, or with specified value(s).<br>
If option is specified without value(s), size selection mean value is calculated as the mode of the fragment’s lognormal distribution. 
In other words, this value is equal to the most frequent value of the lognormal distribution, which for given defaults is 200.<br>
See [Fragment distribution and size selection](#fragments-distribution-and-size-selection).<br>
Range: 50:2-2000:500 (except for the zero values)<br>
Default (if applied): auto:30

`-r|--rd-len <int>`<br>
specifies the fixed length of output read.<br>
Also see `-R|--rd-dist` option.<br>
Notably, **isChIP** does not perform a preliminary check of the read length relative to the parameters of fragments distribution or size selection. 
If the reads are too long, their number at the output can be reduced up to 0 without any warnings.<br>
Range: 20-500<br>
Default: 50

`-R|--rd-dist [<[int]:[int]>]`<br>
specifies mean and standard deviation of variable read normal distribution, according to Ion Torrent/Roche454 protocol.<br>
As can be seen in the ![Read distributions](https://github.com/fnaumenko/bioStat/tree/master/pict/Read_distrs.png), 
the real read distributions do not always follow the normal pattern. 
Nevertheless, it should be admitted as the best approximation.<br>
If this option is applied, `-r|--rd-len` option specifies the minimum limit of read length (20 bp by default).<br>
Enabling this option when MDA is activated leads to an additional increase in the run time.<br>
Default (if applied): 200:30

`--rd-pos`<br>
Adds read’s actual start position to its name in the output files. It is useful for verifying aligners.<br>
In case of PE reads start and end positions of fragment are added.<br>
Read’s name format:  \<app\>:\<chr\>[:pos|:start-end].\<numb\>[/mate]<br>
Mate is only added in BED output.

`--rd-Nlim <int>`<br>
Maximum permitted number of ambiguous code N in read. Reads exceeding this limit are ignored. 
The corresponding source fragments are also ignored for output formats accumulating fragments (`BG`, `FDENS`, `FDIST`; see `-f|--format`).<br>
This option is designed primarily for aligners validation.<br>
In the real read sequences ambiguous code could cause a failure in sequencing as well as undefined regions in the cut fragments, 
while simulated reads are fully defined by reference library. 
Consequently the ambiguous bases in the library are translated into reads. 
Different aligners considered these codes differently: some of them as invalid, some of them as case of overlapping.<br>
Range: 0 - <`-r|--rd-len` value><br>
Default: `-r|--rd-len` value (which means no checkup)

`--rd-lim <long>`<br>
specifies the maximum number of total written reads. The value emulates sequencer’s limit. 
It restricts the number of recorded reads for each chromosome proportionally 
(and the number of their source fragments for `BG`, `FDENS`, `FDIST` output formats as well).<br>
Range: 1e5-1e19<br>
Default: 2e8.

`--rd-ql <char>`<br>
provides a uniform read quality score in FQ and SAM output.<br>
Range: '!'-'\~'<br>
Default: '\~' (decimal 126)

`--rd-ql-patt <name>`<br>
specifies a pattern of the read quality values in FQ and SAM output. 
Parameter `name` is a name of a plain text file containing at least one line encoding the quality values 
for the sequence as described in [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) specification.<br>
If the length of line is less then read’s actual (fixed or variable) length, 
then the rest of pattern is filed by value defined by `--rd-ql` option.<br>
Only the first uncomment line is taken into account.<br>
Сomment lines are lines starting with the character '#'.<br>

`-rd-ql-map <int>`<br>
provides a read mapping quality in SAM and BED output (in the last case it is also called 'score').<br>
Range: 0-255<br>
Default: 255

`-f|--format <FQ,BED,SAM,BG,FDENS,RDENS,FDIST,RDIST>`<br>
specifies output file formats.<br>
Value `FQ` forces to output the sequence. 
In paired-end mode two [FQ](https://en.wikipedia.org/wiki/FASTQ_format) files are generated, with suffixes ‘_1’ and ‘_2’.<br>
Unsorted [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) files represent the immediate (direct) precise alignment. 
They are useful for the peak-callers benchmarking because make mapping stage unnecessary. 
`SAM` format is indispensable for aligner validation as a control test.<br>
Ordered `BG` file represents coverage in [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format with the WIG extension. 
In contrast to the utilities that restore the coverage by extending the read to the average fragment length 
(such as [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html), 
[deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), 
[peakranger wigpe](http://ranger.sourceforge.net/manual1.18.html)), **isChIP** produces an actual coverage. 
The difference can be observed in the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/formal-actual_coverage_legend.png).<br>
It is also possible to generate one BedGraph file per strand (`-x|--strand` option).<br>
`FDENS,RDENS` represent the densities of fragments and reads, respectively, in WIG format with span = 1. 
Unlike coverage, each fragment/read is represented by one point. 
For fragments this point is the center, for reads it is the 5’ position. This view can be useful in some cases, such as ChIP-exo. 
Files have the *fdens.wig*/*rdens.wig* extension.<br>
`FDIST/RDIST` is a conditional format to control the output fragment/read frequency distribution de facto. 
This is a plain text file with an extension *fdens*/*rdens* (*fdens.txt*/*fdens.txt* in Windows), 
representing the obtained distribution as a list of pairs \<fragment length\>\<number of repetitions\>. 
Such a presentation can be visualized as a graph, for example, in Excel.<br>
In addition, the header of the file displays the initial distribution parameters, and parameters restored by points de facto.<br>
Any formats can be set in any order.<br>
Default: `FQ`

`-C|--control`<br>
forces to generate control sequence (‘input’) simultaneously with test within the test mode. 
The resulting control is fully complied with the test in terms of the average specific background density.<br>
The control’s number of 'cells' and sample are adjusted to minimize run-time and are printed on screen in `PAR` and `DBG` verbose level.<br>
The control file has the same name as the test, with the suffix '_input'.<br>
In *control* mode is ignored.

`-x|--strand`<br>
forces to generate two additional WIG files, one per each strand.<br> 
It only matters in `SE` sequencing mode and activated WIG output format.

`-o|--out <file>`<br>
specifies output files location. If option’s value is a file, it is used as a common file name. 
If file has an extension, it is discarded. If value is a directory, the default file name is used.<br>
Default: *test* mode: **mTest.\***, *control* mode: **mInput.\***

`-z|--gzip`<br>
forces to compress output data files (except fragment/read size distribution)

`-V|--verbose <SL|RES|RT|PAR|DBG>`<br>
sets verbose level:<br>
`SL`	silent mode: no output except critical messages<br>
`RES`	print processing total info<br>
`RT`	print processing info for each chromosome and total<br>
`PAR`	print actual parameters value (including number of features in template) and `RT` level<br>
`DBG`	print `PAR` level with additional info, including *template* features statistics 
and info about possible omitted features.<br>
See [Output info](#output-info).<br>
Default: `PAR`


### Output info
An example of the displayed output:

```
  chrom      reads sample   ampl  r/kbp      reads   sample  ampl  r/kbp   gaps  g_excl  mm:ss
───────────────────────────────────────────────────────────────────────────────────────────────
t chr 1:  FG: 1334  100%    83.8   55   BG: 350156     1%    72.3   3.57   2.9%  1.5%    00:01
c chr 1:                    81.0            356366           69.1   3.58                 00:00
...
───────────────────────────────────────────────────────────────────────────────────────────────
t total:  FG: 93156 100%    82.1  55    BG: 944802      1%   70.9   3.54   3.6%  2.4%
total recorded reads:: test: 1037958, control: 993378
```

* `t, c` – test, control (displayed only with `-C|--control` option)
* `FG, BG` – foreground (peak), background data (in test mode, in control mode just one background)
* `reads` – number of recorded reads
* `sample` – actual sample: the number of recorded reads relative to the number of selected ones, in percent. 
Being ‘selected’ means the reads derived from fragments have passed through the size selection filter. 
Note that 'BG: sample' shows the absolute background level, not the relative one as specified by `-G|--ground` option. 
The actual sample may be reduced because of the next reasons: user-defined and/or scalable sample, different template feature’s score (for the foreground reads), exceeding the maximum number of ambiguous code N in read. 
Amplified reads are not included in the sample statistics;
* `ampl` – actual amplification coefficient, calculated as \<amplicons number\>/\<fragments number\>. 
The coefficient takes into account both amplification techniques; displayed only if amplification is set.
* `r/kbp` – actual read density, in read per 1000 base pair;
* `gaps` – total fraction of reference undefined regions (gaps); displayed only in `PAR` and `DBG` verbose mode;
* `g_excl` – fraction of at each end excluded from treatment; displayed only in the `PAR` and `DBG` verbose modes;
* `mm:ss` – wall time (only if `–t|--time` option is set)

In the ‘total’ line, the total number of reads is indicated, as well as the average values of the density and the relative number of gaps and excluded from the treatment gaps for the entire genome.

### Limitation
The total number of chromosomes is restricted by 255.<br>
With the maximal permitted parameter values of lognormal distribution, mean = 9 and sd = 1.0, 
the maximal generated fragment length is about 1200 kbp, the average fragment length is about 13 kbp.

## Model: brief description
The real protocol of ChIP-seq is simulated by repeating the basic cycle, performed for each treated chromosome (twice for autosomes). 
The basic cycle corresponds to one nominal cell and consists of the next stages:
* **'shearing of DNA'**: cutting the reference genome in fragments with size distribution fitted to lognormal distribution with given *mean* and *sigma*. 
See [Fragment distribution and size selection](#fragments-distribution-and-size-selection);
* **'ChIP'** (in *test* mode): extraction of the fragments overlapping with the *template* binding events;
* **'loss fragments in ‘Library construction'**: sampling of selected fragments according to user-defined sample (`-G|--ground` option), 
*template* sample (defined as *template* features score) and automatically adjusted scaling sample (`--bg-all` and `--rd-lim` options);
* **'amplification'** the selected fragments if required:<br>
**MDA**:<br>
A. splitting of each fragment into two random amplicons, copying of the fragment and both amplicons to the output;<br>
B. applying step A to each amplicon until it is longer than the length of read.<br>
The amplification model is simplified due to non-essential details: 
the absence of primer annealing and the fragment debranching into 2 amplicons.<br>
**PCR**:<br>
copying each fragment 2^N times, where N denotes a given number of PCR cycles.<br>
In reality the multiplication coefficient is below 2 in each cycle, it also can be lower for some particular amplicons, e.g. CG-reach. 
These details are excluded from the model as insignificant to assess the effect of amplification on the final generalized result.<br>
* **'size selection'**: selection of fragments fitted to desirable size. See [Fragment distribution and size selection](#fragments-distribution-and-size-selection);
* **'sequencing'**: cutting the 5’end of the fragment of desirable length, or the 3’end (by random choice) in SE mode, 
or both ends in PE mode, and complementary reversing the 3’end read;
* recording reads in an appropriate formats.

The input parameters of simulation process (except number of nominal cells) 
are adjusted correspondingly to those in real ChIP-Seq experiment.

The model in its basic features was developed by [Dr T. Subkhankulova](https://www.linkedin.com/in/tatiana-subkhankulova-0876a240), ICL.

## Fragment distribution and size selection
The lognormal distribution of fragments by shearing chromatin based on sonication is confirmed by many researches. 
The typical result is presented here:<br> 
![Frag Distrib Paper](https://github.com/fnaumenko/isChIP/blob/master/pict/Size_distr_analysis_lbl1.png)<br>
In practice, the distribution parameters can vary widely. 
Examples of original and recovered distributions of experimental datasets from NCBI database are shown 
in the ![Frag Distributions](https://github.com/fnaumenko/bioStat/tree/master/pict/FragPE_distrs.png).

Fragment size selection can be performed in different techniques, e.g. by using magnetic beads or by manual cutting of the gel strip. 
Nevertheless, it is safe to assume the general normal character of size selection:<br>
![Size Selection](https://github.com/fnaumenko/isChIP/blob/master/pict/Mag-Bind_label.png)<br>
This is also confirmed by the real fragment frequency distributions in the 
![ Frag Distributions](https://github.com/fnaumenko/bioStat/tree/master/pict/FragPE_distrs.png). 
In particular, experiments experiments in cases 9-13  were clearly carried out using some size selection technique.<br>
On this basis, the default values of the lognormal *sigma* and *mean* are selected in **isChIP** so 
as to provide the most frequently observed distribution with Mean of 200.<br>
The size selection in **isChIP** is carried out according to the normal distribution. 
By default, its *mean* is automatically adjusted so that it coincides with the Mean of an initial lognormal distribution. 
This provides the least computational loss: 
![Model Distribution](https://github.com/fnaumenko/isChIP/blob/master/pict/isChIP_fragDistr2_lbl.png)<br>
When the user sets its own size selection parameters, one should bear in mind a decrease in output, 
which is the greater, the further the size selection mean is from the lognormal Mean.<br>

To facilitate the adjustment of distribution parameters use of a specialized Windows utility 
[**RandomTest**](https://github.com/fnaumenko/RandomTest-Win) is recommended. 
It visualizes the initial lognormal as well as the final distribution of fragments after size selection, 
and allows you to quickly fit the parameters for the desired distribution.

## Example of single cell simulation
There are no reasons to prohibit the attempt to apply the formal ChIP-seq protocol to the single cell simulation.<br>
Here are the *in silico* experiments on sequencing of formal TFBS using MDA, PCR, and both techniques sequentially, as compared to conventional sequencing.

*Template* consisted of 1000 conditional TFBS with a length of 10.<br>
The first two tests are sequencing 100 nominal cells, which corresponds to 100,000 real cells with a total loss of 99%. 
At a given background level of 1%, such experimental data can be considered satisfactory (test 1 with size selection) 
and good (test 2 without size selection). These tests are given for comparison, as a reference.

At the moment, MDA is usually applied to superlong fragments, tens of kilobase pairs. 
They indeed provide a fairly massive output, however, the specific density remains low, and the accuracy of narrow binding site locating is far from desired. 
Test 4 is generated on the basis of a distribution with a mean fragment length of only 1300 bp, 
which is quite moderate, but it already demonstrates unsatisfactory positioning.<br>
To maintain this issue, contrary to common practice tests were performed on short fragments.<br>

Size selection drastically reduces MDA output, therefore, all tests, except for 1, 
were performed in the absence of size selection.<br>
For successful MDA output, not only the initial fragment length is critical, but also the length of read, which limits further displacement reaction. 
It seems obvious that short amplicons should make the greatest contribution. 
To check this assumption, tests 7-8 were performed on short reads.<br>
Indeed, test 8 can be considered quite successful.

| test | cells | frag distr | read len | MDA | PCR cycles | rel peak dens |
|:---:|:---:|:---:|:---:|:---:|:---:|---:|
| **1** | 100 | A | 50 | | | **33** |
| **2** | 100 | B | 50 | | | **100** |
| **3** | 1 | B | 50 | | 5 | **32** |
| **4** | 1 | C | 50 | yes | | **10** |
| **5** | 1 | B | 50 | yes | | **11** |
| **6** | 1 | B | 50 | yes | 3 | **85** |
| **7** | 1 | B | 25 | yes | | **21** |
| **8** | 1 | B | 25 | yes | 3 | **168** |

**fragment distributions** are marked according to the figure at the bottom of the page.<br>
**rel peak dens** means **relative in-peak density**. 
The peak density was calculated as the average at a distance of +/- the average fragment length from the TFBS boundaries. 
The peak density of test 2 is taken as 100.

Actual coverages in the figure are obtained directly by **isChIP** with the option `–F|--format WIG`.<br>
Peak amplitudes are displayed proportional:

![Coverages](https://github.com/fnaumenko/isChIP/blob/master/pict/test_ampl_col.png) 

Coverages in strands (read - positive, blue - negative):

![Coverages](https://github.com/fnaumenko/isChIP/blob/master/pict/test_ampl_strand_col.png) 


Fragments distributions:

![Distributions](https://github.com/fnaumenko/isChIP/blob/master/pict/distr_small.png "distributions")

Synopsis:
```
test1:  isChIP –g $G –c 19 -f wig –o test1 –tx –n 100 –S :30 tc19_10x1000-1e4.bed
test2:  isChIP –g $G –c 19 -f wig –o test2 –tx –n 100 tc19_10x1000-1e4.bed
test3:  isChIP –g $G –c 19 -f wig –o test3 –txD -a 5 tc19_10x1000-1e4.bed
test4:  isChIP –g $G –c 19 -f wig –o test4 –txD --ln 7.06: tc19_10x1000-1e4.bed
test5:  isChIP –g $G –c 19 -f wig –o test5 –txD tc19_10x1000-1e4.bed
test6:  isChIP –g $G –c 19 -f wig –o test6 –txD -a 3 tc19_10x1000-1e4.bed
test7:  isChIP –g $G –c 19 -f wig –o test7 –txD -r 25 tc19_10x1000-1e4.bed
test8:  isChIP –g $G –c 19 -f wig –o test8 –txD -r 25 -a 3 tc19_10x1000-1e4.bed
```
where `$G` is initialized by reference genome path,<br>
`tc19_10x1000-1e4.bed` contains 1000 TFBS for chromosome 19 with a length of 10, evenly distributed in the range of 10,000,000 - 20,000,000
<br>

---
<br><br>
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write on fedor.naumenko@gmail.com
