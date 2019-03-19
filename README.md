# isChIP
**I**n **S**ilico **ChIP**-seq is a fast realistic ChIP-seq simulator.

The modelling of the chromatin immunoprecipitaion followed by next generation sequencing process is based on presumably on Illumina protocol. 
However, the **isChIP**’s options can be adjusted to the other techniques such as  Ion Torrent, etc.<br>
For more information about model see [Model: brief description](#model-brief-description).

### Performance
On 2.5 GHz RAID HPC by default values of ground samples, in one-thread mode, within 1 minute **isChIP** records:<br>
in *test* mode up to 1,000,000 reads (16,500 read/sec);<br>
in *control* mode up to 60,000,000 reads (1,000,000 read/sec).<br>
The required memory is linearly proportional to the number of threads. For one thread, it does not exceed 300 Mb.

### Navigation:
  [Installation](#installation)<br>
  [Usage](#usage)<br>
  [Details](#details)<br>
  [Model: brief description](#model-brief-description)<br>
  [Fragments distribution and size selection](#fragments-distribution-and-size-selection)

## Installation
### Executable file

[Download Linux version](https://github.com/fnaumenko/isChIP/releases/download/1.0/isChIP-Linux-x64.gz)<br>
[Download Windows version](https://github.com/fnaumenko/isChIP/releases/download/1.0/isChIP-Windows-x64.zip)

Alterative for Linux: type in the desired directory:<br>
```wget -O isChIP.gz https://github.com/fnaumenko/isChIP/releases/download/1.0/isChIP-Linux-x64.gz```<br>
```gzip -d isChIP.gz```<br>
```chmod +x isChIP```

### Compiling in Linux
Required libraries:<br>
g++<br>
pthread<br>
zlib (optionally)

Go to the desired directory and type:<br>
```wget -O isChIP.zip https://github.com/fnaumenko/isChIP/archive/1.0.zip```<br>
```unzip isChIP.zip```<br>
```cd isChIP-1.0```<br>
```make```

If **zlib** is not installed on your system, the linker will display a message.<br>
In that case you can compile the program without the ability to work with .gz files. 
For this open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try ```make``` again.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

### Prepare reference genome
Download the required reference genome from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath)/<br>
For example, to download mouse library **mm9**:

**Linux**<br>
type commands:<br>
```mkdir mm9```<br>
```cd mm9```<br>
```rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/ ./```<br>
The alternative way:<br>
```wget -r ftp://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/```<br>
```mv hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes mm9```<br>
```rm -r hgdownload.soe.ucsc.edu```<br>
If you do not want to keep unused modifications, type<br>
```rm mm9/*_random*```

**Windows**<br>
copy and paste the string *ftp://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/* into Windows browser address bar, 
then copy *.fa.gz files to your local directory.<br>
Alternatively use FTP client, e.g. [FileZilla](https://filezilla-project.org/).

## Usage
```
isChIP [options] -g|--gen <name> [<template>]
	<template> - bed file whose features specify binding sites (BS)
```

### Synopsis
```isChIP -g ref_genome –n 5```<br>
Generates 'input' sequences in FastQ with read length of 50 and average read density of about 9 read/kbs, 
comparable with what is experimentally observed.

```isChIP -g mm9_dir –tz –n 300 –r 36 –f fq,sam templ.bed```<br>
Generates test sequences in zipped FastQ and SAM, with read length of 36, timing, and average foreground/background read density 
comparable with what is experimentally observed in [Series GSE56098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56098).

### Help
```
Processing:
  -g|--gen <name>       reference genome library or single nucleotide sequence. Required
  -n|--cells <int>      number of cells [1]
  -a|--ampl-c <int>     number of PCR cycles [0]
  -b|--bg <float>       number of selected fragments outside BSs (background),
                        in percent of foreground. For the test mode only [1]
  --fg <float>          in test mode the number of selected fragments within BSs, in percent;
                        in control mode the number of selected fragments, in percent [100]
  --bg-all <OFF|ON>     turn on/off generation background for all chromosomes.
                        For the test mode only [ON]
  -c|--chr <name>       generate output for the specified chromosome only
  -N|--let-N            include the regions filled with an ambiguous reference characters 'N'
                        on the beginning and on the end of chromosome
  -m|--smode <SE|PE>    sequencing mode: SE - single end, PE - paired end [SE]
  -u|--ts-uni           uniform template score. For the test mode only
  -p|--threads <int>    number of threads [1]
  -seed <int>           fix random emission with given seed, or 0 if do not fix [0]
Fragment distribution:
  --fr-mean <float>     mean of fragment lognormal distribution [5.46]
  --fr-sd <float>       standard deviation of fragment lognormal distribution [0.4]
  --ss-mean <int>       mean of size selection normal distribution [auto]
  --ss-sd <int>         standard deviation of size selection normal distribution [30]
  --ss <OFF|ON>         turn off/on fragment's size selection [ON]
Reads:
  -r|--rd-len <int>     length of output read [50]
  --rd-name <NONE|NUMB|POS>     info added to read's name in output files:
                        NONE - nothing
                        NUMB - read`s unique number across genome
                        POS  - read`s actual start position [NONE]
  --rd-Nlim <int>       maximum permitted number of ambiguous characters 'N' in read [OFF]
  --rd-lim <long>       maximum permitted number of total recorded reads [2e+08]
  --rd-ql <char>        uniform quality value for the sequence  [~]
  --rd-ql-patt <name>   quality values pattern for the sequence 
  --rd-map-ql <int>     read mapping quality for SAM and BED output [255]
Output:
  -f|--format <FQ,BED,SAM,WIG,FREQ>     format of output data, in any order  [FQ]
  -C|--control          generate control simultaneously with test
  -S|--strand           generate two additional wig files, each one per strand  
  -s|--sep              print number of reads with '1000' separator
  -o|--out <name>       location of output files or existing directory
                        [Test mode: mTest.*, Control mode: mInput.*]
  -z|--gzip             compress the output
Other:
  -t|--time             print run time
  -V|--verbose <SL|RES|RT|PAR|DBG>      set verbose level:
                        SL -    silent mode (show critical messages only)
                        RES -   show result summary
                        RT  -   show run-time information
                        PAR -   show actual parameters
                        DBG -   show debug messages [RT]
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

## Details

### Mode

**isChIP** generates output in one of two modes:<br>
*test* – simulation of site of interest sequencing; the output is test sequences/alignment<br>
*control* – simulation of control production; the output is 'input' sequences/alignment.<br>
They are distinguished only by involvement or elimination of sites of interest in a process. 
The sites are represented by the set of features stated in a single optional parameter – BED file, called *template*.

### Template
*Template* is a file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format whose features correspond to binding events. 
Each data line contains 3 (minimum required) or 5 and more fields. 
First 3 fields define the binding site. 
If there are no more fields in the line, all generating significant enriched regions have the same density. 
If 4th field ('name') and 5th field ('score') are in a line, **isChIP** uses the feature’s score to restrict appropriate enriched density. 
The unit does not matter, only score ratios are meaningful. 
This peculiarity allows the use of BED files generated by peak callers, though different peak callers have different score units. 
At the same time different stated scores can be ignored by applying  the option ```-u|--ts-uni```.<br>
*Template* does not have to be sorted, but the features must be grouped by chromosomes. 
This means that features belonging to the same chromosome must be arranged sequentially, in a single group. 
The simplest way to ensure this is to pre-sort the file.

### Options description

*Note:*<br>
Enumerable option values are case insensitive.<br>
Compressed input files in gzip format (.gz) are acceptable.

```-g|--gen <name>```<br>
Reference genome library, or single nucleotide sequence.<br>
Genome library is a directory contained nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a directory, first **isChIP** searches .fa files in it.
If there are not such files, or the file corresponded to chromosome specified by option ```–c|--chr``` option is absent, isChIP searches .fa.gz files.

There are 3 ways to define desirable reference sequences.<br>
To generation output for the whole genome, ```name``` should be a path, 
and all chromosome files should be in a uniform compressed condition.<br>
To restrict output to some chromosomes, ```name``` should be a path, 
and only desirable chromosome files should be decompressed or even present at all.<br>
To generate output for a single chromosome, ```name``` can be a path and option ```–c|--chr``` 
should be set, or ```name``` should indicate an appropriate reference file.<br>
In *test* mode the target references are determined by *template*. 
See also ```--bg-all``` and ```–c|--chr ``` options.

**isChIP** omits  'random' contigs and haplotype sequences.
One can obtain a genome library in  UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, 
e.g. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```-n|--cells <int>```<br>
The number of conditional cells. Each conditional cell corresponds to one passage through the reference chromosome (dual for the autosomes).<br>
The convention is that this number assigns a lossless output by default (zero value of ```--fg``` option). 
If you want to deal with a real number of cells, you should assign a real level of overall loss. 
The key is that this proportionally increases the runtime. 
For example, assigning 100 conditional cells with zero losses leads to the generation of data equivalent to 100,000 ‘real’ cells with a ‘real’ total loss of 99.9%, 
but in a runtime about 1000 times less (e.g. 10 seconds instead of almost 3 hours per chromosome).

Since the total losses are difficult to estimate in advance, the inverse task, namely, 
determining the overall loss in a real experiment post-factum, looks much more practical. 
To do this, it is enough to select a model output with a peak density corresponding to the peak density in the actually observed data through several fast iterations with different number of conditional cells and zero foreground level. 
Then the overall loss is calculated by the formula (1-n/N)*100%, where n is the number of conditional cells, N is the number of real cells in the experiment.

Besides the number of cells, the program’s output heavily depends on parameters of fragment distribution and size selection 
(see [Fragments distribution and size selection](#fragments-distribution-and-size-selection)). 
In general, by default distribution, a number of conditional cells from 3 to 10 provides an output read’s mean density comparable with real ‘inputs’, 
and a value from 80 to 400 leads to data comparable with real tests in term of density.<br>
Range: 1-2e6<br>
Default: 1

```-a|--ampl-c <int>```<br>
The number of PCR amplification cycles. 
See [Model: brief description](#model-brief-description).<br>
Value 0 means non amplification.<br>
Range: 0-100<br>
Default: 0 

```-b|--bg <float>```<br>
The background level: in *test* the number of selected fragments which are not intersected with the *template* binding events, 
as a percentage of the foreground (see ```--fg``` option).<br>
In *control* mode it is ignored.<br>
In practice a level of 1-3% corresponds with a good experimental data set, 
while a level of more than 6% leads to generate data with a rather low signal-to-noise ratio.<br>
Range: 0-100<br>
Default: 1

```--fg <float>```<br>
The foreground level: in *test* mode the number of selected fragments which are intersected with the *template* binding events, 
as a percentage of the total number of generated fragments.<br>
In *control* mode the number of selected fragments, as a percentage of the total number of generated fragments.<br>
This option is designed to simulate overall loss (see ```-n|--cells``` option).<br>
Range: 0-100<br>
Default: 100

```--bg-all <OFF|ON>```<br>
In *test* mode turn off/on background generation for all chromosomes, whether or not they are presented in *template*.<br>
Mapping one chromosome to the whole reference genome leads to the appearance of short local lacks of alignment, 
‘gaps’, corresponding to low mappability regions. This is true for all aligners. 
Generation background for all chromosomes eliminates this issue.<br>
Note that this process is time consuming, and is justified only in case of using subsequent alignment, 
for example, for comparison with experimental FastQ, or to validate aligners. 
When validating peak detectors on one or some chromosomes, it is more advantageous to use a direct SAM output format (see ```-f|--format```).<br>
The established option ```-c|--chr``` makes this option negligible.<br>
In *control* mode is ignored.<br>
Default: ```ON```

```-c|--chr <name>```<br>
Generate output for the specified chromosome only. 
```name``` means chromosome identifier, i.e. number or  character, for instance ```10```, ```X```. 
This creates the same effect as referencing to the chromosome file instead of directory.<br>
This is a strong option, which forces to ignore all other chromosomes from reference genome and *template*, 
and abolishes the impact of option ```--bg-all```.

```-N|--let-N```<br>
Forces the scanning process to include region consisting of ambiguous characters 'N' at the beginning and at the end of a reference chromosome. 
It makes no difference in data after alignment, but slightly increases the output volume and the runtime.

```-m|--smode <SE|PE>```<br>
Generate reads according to stated sequencing mode: ```SE``` – single end, ```PE``` – paired end 
(see [Model: brief description](#model-brief-description)).<br>
Default: ```SE```

```-u|--ts-uni```<br>
In *test* mode forces to ignore features scores established in *template*. 
All binding sites will be simulated with maximum score.<br>
In other modes is ignored.

```-p|--threads <int>```<br>
Number of threads. The workflow is separated between chromosomes, so the actual number of threads 
can be reduced (if the number of actual treated chromosomes is less then assigned value). 
The actual threads number is displayed in ```PAR``` and ```DBG``` verbose mode (see ```-V|--verbose``` option).<br>
Range: 1-20<br>
Default: 0

```--seed <int>```<br>
Fix random numbers emission to get repetitive results. 
The actual seed equals the option value increased by a certain factor to provides a noticeable difference in the of random number generation option values that differ by 1.<br>
Value 0 means non-recurring random generation.<br>
Range: 0-100<br>
Default: 0

```--fr-mean <float>```<br>
Mean of fragment’s lognormal distribution. 
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection).<br>
Range: 2.0-9.0<br>
Default: 5.46

```--fr-sd <float>```<br>
Standard deviation of fragment’s lognormal distribution. 
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection).<br>
Range: 0.1-1.0<br>
Default: 0.4

```--ss-mean <int>```<br>
Mean of size selection normal distribution. 
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection).<br>
By default this value is calculated as the mode of the fragment’s lognormal distribution. 
In other words, option's value is equal to the most frequent value of the lognormal distribution, which for given defaults is 200.<br>
Range: 50-1000<br>
Default: ```auto```

```--ss-sd <int>```<br>
Standard deviation of fragment's size selection normal distribution. 
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection).<br>
Range: 1-200<br>
Default: 30

```-r|--rd-len <int>```<br>
Length of output read.<br>
Note, that **isChIP** does not perform a preliminary check of the read length relative to the mode of fragments distribution or the mean of size selection, 
and therefore does not issue warnings. 
If the reads are too long, their amount at the output is reduced up to 0.<br>
Range: 20-500<br>
Default: 50

```--rd-name <NONE|NUMB|POS>```<br>
Adds stated info to the name of read in output files:<br>
```NONE``` – the read’s name consists only of the application name;<br>
```NUMB``` – actual chromosome and read’s unique number across genome are added;<br>
```POS```  – actual chromosome and read’s actual start position are added.<br>
Position info is useful for verifying aligners. Unique number is required for the paired-end reads aligning.<br>
Additional info slightly increases the runtime and size of output files.<br>
Default: ```NONE(SE)|NUMB(PE)```

```--rd-Nlim <int>```<br>
Maximum permitted number of ambiguous reference characters ('N') in read.<br>
This option is designed primarily for validation of the aligners.<br>
In the real read sequences ambiguous characters could cause a failure in sequencing as well as undefined regions in the cut fragments, 
while simulated reads are fully defined by reference library. 
Consequently the ambiguous characters in the library are translated into reads. 
Different aligners considered undefined characters in references differently: 
some of them as invalid, and some of them as an overlapping case.<br>
Range: 0 - <```-r|--rd-len``` value><br>
Default: ```-r|--rd-len``` value (what is the same as no checkup)

```--rd-lim <long>```<br>
Maximum number of total written reads. The value emulates sequencer’s limit. 
It restricts the number of recorded reads for each chromosome proportionally.<br>
Range: 1e5-1e19<br>
Default: 2e8.

```--rd-ql <char>```<br>
Uniform quality value for the sequence (read) in FQ and SAM output.<br>
Range: '!'-'\~'<br>
Default: '\~' (decimal 126)

```--rd-ql-patt <name>```<br>
Quality values pattern for the sequence (read) in FQ and SAM output. 
```name``` is a plain text file containing at least one line encodes the quality values for the sequence 
as described in [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) specification.<br>
If the length of line is less then read length, the rest of pattern is filed by value defined by ```--rd-ql``` option.<br>
If the length of line is more then read length, the rest of line is ignored.<br>
The lines starting with the character '#' are ignored.<br>
The second and all the following encoding lines are ignored as well.

```-rd-map-ql <int>```<br>
Read mapping quality in SAM and BED output (in the last case it is called 'score').<br>
Range: 0-255<br>
Default: 255

```-f|--format <FQ,BED,SAM,WIG,FREQ>```<br>
Output file formats.<br>
Value ```FQ``` forces to output the sequence. 
In paired-end mode two [FQ](https://en.wikipedia.org/wiki/FASTQ_format) files are generated, with suffixes ‘_1’ and ‘_2’.<br>
Unsorted [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) files represent the immediate (direct) precise alignment. 
They are useful for the peak-callers benchmarking because make mapping stage unnecessary. 
```SAM``` format is indispensable for aligner validation as a control test.<br>
Ordered ```WIG``` file represents coverage in [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) format. 
In contrast to the utilities that restore the coverage by extending the read to the average fragment length 
(such as [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html), 
[deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), 
[peakranger wigpe](http://ranger.sourceforge.net/manual1.18.html)), **isChIP** produces an actual coverage. 
It is also possible to generate one wig file per strand (```-S|--strand``` option). 
You can see the difference in the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/formal-actual_coverage_label.png).<br>
```FREQ``` is a conditional format to control the output fragment frequency distribution de facto. 
This is a plain text file with an extension .freq (.freq.txt in Windows), 
representing the obtained distribution as a list of pairs \<fragment length\>\<number of repetitions\>. 
Such a presentation can be visualized as a graph, for example, in Excel.<br>
Any formats can be set following in any order.<br>
Default: ```FQ```

```-C|--control```<br>
In *test* mode forces to generate control sequence ('input') simultaneously with test. 
The resulting control is fully complied with the test according to the criterion of the average specific background density.<br>
The control’s number of 'cells' and sample are adjusted to minimize run-time and are printed on screen in ```PAR``` and ```DBG``` verbose level.<br>
The control file has the same name as the test, with the suffix '_input'.<br>
In *control* mode is ignored.

```-S|--strand```<br>
In addition to the WIG format forces to generate two additional wig files, each one per strand. 
Is valid only in ```SE``` sequencing mode.

```-o|--out <file>```<br>
Output files location. If option’s value is a file, it is used as a common file name. 
If file has an extension, it is discarded. If value is a directory, the default file name is used.<br>
Default: *test* mode: **mTest.\***, *control* mode: **mInput.\***

```-V|--verbose <SL|RES|RT|PAR|DBG>```<br>
Set verbose level:<br>
```SL```	silent mode: no output except critical messages<br>
```RES```	print processing total info<br>
```RT```	print processing info for each chromosome and total<br>
```PAR```	print actual parameters value (including number of features in template) and ```RT``` level<br>
```DBG```	print ```PAR``` level with additional info, including *template* features statistics.<br>
See [Output info](#output-info).<br>
Default: ```RT```


### Output info
The next information is displaying for each chromosome and total (field values are given as an example):

```
  chrom      reads sample  r/kbp        reads   sample  r/kbp     N   N_excl  mm:ss
───────────────────────────────────────────────────────────────────────────────────
t chr 1:  FG: 1334  100%      55   BG: 350155      1%    3.57   5.2%  4.9%    00:01
c chr 1:                               356366            3.58   5.2%  4.9%    00:00
...
───────────────────────────────────────────────────────────────────────────────────
t total:  FG: 93156 100%      55   BG: 944801      1%    3.54
total recorded reads:: test: 474931, control: 993378
```

* ```t, c``` – test, control (only with ```-C|--control``` option)
* ```FG, BG``` – foreground (peak), background data (in test mode, in control mode just one background)
* ```reads``` – number of recorded reads
* ```sample``` – actual sample: the number of recorded reads relative to the number of selected ones, in percent. Selected means the reads derived from fragments passed through the size selection filter. 
The actual sample may be reduced because of the next reasons: user-defined and/or scalable sample, different template feature’s score (for the foreground reads), exceeding the maximum number of ambiguous reference characters ‘N’ in read. 
Amplified reads are not included in the sample statistics;
* ```r/kbp``` – actual read density, in read per 1000 base pair;
* ```N``` – total fraction of ambiguous reference characters 'N' in chromosome (only in ```PAR``` and ```DBG``` verbose mode);
* ```N_excl``` – fraction of excluded from treatment ambiguous reference characters 'N' in chromosome (only in ```PAR``` and ```DBG``` verbose mode);
* ```mm:ss``` – wall time (only if ```–t|--time``` option is set)

### Limitation
The total number of chromosomes is restricted by 255.<br>
By permitted maximum values of lognormal mean = 9 and sd = 1.0 the maximum generated fragment length is about 1200 kb, 
average length is about 13 kbp.

## Model: brief description
The real protocol of ChIP-seq is simulated by repeating the basic cycle, performed for each treated chromosome (twice for autosomes). 
The basic cycle consists of the next phases:
* **'shearing of DNA'**: cutting the reference genome in fragments with size distribution fitted to lognormal distribution with given *mean* and *sigma*;
* **'ChIP'** (in *test* mode): extraction of the fragments overlapping with the *template* binding events;
* **'loss fragments in ‘Library construction’'**: sampling of selected fragments according to user-defined sample (```--fg``` and ```-b|--bg``` options), 
*template* sample (defined as *template* features score) and automatically adjusted scaling sample (```--bg-all``` and ```--rd-lim``` options);
* **'amplification'** the selected fragments if required:<br>
copying each fragment 2^N times, where N denotes a given number of PCR cycles;
* **'size selection'**: selection of fragments fitted to desirable size. See [Fragments distribution and size selection](#fragments-distribution-and-size-selection);
* **'sequencing'**: cutting the 5’end of the fragment of desirable length, or the 3’end (by random choice) in SE mode, 
or both ends in PE mode, and complementary reversing the 3’end read;
* recording reads in an appropriate formats.

The input parameters of simulation process (except number of conditional cells) 
are adjusted correspondingly to those in real ChIP-Seq experiment.

The model was developed by [Dr. Tatiana Subkhankulova](https://www.linkedin.com/in/tatiana-subkhankulova-0876a240), University of Cambridge.

## Fragments distribution and size selection
The lognormal distribution of fragments by shearing chromatin based on sonication confirmed by many researches.<br>
In practice, the distribution parameters can vary widely:<br>
![Real Distributions](https://github.com/fnaumenko/isChIP/blob/master/pict/fragDistr_ChIP-seq_label_small.png) 

Fragment size selection can be performed in different techniques, e.g. by using magnetic beads or by manual cutting of the gel strip. 
Nevertheless, it is safe to assume the general normal character of size selection:<br>
![Size Selection](https://github.com/fnaumenko/isChIP/blob/master/pict/Mag-Bind_label.png)<br>
This is also confirmed by the real fragment frequency distributions in the first figure. 
In particular, experiments SRR408580 (in green) and SRR965509 (in yellow) were clearly carried out using some size selection technique.<br>

On this basis, the default values of the lognormal *sigma* and *mean* are selected in **isChIP** so 
as to provide the most frequently observed distribution with mode of 200.<br>
The size selection in **isChIP** is carried out according to the normal distribution. 
By default, its *mean* is automatically adjusted so that it coincides with the mode of an initial lognormal distribution. 
This provides the least computational loss. Of course the size selection mean can be set by the user. 
But in this case, one should bear in mind the decrease in output, the more, the further the size selection mean is from the lognormal mode (yellow graph):<br>
![Model Distribution](https://github.com/fnaumenko/isChIP/blob/master/pict/isChIP_fragDistr_label_.png)

To facilitate the adjustment of distribution parameters use a specialized Windows utility [**RandomTest**](https://github.com/fnaumenko/RandomTest-Win). 
It visualizes the initial lognormal as well as the final distribution of fragments after size selection, 
and allows you to quickly fit the parameters for the desired distribution.



##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
