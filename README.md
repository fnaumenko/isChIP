# isChIP
**I**n **S**ilico **ChIP**-seq is a fast realistic ChIP-seq simulator.

The modelling of the chromatin immunoprecipitaion followed by next generation sequencing process is based on presumably on Illumina protocol. 
However, the **isChIP**’s options can be adjusted to the other techniques such as  Ion Torrent, etc.

Simulated binding events are specified by optional parameter, called [*template*](#template).<br>
For more information about model see [Model: brief description](#model-brief-description).

The program runs on the command line under Linux and Windows.

## Installation
### Executable file

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O isChIP.gz https://github.com/fnaumenko/isChIP/releases/download/1.0/isChIP-Linux-x64.gz```<br>
```gzip -d isChIP.gz```<br>
```chmod +x isChIP```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/isChIP/releases/download/1.0/isChIP-Windows-x64.zip) and unzip by any archiver, for instance [WinRar](https://www.win-rar.com/download.html?&L=0).

### Compiling in Linux
Required libraries:<br>
g++<br>
pthread<br>
zlib (optionally)

Go to the desired directory and type commands:<br>
```wget -O isChIP.zip https://github.com/fnaumenko/isChIP/archive/1.0.zip```<br>
```unzip isChIP.zip```<br>
```cd isChIP-1.0```<br>
```make```

If **zlib** is not installed on your system, a linker message will be displayed.<br>
In that case you can compile the program without the ability to work with .gz files: 
open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try ```make``` again.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

### Prepare reference genome
Download the required reference genome from UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/.<br>
For example, to download mouse library **mm9** go to the desired directory and...<br>

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
Copy and paste the string ftp://hgdownload.soe.ucsc.edu/goldenPath/mm9/chromosomes/ into Windows browser address bar, 
then copy *.fa.gz files to your local directory.<br>
The alternative way:<br>
use FTP client, f.e. [FileZilla](https://filezilla-project.org/).

## Synopsis
```isChIP -g ref_genome –n 5```<br>
Generates 'input' sequences in FastQ with read length of 50 and average read density of about 9 read/kbs, comparable with what is experimentally observed.

```isChIP -g mm9_dir –n 250 –b 4 –r 36 –f fq,sam templ.bed```<br>
Generates test sequences in FastQ and direct alignment in SAM, with read length of 36 and average foreground/background read density comparable with what is experimentally observed in [Series GSE56098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56098).<br>
```templ.bed``` consists of lines with features defining well-proven Oct4 binding motif in the centre of some peaks in the referenced experimental data:
```
chr1   9899347   9899362  ...
chr1  16120224  16120239  ...
chr1  40360157  40360172  ...
. . .
```
Please note that correlation between the whole experimental and model alignment is incorrect due to the abundance of non-simulated peaks and artefacts.

## Usage
```isChIP [options] -g|--gen <name> [template]```

### Help
```
Input:
  -g|--gen <name>       reference genome library or single nucleotide sequence. Required
Processing:
  -a|--amplif <int>     coefficient of amplification [1]
  -b|--bg-level <float> number of selected fragments outside the features,
                        in percent of foreground. For the test mode only [1]
  --fg-level <float>    in test mode the number of selected fragments within the features, in percent;
                        in control mode the number of selected fragments, in percent [100]
  -n|--cells <long>     number of cells [1]
  -c|--chr <name>       generate output for the specified chromosome only
  --bg-all <OFF|ON>     turn on/off generation background for all chromosomes. For the test mode only [ON]
  --bind-len <int>      minimum binding length. For the test mode only [1]
  --flat-len <int>      boundary flattening length. For the test mode only [0]
  --let-N               include the ambiguous reference characters (N) on the beginning
                        and on the end of chromosome
  -m|--smode <SE|PE>    sequencing mode: SE - single end, PE - paired end [SE]
  --strand-admix <OFF|ON>       turn on/off opposite strand admixture at the bound of binding site.
                        For the test mode only [OFF]
  --ts-uni              uniform template score. For the test mode only
  -p|--threads <int>    number of threads [1]
  --fix                 fix random emission to get repetitive results
  -R|--regular <int>    regular mode: write each read on starting position increased by stated shift
Fragment:
  --frag-len <int>      average size of selected fragments [200]
  --frag-dev <int>      deviation of selected fragments [20]
  --sz-sel <OFF|ON>     turn on/off fragment's size selection [ON]
  --sz-sel-sigma <int>  standard deviation of the size selection normal distribution [20]
Fragment's size distribution:
  --mean <int>          expectation of the based normal distribution [200]
  --sigma <int>         standard deviation of the based normal distribution [200]
  --ln-factor <int>     power multiplication factor in lognormal distribution [500]
  --ln-term <float>     power summand in lognormal distribution [5.1]
Reads:
  -r|--rd-len <int>     length of generated read [50]
  --rd-name <NMB|POS>   name of read in output files includes:
                        NMB - read`s unique number within chromosome
                        POS - read`s true start position [POS]
  --rd-Nlimit <int>     maximum permitted number of ambiguous characters (N) in read [--read-len]
  --rds-limit <long>    maximum permitted number of total written reads [2e+08]
  --rd-ql <char>        uniform quality value for the sequence  [~]
  --rd-ql-patt <name>   quality values pattern for the sequence 
  --rd-map-ql <int>     read mapping quality for SAM and BED output [255]
Output:
  -f|--format <FQ,BED,SAM>      format of output sequences/alignment, in any combination [FQ]
  -o|--out <name>       location of output files or existing directory
                        [Test mode: mTest.*, Control mode: mInput.*, Regular mode: mRegular.*]
  -z|--gzip             compress output files with gzip
Other:
  -t|--time             print run time
  -V|--verbose <CRIT|RES|RT|PAR|DBG>    set verbose level:
                        CRIT -  show critical messages only (silent mode)
                        RES -   show result summary
                        RT  -   show run time information
                        PAR -   show process parameters
                        DBG -   show debug messages [RT]
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

## Details

### Mode

**isChIP** generates output in one of three modes:<br>
*test* – simulation of site of interest sequencing; the output is test sequences/alignment<br>
*control* – simulation of control production; the output is 'input' sequences/alignment<br>
*regular* – simple regular cutting of reference chromosome, an auxiliary mode for special use.<br>
*Test* and *control* modes are distinguished only by involvement or elimination of sites of interest in a process. 
The sites are represented by the set of features stated in a BED file, called *template*. 
*Template* is a single optional parameter.

### Template
*Template* is a file in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format whose features correspond to binding events. 
Each data line contains 3 (minimum required) or 5 and more fields. 
First 3 fields define the binding site. 
If there are no more fields in the line, all generating significant enriched regions have the same density. 
If 4th field ('name') and 5th field ('score') are in a line, **isChIP** uses the feature’s score to restrict appropriate enriched density. 
The unit does not matter, only score ratios are meaningful. 
This peculiarity allows the use of BED files generated by peak callers, since different peak callers have different score units. 
Nevertheless it is possible to ignore stated scores by using the option ```--ts-uni```.

Compressed files in gzip format (.gz) are acceptable.

### Options description
Non-numeric option values are case insensitive.

```-g|--gen <name>```<br>
Reference genome library, or single nucleotide sequence.<br>
Reference sequence or directory contained reference sequences in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a directory, first **isChIP** searches .fa files in it.
If there are not such files in this directory, or file specified by ```–c|--chr``` option is absent, isChIP searches .fa.gz files.

There are 3 ways to define desirable reference sequences.<br>
To generation output for the whole genome, ```name``` should be a path, and all chromosome files should be in a uniform compressed condition.<br>
To restrict output to some chromosomes, ```name``` should be a path, and only desirable chromosome files should be decompressed or even present at all.<br>
To generate output for a single chromosome, ```name``` can be a path and option ```–c|--chr``` should be set, or ```name``` should indicate an appropriate reference file.<br>
In *test* mode the target references are determined by *template*. See also ```--bg-all``` and ```–c|--chr ``` options.

**isChIP** omits  'random' contigs and haplotype sequences.
One can obtain a genome library in  UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```-a|--amplify <int>```<br>
The coefficient of MDA-amplification of fragments passed through the size selection filter.<br>
The default value of 1 means non amplification.

```-b|--bg-level <float>```<br>
The background level: in *test* the number of selected fragments which are not intersected with the *template* binding events, 
as a percentage of the foreground (see ```--fg-level``` option).<br>
In *control* mode it is ignored.<br>
In practice a level of 1-5% corresponds with a good experimental data set, while a level of more than 10% would be a bad data set.<br>
Default: 1

```--fg-level <float>```<br>
The foreground level: in *test* mode the number of selected fragments which are intersected with the *template* binding events, 
as a percentage of the total number of generated fragments.<br>
In *control* mode the number of selected fragments, as a percentage of the total number of generated fragments.<br>
This option is constructed to reduce output sensitivity to number of cells, or to estimate the overall loss in the real process.<br>
Default: 100

```-n|--cells <long>```<br>
The number of 'cells'. Each 'cell' corresponds to one passage through the reference chromosome (dual for the numeric ones).<br>
It is NOT the equivalent of the number of cells in the real experiment. 
In the last case the overall loss can be as high as 99%. 
Using the default value of ```--fg-level``` option the model does not lose any related fragment.<br>
The program’s output heavily depends on correlation between established fragment size and parameters of fragment distribution and size selection. 
In general, by default a distribution of numbers from 3-10 provides an output read mean density comparable with the actual ‘inputs’, 
and a distribution of numbers from 100-500 leads to data comparable with the actual tests in term of density.<br>
Default: 1

```-c|--chr <name>```<br>
Generate output for the specified chromosome only. 
```name``` means short chromosome name, i.e. number or  character, for instance ```–c 10```, ```--chr X```. 
This creates the same effect as referencing to the chromosome file instead of directory. 
This is a strong option, which abolishes the impact of option ```--bg-all``` and all other chromosomes from *template*.

```--bg-all <OFF|ON>```<br>
In *test* mode turn on/off background generation for all chromosomes, irrespective of chromosomes included in *template*.<br>
As we discovered, mapping by any aligner one or several chromosomes to the whole reference genome leads to short local lacks of alignment, 'gaps', corresponding to low mappability regions. 
Generation background for all chromosomes eliminates this issue.<br>
Please note that including background for all chromosomes is time consuming. 
Inclusion of background for all chromosomes can be avoided by generating direct alignments in SAM or BED (see ```-f|--format```), 
because it excludes the external mapping process.<br>
Default: ```ON```

```--bind-len <int>```<br>
In *test* mode the minimum binding length. 
That is a minimum number of nucleotides that ensures the binding while fragment intersects binding site. 
If some of the features in *template* are less then given binding length, 
they will be skipped and a summary warning message will appear. 
To show the detailed warning message about each skipped feature, set verbose level to ```DBG``` (see ```--verbose```).<br>
Default: 1

```--flat-len <int>```<br>
In *test* mode, the boundary flattening length is a distance from the boundary of a binding site on which the probability of a fragment associating is increased from zero to maximum. 
As such it simulates the smoothing of enriched regions.<br>
Default: 0

```--let-N```<br>
As a rule, the first (and sometimes the last) tens or hundreds of kilobases in the reference chromosomes are meaningless. 
i.e. filled with ambiguous reference characters 'N'. 
By default **isChIP** excludes these initial and final regions from the generation.<br>
This option forces to scan the entire chromosome. It makes no difference in data after alignment, 
but increases a little a run time and a quality of random number distribution at the beginning of the process.

```--smode <SE|PE>```<br>
Generation reads according to stated sequencing mode: ```SE``` – single end, ```PE``` – paired end.<br>
Default: ```SE```

```--strand-admix <OFF|ON>```<br>
In *test* mode turn on/off opposite strand admixture at the bound of binding site.<br>
In accordance with the nature of the binding of proteins, at the boundaries of the binding site, there should be reads with the same strand only. 
However, in practice, often the experimental sequences demonstrate the presence of small admixtures of reads with the opposite strand (see ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/mix-strand.png)). 
This option allows to simulate this effect.<br>
Default: ```OFF```

```--ts-uni```<br>
In *test* mode ignores features scores if they are established in *template*. 
All features will be simulated with maximum score.

```-p|--threads <int>```<br>
Number of threads. 
The workflow is separated between chromosomes, so this option takes no effect by processing the single chromosome.

```--fix```<br>
Fix random numbers emission to get repetitive results.

```-R|--regular <int>```<br>
*Regular* mode: write each read on starting position increased by stated shift.<br>
This mode is used for specific tasks.

```--frag-len <int>```<br>
Average size of selected fragments.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 200

```--frag-dev <int>```<br>
Deviation of selected fragments.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 20

```--sz-sel <OFF|ON>```<br>
Turn on/off fragment's size selection.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: ```ON```

```--sz-sel-sigma <int>```<br>
Standard deviation of the fragment's size selection normal distribution.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 20

```--mean <int>```<br>
Expectation of the fragment's size based normal distribution.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 200

```--sigma <int>```<br>
Standard deviation of the fragment's size based normal distribution.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 200

```--ln-factor <int>```<br>
Power multiplication factor in fragment's size lognormal distribution.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 500

```--ln-term <float>```<br>
Power summand in fragment's size lognormal distribution.<br>
For more information see [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 5.1

```--rd-name <NMB|POS>```<br>
Forces to include in the name of each read its unique number (```NMB```) or its true start position (```POS```).<br>
Default: ```POS```

```--rd-Nlimit <int>```<br>
Maximum permitted number of ambiguous reference characters ('N') in output reads.<br>
In the real read sequences ambiguous characters could cause a failure in sequencing as well as undefined regions in the cut fragments, 
while simulated reads are fully defined by reference library. 
Consequently the ambiguous characters in the library are translated into reads. 
Different aligners considered undefined characters in references differently: some of them as invalid, and some of them as an overlapping case.<br>
Default: length of read (all characters could be ambiguous)

```--rds-limit <long>```<br>
Maximum number of total written reads. The value emulates sequencer’s limit.<br>
This value restricts the number of written reads for each chromosome proportionally.<br>
Default: 200 000 000. In practical simulation the default value is never achieved.

```--rd-ql <char>```<br>
Quality value for all positions in the read for the sequence (for ```FQ``` and ```SAM``` output).<br>
Default: '~' (decimal 126, maximum)

```--rd-ql-patt <name>```<br>
Set pattern of the read quality values for the sequence (for ```FQ``` and ```SAM``` output). 
```name``` is a plain text file containing at least one line encodes the quality values for the sequence exactly as described in [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) format.<br>
If the length of line is less then read length, the rest of pattern is filed by value defined by ```--rd-ql``` option.<br>
If the length of line is more then read length, the rest of line is ignored.<br>
The lines starting with the character '#' are ignored.<br>
The second and all the following encoding lines are also ignored.

```--rd-map-ql <int>```<br>
Read mapping quality for ```SAM``` and ```BED``` output (in the last case it is called 'score').<br>
Default: 255 (maximum)

```-f|--format <FQ,BED,SAM>```<br>
Output files formats. 
Value ```FQ``` forces to output the sequence. 
In paired end mode two [FQ](https://en.wikipedia.org/wiki/FASTQ_format) files are generated, with suffixes ‘_1’ and ‘_2’. 
[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [SAM](https://en.wikipedia.org/wiki/SAM_(file_format)) files contain the immediate (direct) alignment. 
Any formats can be set, but a minimum of one is allowed.<br>
Default: ```FQ```

```-o|--out <file>```<br>
Output files location. If option’s value is a file, it is used as a common file name. 
If file has an extension, it is discarded. 
If value is a directory, the default file name is used.<br>
Default: *test* mode: **mTest.\***, *control* mode: **mInput.\***, *regular* mode: **mRegular.\***

## Model: brief description
The real protocol of ChIP-seq is simulated by repeating the basic cycle. 
Each basic cycle corresponds to single cell simulation, and consists of the next phases:
* 'shearing' the chromatin, or random cutting the reference genome in fragments of correspondent length; the distribution of fragments by default corresponds to log-normal with correspondent parameters;
* 'extraction' of the fragments overlapping with the binding events;
* amplification the selected fragments if required in number of cycles;
* 'loss' of selected fragments according to desired percentage;
* 'contamination' with background fragments; addition of the random reference genome fragments with correspondent length and desired relative number;
* size selection: selection of fragments fitted to desirable size;
* sequencing of the fragments from positive and negative strands: 
cutting the 5’end of the fragment of desirable length, or the 3’end and reversing it (by random choice, in single end mode), or both ends (in paired end mode).

The input parameters of simulation process (fragment size distribution, background levels, amplification coefficient, type of fragment library etc.) are adjusted correspondently to those in real ChIP-Seq experiment.

The model was developed by [Dr. Tatiana Subkhankulova](https://www.linkedin.com/in/tatiana-subkhankulova-0876a240), Imperial College London.

## Fragments distribution and size selection
The reference sequence is cut into fragments, which have their sizes drawn from a lognormal distribution. 
Then each fragment goes through the size selection filter. 
According to some publications, in reality size selection has normal distributed view. 
But the strong implementation of this idea leads to great loss of fragments, and as a consequence, an increase in run time without any sense in output. 
For this reason, **isChIP** implemented an 'expanded' pseudo-normal distribution of size selection, 
when reads are filtered only at the offset edges of the distribution, 
as it is shown in the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/sizeSelFilter.png). 
By using default distribution parameters, in dark blue – distribution of generated fragments, 
in light blue – hypothetical real size selection, in green – pseudo size selection.<br>
Lognormal distribution is implemented as X=e^(Y*factor+term), 
where Y is stated as a normal distributed value.<br>
Accordingly, it is managed by 4 options: ```--mean``` and ```--sigma``` are defined the normal random generator, 
and ```--ln-factor``` and ```--ln-term``` are specified the lognormal outlet.<br>
Size selection filter is managed by 3 options: 
```--sz-sel-sigma``` and ```--frag-len``` response to standard deviation and mean in standard normal distribution, 
and ```--frag-dev``` is half-width on which the standard distribution is 'expanded' (**d** on the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/sizeSelFilter.png)).

To visualize the lognormal distribution for different values of these parameters, 
use the [RandomTest](https://github.com/fnaumenko/RandomTest-Win) application (so far only under Windows).

##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
