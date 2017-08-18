# isChIP
**I**n **S**ilico **ChIP**-seq is a fast realistic ChIP-seq simulator.

The modelling of the chromatin immunoprecipitaion followed by next generation sequencing process is based on presumably on Illumina protocol. However, the **isChIP**’s options can be adjusted to the other techniques such as  Ion Torrent, etc.

The real protocol of ChIP-seq is simulated by repeating the basic cycle. Each basic cycle corresponds to single cell simulation, and consists of the next phases:
* “shearing” the chromatin, or random cutting the reference genome in fragments of correspondent length; the distribution of fragments by default corresponds to log-normal with correspondent parameters;
* “extraction” of the fragments overlapping with the binding events;
* amplification the selected fragments if required in number of cycles;
* “loss” of selected fragments according to desired percentage;
* “contamination” with background fragments; addition of the random reference genome fragments with correspondent length and desired relative number;
* size selection: selection of fragments fitted to desirable size;
* sequencing of the fragments from positive and negative strands: cutting the 5’end of the fragment of desirable length; addition to the output file.

Simulated binding events are specified by optional parameter, called [*template*](#template). (Also see [Modes](#modes) section)

The input parameters of simulation process (fragment size distribution,  total number of sequenced reads,  background levels, amplification coefficient, type of fragment library (single-end or paired-end, etc.) are adjusted correspondently to those in real ChIP-Seq experiment.

The model was developed by [Dr. Tatiana Subkhankulova](https://www.linkedin.com/in/tatiana-subkhankulova-0876a240), Imperial College London.

**isChIP** runs on the command line under Windows, Linux and Mac OS X.<br>

## Synopsis
```isChIP -g ref_genome –n 5```<br>
Generates “input” sequences in FastQ with read length of 50 and average read density about 9 read/kbs, comparable with experimentally one observed.

```isChIP -g mm9_dir –n 250 –b 4 –r 36 –f fq,sam templ.bed```<br>
Generates test sequences in FastQ and direct alignment in SAM, with read length of 36 and average foreground/background read density comparable with experimentally observed in [Series GSE56098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56098).<br>
templ.bed consists of lines with features defining well-proven Oct4 binding motif in the centre of some peaks in the referenced experimental data:
```
chr1   9899347   9899362  ...
chr1  16120224  16120239  ...
chr1  40360157  40360172  ...
. . .
```
Please note that correlation between the whole experimental and model alignment is incorrect due to the abundance of not simulated peaks and artifacts.

## Usage
```isChIP [options] -g|--gen <file> [template]```

## Help
```
Input:
  -g|--gen <name>       reference genome (file name or directory). Required
Processing:
  -a|--amplif <int>     coefficient of amplification [1]
  -b|--bg-level <float> background of out-of-features fragments, in percent of foreground.
                        For test mode only [0]
  --fg-level <float>    foreground of in-features fragments, in percent. For test mode only [100]
  -n|--cells <long>     number of cells [1]
  -c|--chr <chars>      generate output for stated chromosome only
  --flat-len <int>      boundary flattening length. For test mode only [0]
  --bg-all <OFF|ON>     turn on/off generation background for all chromosomes.
                        For test mode only [ON]
  --bind-len <int>      minimum binding length. For test mode only [1]
  --let-N               include the ambiguous reference characters (N) on the beginning
                        and on the end of chromosome
  --smode <SE|PE>       sequencing mode: SE - single-end, PE - paired-end [SE]
  --strand-admix <OFF|ON>       turn on/off opposite strand admixture at the bound of binding site.
                        For test mode only [OFF]
  --ts-uni              uniform template score. For test mode only
  -p|--threads <int>    number of threads [1]
  --debug               fix random emission to get repetitive results
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
  -r|--read-len <int>   length of output read [50]
  --read-name <NMB|POS> name of read in output files includes:
                        NMB - read`s unique number within chromosome
                        POS - read`s true start position [POS]
  --read-Nlimit <int>   maximum permitted number of ambiguous characters (N) in read [--read-len]
  --reads-limit <long>  maximum permitted number of total written reads [2e+08]
  --fq-qual <chars>     the quality values for the read in FQ output [~]
  --map-qual <int>      the mapping quality for the read in SAM output [42]
Output:
  -f|--format <FQ,BED,SAM>      format of output sequences/alignment, in any combination [FQ]
  -o|--out <name>       location of output files or existing directory
                        [Test mode: mTest.*, Control mode: mInput.*, Regular mode: mRegular.*]
  -z|--gzip             compress output files with gzip
Other:
  -t|--time             print run time
  -V|--verbose <CRIT|RES|RT|PAR|DBG>    set verbose level:
                        CRIT - show critical messages only (silent mode)
                        RES -  show result summary
                        RT  -  show run time information
                        PAR -  show process parameters
                        DBG -  show debug messages [RT]
  -v|--version          print program and ZLib version, and quit
  -h|--help             print usage information and quit
```

## Details

### Modes

**isChIP** generates output in one of three modes:<br>
*test* – simulation of site of interest sequencing; the output is test sequences/alignment<br>
*control* – simulation of control production; the output is “input” sequences/alignment<br>
*regular* – simple regular cutting of reference chromosome, an auxiliary mode for special use.<br>
*Test* and *control* modes are distinguished only by involvement or elimination of sites of interest in a process, stated in optional parameter called [*template*](#template).

### Template
*Template* represents a bunch of simulated binding events. It is a file in BED format.
Each data line contains 3 (minimum required) or 5 and more fields.
First 3 fields define the binding site.
If there are no more fields in the line, all generating significant enriched regions have the same density.
If 4th field (name) and 5th field (score) are in a line, **isChIP** uses feature’s score to restrict appropriate enriched density.
Unit does not matter, only scores ratio is meaningful.
This peculiarity allows to use BED files generated by peak callers, since different peak callers have different score units.
Nevertheless it is possible to ignore stated scores by option ```--ts-uni```.


Input files (reference genome and template) can be zipped by gzip.<br>

### Options description
Non-numeric option values are case insensitive.

```-g|--gen <name>```<br>
Reference sequence or directory contained reference sequences in FASTA format.<br>
If ```name``` is a directory, first **isChIP** searches .fa files in it.
If there are not such files in this directory, or file pointed by ```–c|--chr``` option is absent, isChIP searches .fa.gz files.

There are 3 ways to define desirable reference sequences.<br>
To generation output for the whole genome, name should be a path, and all chromosome files should be in a uniform compressed condition.<br>
To limit the generation of several chromosomes, name should be a path, and only desirable chromosome files should be decompressed or even present at all.<br>
To generate output for a single chromosome, name can be a path and option ```–c|--chr``` should be set, or name should indicate an appropriate reference file.<br>
In test mode the target references are determined by *template*. See also ```--bg-all``` and ```–c|--chr ``` options.

**isChIP** omits  'random' contigs and haplotype sequences.

UCSC genome libraries are available on ftp://hgdownload.soe.ucsc.edu/goldenPath.<br>
Ensembl genome libraries can be obtained at ftp://ftp.ensembl.org/pub/release-73/fasta.<br>
This option is required.

```-a|--amplify <int>```<br>
The coefficient of MDA-amplification of fragments passed through the size selection filter.<br>
The default value of 1 means non amplification.

```-b|--bg-level <float>```<br>
In *test* mode the level of background, in percent of foreground. The background is a set of reads cut from fragments that are not intersected with the *template* binding events. Fragments intersected with binding sites form the foreground.<br>
In practice the level of 2-5% corresponds with a good experimental data, while level more than 10% - with a foul one.<br>
Default: 1

```--fg-level <float>```<br>
The level of foreground, in percent of maximum possible. Used to reduce output sensitivity to number of cells, or to estimate the overall loss in the real process.<br>
Default: 100

```-n|--cells <long>```<br>
The number of “cells”. Each “cell” corresponds to one passage through the reference chromosome (dual for the numeric ones).<br>
It is NOT the equivalent of number of cells in the real experiment. In the real experiment the overall loss can be reached to 99%. By default value of ```--fg-level``` option the model does not lose any related fragment.<br>
The total program’s output heavy depends on correlation between established fragment size and parameters of fragment distribution and size selection. In general, by default distribution the number of 3-10 provides the output read mean density comparable with the actual “inputs”, and the number of 100-500 leads to data comparable with the actual tests in term of density.<br>
Default: 1

```-c|--chr <chars>```<br>
Generate output for stated chromosome only, for instance ```–c 10```, ```--chr X```. It takes the same effect as the referencing to the chromosome file instead of directory. It is a strong option, it abolishes the impact of option ```--bg-all``` and all another chromosomes from *template</i>.<br>
Default: all

```--flat-len <int>```<br>
In *test* mode the boundary flattening length is a distance from the boundary of binding site on which the probability of fragment associating is increased from zero to maximum. Such it simulates the smoothing of enriched regions.<br>
Default: 0

```--bg-all <OFF|ON>```<br>
In *test* mode turn on/off background generation for all chromosomes, irrespective of chromosomes included in *template*. As we discovered, mapping by any aligner one or several chromosomes to the whole reference genome leads to short local lacks of alignment, “gaps”. Generation background for all chromosomes eliminates this issue.
Please note that including background for all chromosomes is time consuming. There is no need to include all background by generating direct alignment in SAM or BED.<br>
Default: ON

```--bind-len <int>```<br>
In *test* mode the minimum binding length. That is a minimum number of nucleotides that ensures the binding while fragment intersects binding site. If some of the features in *template* are less then given binding length, they will be skipped and a summary warning message will appear. To see the detailed warning message about each skipped feature, set verbose level to debug (see ```--verbose```).<br>
Default: 1

```--let-N```<br>
Forces to include into scanning process the ambiguous reference characters (N) on the beginning and on the end of chromosome. It makes no difference in data after alignment, but increases a little a run time and a quality of random number distribution at the beginning of the process.

```--smode <SE|PE>```<br>
Generation reads according to sequencing mode: SE - single-end, PE - paired-end.<br>
Default: SE

```--strand-admix <OFF|ON>```<br>
In *test* mode turn on/off opposite strand admixture at the bound of binding site.<br>
In accordance with the nature of the binding of proteins, at the boundaries of the binding site, there should be reads with the same strand only. However, in practice, all the experimental sequences demonstrate at the boundaries of the binding sites the presence of small admixtures of reads with the opposite strand (see ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/mix-strand.png)). This option allows you to specify a "theoretical" or "practical" simulation.<br>
Default: OFF

```--ts-uni```<br>
In *test* mode ignores features scores if they are established in *template*. All features will be simulated with maximum score.

```-p|--threads <int>```<br>
Number of threads. The workflow is separated between chromosomes, so this option takes no effect by processing the single chromosome.

```--debug```<br>
Fix random emission to get repetitive results.

```-R|--regular <int>```<br>
*Regular* mode: write each read on starting position increased by stated shift.<br>
Used for specific tasks.

```--frag-len <int>```<br>
Average size of selected fragments.<br>
Default: 200

```--frag-dev <int>```<br>
Deviation of selected fragments. See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 20

```--sz-sel <OFF|ON>```<br>
Turn on/off fragment's size selection. See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: ON

```--sz-sel-sigma <int>```<br>
Standard deviation of the fragment's size selection normal distribution.<br>
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 20

```--mean <int>```<br>
Expectation of the fragment's size based normal distribution.<br>
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>Default: 200

```--sigma <int>```<br>
Standard deviation of the fragment's size based normal distribution.<br>
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 200

```--ln-factor <int>```<br>
Power multiplication factor in fragment's size lognormal distribution.<br>
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 500

```--ln-term <float>```<br>
Power summand in fragment's size lognormal distribution.<br>
See [Fragments distribution and size selection](#fragments-distribution-and-size-selection) section.<br>
Default: 5.1

```--read-name <NMB|POS>```<br>
Forces to include in the name of each read its unique number (NMB) or its true start position (POS).<br>
Default: POS

```--read-Nlimit <int>```<br>
Maximum permitted number of ambiguous reference characters (‘N’) in output reads.<br>
In the real sequences ambiguous characters could mind a failure of sequencing as well as undefined regions in the cut fragments, while simulated reads are full defined by reference library. Consequently the ambiguous characters in library are translated into reads. Different aligners considered undefined characters in references differently: some of them as invalid, some of them as overlapping case.<br>
This option allows to manage a level of this ambiguity.<br>
Default: length of read (all characters could be ‘N’)

```--reads-limit <long>```<br>
Maximum number of total written reads. Emulates sequenator’s limit.<br>
This value restricts number of written reads for each chromosome proportionally.<br>
Default: 200 000 000. In practical simulation this value is never achieved.

```--fq-qual <chars>```<br>
Quality value for all the positions in the read in FQ output.<br>
Default: ~ (maximal)

```--map-qual <int>```<br>
Mapping quality for the read in SAM output.<br>
Default: 42 (maximal)

```-f|--format <FQ,BED,SAM>```<br>
Output files formats. FQ states the read file. In paired-end mode two FQ files are generated, with suffixes ‘_1’ and ‘_2’.<br>
BED and SAM files state the immediate (direct) alignment. Any formats can be set, but as minimum one.<br>
Default: FQ

```-o|--out <file>```<br>
Output files location. If option’s value is a file, it is used as a common file name. If file has an extension, it is discarded. If value is a directory, the default file name is used.<br>
Default: *test* mode: **mTest.\***, *control* mode: **mInput.\***, *regular* mode: **mRegular.\***

## Fragments distribution and size selection
Reference sequence is cutting into fragments with lognormal distributed size. Then each fragment goes through the size selection filter.
According to the some publications, in reality size selection has normal distributed view.
But the strong implementation of this idea leads to great loss of fragments, i.e. run time increase, without any sense in output.
For this reason in **isChIP** the size selection has a pseudo normal distribution view, as it is shown on the ![figure](https://github.com/fnaumenko/isChIP/tree/master/auxil/sizeSelFilter.png).
By using default distribution parameters, in dark blue – fragment distribution, in light blue – hypothetical real size selection, in green – pseudo size selection.<br>
Lognormal distribution is implemented as X=e^(Y*factor+term)<br>
where Y is stated as a normal distributed value.<br>
Accordingly, it is managed by 4 options: ```--mean``` and ```--sigma``` are defined the normal random generator, and ```--ln-factor``` and ```--ln-term``` are specified the lognormal outlet.<br>
Size selection filter is managed by 3 options: ```--sz-sel-sigma``` and ```--frag-len``` response to standard deviation and mean in standard normal distribution, and   frag-dev is half-width on which the standard distribution is “expanded” (d on the ![figure](https://github.com/fnaumenko/isChIP/tree/master/pict/sizeSelFilter.png)).
