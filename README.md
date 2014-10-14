## Swiss

* [Synopsis](#synopsis)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
  * [Simple example](#simple-example)
  * [Genome build](#genome-build)
  * [Association result formats](#association-result-formats)
    * [Simple format](#simple-format)
    * [EPACTS multi-assoc format](#epacts-multi-assoc-format)
  * [LD sources ](#ld-sources)
  * [Filtering results](#filtering-results)
  * [Clumping options](#clumping-options)
    * [LD based clumping](#ld-based-clumping)
    * [Distance based clumping](#distance-based-clumping)
  * [GWAS catalogs](#gwas-catalogs)
  * [GWAS catalog lookups](#gwas-catalog-lookups)
  * [Output from Swiss](#output-from-swiss)
  * [Common command-lines used](#common-command-lines-used)
* [Options](#options)
* [Limitations](#limitations)
* [License](#license)

## Synopsis

Swiss is a tool for pruning association scan results from a GWAS or sequencing study, and identifying regions near or in LD with previously reported GWAS signals. 

Swiss implements the following procedure:

* Prune a list of variants using LD or distance, keeping the best variant by p-value (very similar to PLINK's method.) 
* Identify which of the pruned variants are near, or in LD, with previously reported GWAS signals

Swiss supports two main formats: 

* A tab-delimited file of association results with the usual columns (CHROM, POS, SNP, PVAL)
* An EPACTS multi-assoc file containing association p-values across a number of traits

## Requirements

* Python 2.7 (not 3.x)
* Linux (tested on Ubuntu)

## Download

The latest "binary" downloads are here:

| Version | Date       | Data for LD          | Size | File                                                        |
|---------|------------|----------------------|------|-------------------------------------------------------------|
| 0.9.2   | 10/12/2014 | Yes                  | 3.8G | http://csg.sph.umich.edu/boehnke/swiss/swiss_0.9.2.tgz      |
| 0.9.2   | 10/12/2014 | No                   | 473K | http://csg.sph.umich.edu/boehnke/swiss/swiss_nold_0.9.2.tgz |

This Github site will always contain the most recent source code, and will usually be slightly ahead of the binary/packaged versions listed above. 

## Installation

* Extract the tarball wherever you'd like
* Run the bin/setup.py script

## Usage

### Simple example

```bash
swiss --assoc my_file.txt --ld-clump --clump-p 5e-08 --out my_results
```

### Genome build

You should always specify which genome build you're working in by using `--build`. By default, the build is hg19. 

Additionally, if you specify your own GWAS catalog, or VCF files for calculating LD, you should verify that the positions for these match the genome build of your association results. 

### Association result formats

#### Simple format

The simplest format looks like your typical association results: 

| CHROM | POS | SNP       | PVALUE |
|-------|-----|-----------|--------|
| 1     | 1   | rs123     | 5e-08  |
| X     | 2   | X:2_A/C   | 6e-09  |
| Y     | 314 | Y:314_C/G | 1e-11  |

You can specify the delimiter with `--delim` and the names of the columns with `--snp-col`, `--chrom-col`, `--pos-col`, `--pval-col`. The defaults are listed below under options. If you're analyzing multiple files, 1 per trait, you may want to tell swiss the name of your trait using `--trait <trait>`. This will include a TRAIT column in your output, which can be useful for joining results together later. 

The file can be gzipped. 

#### EPACTS multi-assoc format

Additionally, you can tell Swiss that your file is an EPACTS multi-assoc file with the `--multi-assoc` flag. This type of file looks like the following: 

| #CHROM | BEG   | END   | MARKER_ID      | NS   | AC       | CALLRATE | GENOCNT  | MAF     | TRAIT1.P | TRAIT1.B | TRAIT2.P | TRAIT2.B |
|--------|-------|-------|----------------|------|----------|----------|----------|---------|----------|----------|----------|----------|
| 1      | 15903 | 15903 | 1:15903_G/GC   | 8448 | 14459.66 | 1        | 0/3/8445 | 0.1442  | 0.5      | 0.195    | 0.659    | 0.128    |
| 1      | 19190 | 19191 | 1:19190_GC/G   | 8448 | 98.23    | 1        | 8448/0/0 | 0.00581 | 0.703    | 0.266    | 0.588    | -0.379   |
| 1      | 20316 | 20317 | 1:20316_GA/G   | 8448 | 120.46   | 1        | 8448/0/0 | 0.00713 | 0.714    | -0.512   | 0.645    | 0.644    |
| 1      | 30967 | 30970 | 1:30967_CCCA/C | 8448 | 47.35    | 1        | 8448/0/0 | 0.0028  | 0.322    | 3.15     | 0.296    | 3.32     |
| 1      | 51972 | 51975 | 1:51972_GGAC/G | 8448 | 268.34   | 1        | 8448/0/0 | 0.01588 | 0.673    | 0.301    | 0.866    | -0.121   |
| 1      | 53138 | 53140 | 1:53138_TAA/T  | 8448 | 402.05   | 1        | 8448/0/0 | 0.0238  | 0.368    | -0.768   | 0.905    | -0.103   |
| 1      | 54421 | 54421 | 1:54421_A/G    | 8448 | 422.81   | 1        | 8448/0/0 | 0.02502 | 0.367    | -0.776   | 0.98     | -0.0215  |
| 1      | 66221 | 66221 | 1:66221_A/AT   | 8448 | 338.19   | 1        | 8448/0/0 | 0.02002 | 0.0378   | 1.24     | 0.211    | 0.747    |
| 1      | 66222 | 66223 | 1:66222_TA/T   | 8448 | 298.81   | 1        | 8448/0/0 | 0.01769 | 0.0653   | 1.13     | 0.314    | 0.615    |

There are a set of columns (.P, .B) for each trait that was analyzed. The file is tab-delimited, and gzipped. 

Example command line: 

```bash
swiss --assoc results.epacts.gz --multi-assoc --out my_results 
```

By default, swiss will try to run on every single trait given in the file. However, if you only wish to look at a single trait, you can use `--trait` instead: 

```bash
swiss --assoc results.epacts.gz --multi-assoc --out my_results --trait TRAIT1
```

If you're running on a machine with multiple CPU cores, you can ask swiss to do multiple traits from the multi-assoc file at the same time by telling it how many to run with `-T <num of parallel jobs>`. Please remember these run on the same machine, and not on the cluster - **do not overwhelm the machine!**

### LD sources 

Swiss comes with a few built-in sources of LD information: 

```bash
swiss --list-ld-sources

Build      LD Sources
-----      ----------
hg19       1000G_2012-03_AFR, 1000G_2012-03_AMR, 1000G_2012-03_ASN, 1000G_2012-03_EUR, GOT2D_2011-11
```

You can select different sources to use when LD pruning results, and when looking for GWAS catalog variants in LD. For example, you may wish to use your own genotypes for pruning (since they will cover all of your markers), but when looking for GWAS catalog variants in LD, it may be better to use a reference panel such as GoT2D for better coverage of your novel variants + known GWAS variants. 

* For the pruning step, use: `--ld-clump-source <name>`. 
* For the GWAS catalog LD lookup step, use: `--ld-gwas-source <name>`. 

Both options can be the same (and in fact, if you only specify one of them, *it assumes you meant to use that source for both.*)

You can always provide a VCF directly to use instead of selecting a built-in one: 

```bash
swiss --ld-clump-source /path/to/vcf.gz
```

If you have multiple VCF files split up across chromosomes, you can specify a .json file that maps chromosomes to VCF files: 

```bash
swiss --ld-clump-source /path/to/vcfmap.json
```

Where the `vcfmap.json` file looks like: 

```
{
  "1": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr1.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "10": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr10.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "11": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr11.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "12": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr12.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "13": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr13.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "14": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr14.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "15": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr15.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "16": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr16.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "17": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr17.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "18": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr18.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "19": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr19.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "2": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr2.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "20": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr20.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "21": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr21.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "22": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr22.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "3": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr3.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "4": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr4.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "5": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr5.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "6": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr6.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "7": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr7.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "8": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr8.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "9": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chr9.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz",
  "X": "/net/got2d/cfuchsb/T2Dgo/paper/data/2657/GoT2D.chrX.final_integrated_snps_indels_sv_beagle_thunder.qc.vcf.gz"
}
```

JSON format is a little fussy, so be careful. Make sure to use **double quotes** like above. 

### Filtering results

If you provided an imputation quality column in your association results (specified with `--rsq-col`), swiss can remove variants below a certain threshold using `--rsq-filter <threshold>`. 

### Clumping options

#### LD based clumping

Swiss can clump your association results using LD. The result being that only the best variants by p-value are kept first, and the remaining variants in LD with it are dropped. 

```bash
swiss --ld-clump --ld-clump-source GOT2D_2011-11 --clump-ld-thresh 0.8 --clump-p 4e-09
```

In the example above, variants in LD (r2) > 0.8 with the top variant per region are removed, and only variants with a p-value < 4e-09 are  considered at all. 

#### Distance based clumping

Similarly, you can prune based on distance. The best variants by p-value are retained, the remaining variants within X distance are dropped, and this process is continued until no variants remain to be considered. 

```bash
swiss --dist-clump --clump-dist 250000
```

In the example, variants within 250kb of the best p-value variant are removed, and so forth. 

### GWAS catalogs

Swiss supports two types of GWAS catalogs: built-in ones that come with the program, and user-supplied catalogs. 

The built-in catalogs can be found by doing: 

```bash
swiss --list-gwas-cats

Build      Catalog
-----      -------
hg18       fusion
hg19       fusion
```

Then you can select the catalog to use by `--gwas-cat fusion`, for example. Build is selected with `--build hg19`. 

The fusion catalog is an internal one maintained by our group here. 

If you'd like a list of traits contained in a particular catalog: 

```bash
swiss --list-gwas-traits

Available traits for GWAS catalog 'fusion':

APOA1B
------

ApoA1
ApoB
ApoB/ApoA1

Amino acids clumped
-------------------

2-aminobutyrate
2-hydroxyisobutyrate
3-(4-hydroxyphenyl)lactate
3-(4-hydroxyphenyl)lactate/ alpha-hydroxyisovalerate
3-phenylpropionate (hydrocinnamate)
4-acetamidobutanoate/ X-03056
5-oxoproline
```

You can also specify your own GWAS catalog by giving a filename instead of a codename for the catalog, like: `--gwas-cat /path/to/my/gwascat.tab`. 

The GWAS catalog format looks like the following (tab-delimited):

| SNP         | CHR | POS      | Group       | PHENO  | P_VALUE  |
|-------------|-----|----------|-------------|--------|----------|
| rs59793352  | 11  | 61420989 | Fatty acids | Bis/FA | 1.11E-11 |
| rs3932169   | 11  | 61420990 | Fatty acids | Bis/FA | 1.29E-15 |
| 11-61220700 | 11  | 61464124 | Fatty acids | Bis/FA | 4.38E-11 |
| rs2240286   | 11  | 61506454 | Fatty acids | Bis/FA | 1.16E-13 |
| 11-61288386 | 11  | 61531810 | Fatty acids | Bis/FA | 2.09E-10 |
| rs149803    | 11  | 61539020 | Fatty acids | Bis/FA | 9.30E-17 |

It can contain additional columns, for example you may have citations along with each hit or other supporting information: 

| SNP         | CHR | POS      | CITATION                         | EFFECT_ALLELE | EFFECT_SIZE | Group       | PHENO  | POPULATION         | P_VALUE  |
|-------------|-----|----------|----------------------------------|---------------|-------------|-------------|--------|--------------------|----------|
| rs59793352  | 11  | 61420989 | Kettunen et al. (Nat Genet 2012) | G             | -0.14       | Fatty acids | Bis/FA | European (Finnish) | 1.11E-11 |
| rs3932169   | 11  | 61420990 | Kettunen et al. (Nat Genet 2012) | G             | -0.17       | Fatty acids | Bis/FA | European (Finnish) | 1.29E-15 |
| 11-61220700 | 11  | 61464124 | Kettunen et al. (Nat Genet 2012) | A             | -0.34       | Fatty acids | Bis/FA | European (Finnish) | 4.38E-11 |
| rs2240286   | 11  | 61506454 | Kettunen et al. (Nat Genet 2012) | C             | 0.15        | Fatty acids | Bis/FA | European (Finnish) | 1.16E-13 |
| 11-61288386 | 11  | 61531810 | Kettunen et al. (Nat Genet 2012) | G             | -0.29       | Fatty acids | Bis/FA | European (Finnish) | 2.09E-10 |
| rs149803    | 11  | 61539020 | Kettunen et al. (Nat Genet 2012) | G             | 0.23        | Fatty acids | Bis/FA | European (Finnish) | 9.30E-17 |

The extra columns will be included with the output from Swiss. 

### GWAS catalog lookups

After LD or distance based clumping, Swiss will look for GWAS catalog hits that are near, or in LD, with your clumped/top variants. It does both and generates two files, one for each: 

* prefix.ld-gwas.tab - file contains GWAS catalog variants that were in LD with your top variants after clumping
* prefix.near-gwas.tab - contains GWAS catalogs near your top variants by distance

You can control the LD threshold using `--gwas-cat-ld <threshold>` and distance threshold using `--gwas-cat-dist <threshold>`. 

Swiss normally only includes columns from the GWAS catalog (as well as a few relevant columns from your association results) in these files. If you want to include additional columns from your assoc file:

```bash
swiss --assoc my_assoc.txt --include-cols "RSQ,EFF_AL,EFF_FREQ"
```

### Output from Swiss

Swiss generates the two GWAS catalog lookup files (listed above), and a third .clump file containing your top variants after clumping. The files are named starting with a prefix given by `--out`, for example: 

```bash
swiss --assoc my_assoc.txt --ld-clump --out prefix
```

Will create: 

* prefix.ld-gwas.tab
* prefix.near-gwas.tab
* prefix.clump

The .clump file looks like this: 

| #CHROM | BEG      | END      | MARKER_ID                | PVALUE   | BETA    | MRSQ    | TRAIT  | ld_with                                                                  | ld_with_values | failed_clump |
|--------|----------|----------|--------------------------|----------|---------|---------|--------|--------------------------------------------------------------------------|----------------|--------------|
| 11     | 60784275 | 60784275 | 11:60784275_G/A          | 4.47E-08 | -0.0992 | 0.98842 | otPUFA |                                                                          |                | pass         |
| 11     | 60786289 | 60786289 | 11:60786289_C/T          | 3.64E-10 | -0.307  | 0.93654 | otPUFA |                                                                          |                | pass         |
| 11     | 60859791 | 60859791 | 11:60859791_C/T_rs175133 | 9.51E-11 | 0.118   | 0.99901 | otPUFA | 11:60899767_A/G_exm915580,11:60853986_A/G,11:60859624_A/C_SNP11-60616200 | 0.40,0.60,0.61 | pass         |
| 11     | 60866519 | 60866519 | 11:60866519_A/ACCCAG     | 1.49E-11 | -0.246  | 0.94861 | otPUFA |                                                                          |                | fail         |

The `ld_with` column gives a comma separated list of variants that were pruned away (if LD clumping was used.) The r2 values are given for each variant (in the same order) in the `ld_with_values` column. 

If a variant failed LD calculation for some reason (not present in the VCF file, variant was an indel, etc.) the `failed_clump` column will say **fail**. The program will also generate a warning while running. 

The .ld-gwas.tab and .near-gwas.tab files are very similar (removing some columns for brevity): 

| ASSOC_MARKER    | ASSOC_CHRPOS | ASSOC_TRAIT | GWAS_SNP   | GWAS_CHRPOS | ASSOC_GWAS_LD | GWAS_GENE_LABEL | GWAS_Group | GWAS_PHENO | GWAS_P_VALUE |
|-----------------|--------------|-------------|------------|-------------|---------------|-----------------|------------|------------|--------------|
| 15:58683366_A/G | 15:58683366  | TotFA       | rs4775041  | 15:58674695 | 0.54800787    | LIPC            | Lipids     | HDL        | 3.20E-20     |
| 15:58683366_A/G | 15:58683366  | TotFA       | rs4775041  | 15:58674695 | 0.54800787    | LIPC            | Lipids     | TG         | 1.60E-08     |
| 15:58683366_A/G | 15:58683366  | TotFA       | rs10468017 | 15:58678512 | 0.636239711   | LIPC            | Lipids     | HDL        | 8.00E-23     |
| 15:58683366_A/G | 15:58683366  | TotFA       | rs1532085  | 15:58683366 | 1             | LIPC            | Lipids     | HDL        | 1.00E-188    |

* ASSOC_MARKER: Variant from your clumped association results (the top hit.) 
* ASSOC_CHRPOS: CHR:POS naming for the variant
* ASSOC_TRAIT: Either taken from the multi-assoc file, or specified with `--trait`. 
* GWAS_SNP: The GWAS catalog variant that your ASSOC_MARKER is in LD with. 
* ASSOC_GWAS_LD: The r2 between the GWAS_SNP and the ASSOC_MARKER. 
* GWAS_PHENO: The phenotype associated with the GWAS_SNP according to the GWAS catalog. 
* GWAS_P_VALUE: P-value reported in the GWAS catalog. 

The .near-gwas.tab file has ASSOC_GWAS_DIST instead of ASSOC_GWAS_LD, and denotes the distance between the ASSOC_MARKER and the GWAS_SNP. 

### Common command-lines used

```bash
swiss --assoc example.multiassoc.epacts.gz --multi-assoc \
--build hg19 --ld-clump-source /net/snowwhite/home/welchr/projects/FFA/metsim_got2d_exomechip.json \
--ld-gwas-source /net/snowwhite/home/welchr/projects/FFA/metsim_got2d_exomechip.json \
--gwas-cat fusion --ld-clump --clump-p 5e-08 --out example
```

The command above will: 

* Run on an EPACTS multiassoc file (and do all traits. To do a single trait, use `--trait`). 
* Use LD clumping to prune variants, and use VCF files specified by metsim_got2d_exomechip.json to do it
* Remove any variant with p-value > 5e-08
* Use the FUSION GWAS catalog (our internal catalog) for looking up GWAS variants in LD with top signals
* Again use the VCFs specified by metsim_got2d_exomechip.json to find GWAS variants in LD with top signals

---

```bash
swiss --assoc my_results.tab --delim tab --chrom-col CHROM --pos-col POS --pval-col PVAL --snp-col SNP \
--rsq-col RSQ --rsq-filter 0.3 \
--build hg19 --ld-clump-source 1000G_2012-03_EUR --ld-gwas-source 1000G_2012-03_EUR \
--gwas-cat fusion --dist-clump --clump-p 5e-08 --clump-dist 500000 --out example
```

The command above will: 

* Run on a simple tab-delimited format of GWAS association results, specifying the column names directly
* Filter variants on imputation quality 0.3
* Clump results using distance of 500kb, and also remove variants with p > 5e-08
* Use 1000G EUR to both LD clump AND find GWAS variants in LD with top signals

## Options

```
  -h, --help
    show this help message and exit

  --assoc <string>
    [Required] Association results file.

  --multi-assoc
    Designate that the results file is in EPACTS multi-assoc format.
    Default value is: False

  --trait <string>
    Description of phenotype for association results file. E.g. 'HDL' or 'T2D'

  --delim <string>
    Association results delimiter.
    Default value is: tab

  --build <string>
    Genome build your association results are anchored to.
    Default value is: hg19

  --snp-col <string>
    SNP column name in results file.
    Default value is: MarkerName

  --pval-col <string>
    P-value column name in results file.
    Default value is: P-value

  --chrom-col <string>
    Chromosome column name in results file.
    Default value is: CHR

  --pos-col <string>
    Position column name in results file.
    Default value is: POS

  --rsq-col <string>
    Imputation quality column name.
    Default value is: RSQ

  --rsq-filter <string>
    Remove variants below this imputation quality.
    Default value is: None

  --filter <string>
    Give a general filter string to filter variants.
    Default value is: None

  --out <string>
    Prefix for output files.
    Default value is: swiss_output

  --ld-clump
    Clump association results by LD.
    Default value is: False

  --clump-p <float>
    P-value threshold for LD and distance based clumping.
    Default value is: 5e-08

  --clump-ld-thresh <float>
    LD threshold for clumping.
    Default value is: 0.2

  --clump-ld-dist <int>
    Distance from each significant result to calculate LD.
    Default value is: 1000000

  --dist-clump
    Clump association results by distance.
    Default value is: False

  --clump-dist <int>
    Distance threshold to use for clumping based on distance.
    Default value is: 250000

  --ld-clump-source <string>
    Name of pre-configured LD source, or a VCF file from which to compute LD.
    Default value is: GOT2D_2011-11

  --list-ld-sources
    Print a list of available LD sources for each genome build.
    Default value is: False

  --gwas-cat <string>
    GWAS catalog to use.
    Default value is: fusion

  --ld-gwas-source <string>
    Name of pre-configured LD source or VCF file to use when calculating LD with GWAS variants.
    Default value is: GOT2D_2011-11

  --list-gwas-cats
    Give a listing of all valid GWAS catalogs and their descriptions.
    Default value is: False

  --list-gwas-traits
    List all of the available traits in a selected GWAS catalog.
    Default value is: False

  --list-gwas-trait-groups
    List all of the available groupings of traits in a selected GWAS catalog.
    Default value is: False

  --gwas-cat-p <float>
    P-value threshold for GWAS catalog variants.
    Default value is: 5e-08

  --gwas-cat-ld <float>
    LD threshold for considering a GWAS catalog variant in LD.
    Default value is: 0.1

  --gwas-cat-dist <int>
    Distance threshold for considering a GWAS catalog variant 'nearby'.
    Default value is: 250000

  --include-cols <string>
    List of columns to merge in from association results (grouped by variant.)
    Default value is: None

  --skip-overlap-check
    Skip the check of whether the GWAS catalog has variants that are not in your --ld-gwas-source.
    Default value is: False

  --cache <string>
    Prefix for LD cache.
    Default value is: ld_cache

  -T, --threads <int>
    Number of parallel jobs to run. Only works with --multi-assoc currently.
    Default value is: 1
```

## Limitations

Currently does not support indels (or to be more exact, indels are simply skipped and left alone in your results, and flagged with a **fail**. 

The latest human genome build (hg38) is not supported. 

## License

Copyright (C) 2014 Ryan Welch, The University of Michigan

Swiss is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Swiss is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

