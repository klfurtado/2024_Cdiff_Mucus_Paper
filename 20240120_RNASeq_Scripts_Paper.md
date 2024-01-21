# RNA-Seq Workflow for *C. difficile*

RNA-Seq read processing was performed using a SLURM scheduler and the remote computing services offered by UNC-Chapel Hill. Below are the scripts used to generate the read count matrix used in the 2024 paper. Read cleaning was initially performed in 2020-2021, so scripts may need to be modified for newer versions of bioinformatic tools.

Files moved into and out of the remote workspace using `scp -r`

## 1. Check Sequence Quality with FASTQC

Script for FASTQC:

```
#!/bin/bash
##

#SBATCH --ntasks=1
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_202006/FastQC/run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_202006/FastQC/run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 1-00:00:00

echo "Loading FastQC Module"

module load fastqc/0.11.8

echo "FastQC loaded."
echo "Running FastQC."

fastqc /pine/scr/k/f/kfurtado/RNASeq_202006/seqs/RNASeq_202003/* -o /pine/scr/k/f/kfurtado/RNASeq_202006/FastQC/

echo "Finished. Check relevant output and error files in pwd."
echo "Transfer FastQC files to local computer using scp -r."
```
## 2. Trim Poor Quality Ends and Illumina Adapters with Trimmomatic
Before running script, ensure the paths are set to the correct location for the sequence files, or copy the sequences files to the current directory. Also ensure that the FASTA adapter file from the [Trimmomatic Github page](https://github.com/usadellab/Trimmomatic) is loaded in the correct directory.

Content of my adapters.fa file:

```
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

Script to run Trimmomatic:

```
#!/bin/bash
##

#SBATCH --ntasks=1
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Trimmomatic/run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Trimmomatic/run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 1-00:00:00

echo "Loading Trimmomatic Module"

module load trimmomatic/0.36

echo "Trimmomatic loaded."
echo "Running Trimmomatic to loop through all 9 pairs of files."
echo "Trimmomatic being run in PE mode, using palindromic adapter trimming to remove TruSeq3 adapters. Adapter file available from Trimmomatic GitHub page."

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C
do
	time trimmomatic PE -phred33 -basein /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Sequences/${basename}_R1_001.fastq.gz -baseout ${basename}_filtered.fastq.gz ILLUMINACLIP:/pine/scr/k/f/kfurtado/RNASeq_202006/Trimmomatic/adapters.fa:1:30:10:2:true SLIDINGWINDOW:5:25 MAXINFO:50:0.8 LEADING:30 TRAILING:25 MINLEN:40
done

echo "Finished running trimmomatic. Check .out and .err files in pwd for information."
echo "For each pair of input forward and reverse reads, there should be 4 output files, for paired and unpaired from each direction." 
```

### 2b. Test effectiveness of trimming using FASTQC
Use scripts similar to those step 1 to check.

## 3. Create an environment for running Fastq_screen. 

```
module load anaconda
conda create -n fastq-screen-env -c bioconda fastq-screen
``` 
**Note**: Only Fastq screen ended up being used in the final workflow. This step is only required if the remote computing cluster does not have Fastq_screen.

## 4. Remove Potential Contaminating genomes with Fastq_screen
To use Fastq_Screen, must first build indices of each potential contaminating genome using Bowtie2-build.

For Mucin RNA-Seq project, potential contaminating genomes are:

* *Homo sapiens* mitochondrion: NC_012920.1
* PhiX 174 Genome: NC_001422.1
* *Sus scrofa* mitochondrion: NC_000845.1

**Note**: For 2024 paper, only sequences containing human mucus were considered, so removing reads from *Sus scrofa* is technically not required. 

Genomes should be downloaded from NCBI in FASTA format and then transferred to the working directory in Longleaf via `scp -r`

Then, build genomes using bowtie2-build:

```
#!/bin/bash
##

#SBATCH --ntasks=1
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/build_run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/build_run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Loading bowtie2 to build indexes."

module load bowtie2/2.4.1

echo "Building index using phage genome from NCBI Genomes database."

time bowtie2-build -f ./nc001422.1_phix_174.fasta phix174

echo "Building index using pig mitochondrial genome from NCBI Genomes database."

time bowtie2-build -f ./nc000845.1_sus_scrofa_mt.fasta pig_mt

echo "Building index using human mitochondrial genome from NCBI Genomes database."

time bowtie2-build -f ./nc012920.1_homo_sapiens_mt.fasta human_mt

echo "Indexing complete. 6 basefiles should be created for each genome."
```

Then, created and moved all basefiles into separate directories for each reference genome. Example for pig_mt subdirectory: 

```
(base) [kfurtado@longleaf-login2 pig_mt]$ ls -l
total 8352
-rw-r--r-- 1 kfurtado users 4200087 Jul 28 14:21 pig_mt.1.bt2
-rw-r--r-- 1 kfurtado users    4160 Jul 28 14:21 pig_mt.2.bt2
-rw-r--r-- 1 kfurtado users      17 Jul 28 14:21 pig_mt.3.bt2
-rw-r--r-- 1 kfurtado users    4154 Jul 28 14:21 pig_mt.4.bt2
-rw-r--r-- 1 kfurtado users 4200087 Jul 28 14:21 pig_mt.rev.1.bt2
-rw-r--r-- 1 kfurtado users    4160 Jul 28 14:21 pig_mt.rev.2.bt2
```

Script for Fastq_screen:

A few notes:

 * The .conf file *must* be named "fastq_screen.conf" (note the underscore, not a dash).
 * The .conf file *must* be in the program directory (within the environment), not in the working directory.  

```
#!/bin/bash
##

#SBATCH --ntasks=1
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Load Fastq_screen environment from conda."

source activate fastq-screen-env

echo "Environment loaded."
echo "Loading bowtie2 v. 2.4.1 module required for Fastq_screen."

module load bowtie2/2.4.1

echo "Running fastq_screen to check trimmed reads for hits within PhiX and Sus Scrofa and Human mitochondrial genomes."

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C 
	do
	for suffix in 1P 1U 2P 2U
		do time /nas/longleaf/home/kfurtado/apps/miniconda3/envs/fastq-screen-env/share/fastq-screen-0.14.0-1/fastq_screen --threads 8 --subset 0 --aligner bowtie2 --force --nohits /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Trimmomatic/${basename}_filtered_${suffix}.fastq.gz
	done
done

echo "Reads should be screened, check output files."
echo "Deactivate Fastq_screen environment."

source deactivate fastq-screen-env

echo "Next, run FastQC for final assessment of quality."
```

Original .conf file: 

```
# This is an example configuration file for FastQ Screen

############################
## Bowtie, Bowtie 2 or BWA #
############################
## If the Bowtie, Bowtie 2 or BWA binary is not in your PATH, you can set 
## this value to tell the program where to find your chosen aligner.  Uncomment 
## the relevant line below and set the appropriate location.  Please note, 
## this path should INCLUDE the executable filename.

#BOWTIE	/usr/local/bin/bowtie/bowtie
#BOWTIE2 /usr/local/bowtie2/bowtie2
#BWA /usr/local/bwa/bwa



############################################
## Bismark (for bisulfite sequencing only) #
############################################
## If the Bismark binary is not in your PATH then you can set this value to 
## tell the program where to find it.  Uncomment the line below and set the 
## appropriate location. Please note, this path should INCLUDE the executable 
## filename.

#BISMARK	/usr/local/bin/bismark/bismark



############
## Threads #
############
## Genome aligners can be made to run across multiple CPU cores to speed up 
## searches.  Set this value to the number of cores you want for mapping reads.

THREADS		8



##############
## DATABASES #
##############
## This section enables you to configure multiple genomes databases (aligner index 
## files) to search against in your screen.  For each genome you need to provide a 
## database name (which can't contain spaces) and the location of the aligner index 
## files.
##
## The path to the index files SHOULD INCLUDE THE BASENAME of the index, e.g:
## /data/public/Genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
## Thus, the index files (Homo_sapiens.GRCh37.1.bt2, Homo_sapiens.GRCh37.2.bt2, etc.) 
## are found in a folder named 'GRCh37'.
##
## If, for example, the Bowtie, Bowtie2 and BWA indices of a given genome reside in 
## the SAME FOLDER, a SINLGE path may be provided to ALL the of indices.  The index 
## used will be the one compatible with the chosen aligner (as specified using the 
## --aligner flag).  
##
## The entries shown below are only suggested examples, you can add as many DATABASE 
## sections as required, and you can comment out or remove as many of the existing 
## entries as desired.  We suggest including genomes and sequences that may be sources 
## of contamination either because they where run on your sequencer previously, or may 
## have contaminated your sample during the library preparation step.
##
## Human - sequences available from
## ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/
#DATABASE	Human	/data/public/Genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
##
## Mouse - sequence available from
## ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/
#DATABASE	Mouse	/data/public/Genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37
##
## Ecoli- sequence available from EMBL accession U00096.2
#DATABASE	Ecoli	/data/public/Genomes/Ecoli/Ecoli
##
## PhiX - sequence available from Refseq accession NC_001422.1
#DATABASE	PhiX	/data/public/Genomes/PhiX/phi_plus_SNPs
##
## Adapters - sequence derived from the FastQC contaminats file found at: www.bioinformatics.babraham.ac.uk/projects/fastqc
#DATABASE	Adapters	/data/public/Genomes/Contaminants/Contaminants
##
## Vector - Sequence taken from the UniVec database
## http://www.ncbi.nlm.nih.gov/VecScreen/UniVec.html
#DATABASE	Vectors		/data/public/Genomes/Vectors/Vectors
```

Edited .conf file for genomes of interest and to specify path to the aligner:

```
# This is the EDITED configuration file for Fastq_screen
# Lines to be used have the # removed from the first line of the file.

############################
## Bowtie, Bowtie 2 or BWA #
############################
## If the Bowtie, Bowtie 2 or BWA binary is not in your PATH, you can set 
## this value to tell the program where to find your chosen aligner.  Uncomment 
## the relevant line below and set the appropriate location.  Please note, 
## this path should INCLUDE the executable filename.

#BOWTIE	/usr/local/bin/bowtie/bowtie
BOWTIE2 /nas/longleaf/home/kfurtado/apps/miniconda3/envs/flash2env/bin/bowtie2
#BWA /usr/local/bwa/bwa



############################################
## Bismark (for bisulfite sequencing only) #
############################################
## If the Bismark binary is not in your PATH then you can set this value to 
## tell the program where to find it.  Uncomment the line below and set the 
## appropriate location. Please note, this path should INCLUDE the executable 
## filename.

#BISMARK	/usr/local/bin/bismark/bismark



############
## Threads #
############
## Genome aligners can be made to run across multiple CPU cores to speed up 
## searches.  Set this value to the number of cores you want for mapping reads.

THREADS		8



##############
## DATABASES #
##############
## This section enables you to configure multiple genomes databases (aligner index 
## files) to search against in your screen.  For each genome you need to provide a 
## database name (which can't contain spaces) and the location of the aligner index 
## files.
##
## The path to the index files SHOULD INCLUDE THE BASENAME of the index, e.g:
## /data/public/Genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
## Thus, the index files (Homo_sapiens.GRCh37.1.bt2, Homo_sapiens.GRCh37.2.bt2, etc.) 
## are found in a folder named 'GRCh37'.
##
## If, for example, the Bowtie, Bowtie2 and BWA indices of a given genome reside in 
## the SAME FOLDER, a SINLGE path may be provided to ALL the of indices.  The index 
## used will be the one compatible with the chosen aligner (as specified using the 
## --aligner flag).  
##
## The entries shown below are only suggested examples, you can add as many DATABASE 
## sections as required, and you can comment out or remove as many of the existing 
## entries as desired.  We suggest including genomes and sequences that may be sources 
## of contamination either because they where run on your sequencer previously, or may 
## have contaminated your sample during the library preparation step.
##
## Human - sequences available from
## ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/
#DATABASE	Human	/data/public/Genomes/Human_Bowtie/GRCh37/Homo_sapiens.GRCh37
##
## Mouse - sequence available from
## ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/
#DATABASE	Mouse	/data/public/Genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37
##
## Ecoli- sequence available from EMBL accession U00096.2
#DATABASE	Ecoli	/data/public/Genomes/Ecoli/Ecoli
##
## PhiX - sequence available from Refseq accession NC_001422.1
DATABASE phix174	/pine/scr/k/f/kfurtado/RNASeq_202006/Fastq_screen/phix174/phix174
##
## Adapters - sequence derived from the FastQC contaminats file found at: www.bioinformatics.babraham.ac.uk/projects/fastqc
#DATABASE	Adapters	/data/public/Genomes/Contaminants/Contaminants
##
## Vector - Sequence taken from the UniVec database
## http://www.ncbi.nlm.nih.gov/VecScreen/UniVec.html
#DATABASE	Vectors		/data/public/Genomes/Vectors/Vectors
##
## Pig Mitochondria: NCBI Accession NC_000845.1
## Downloaded locally and transferred to Fastq_screen directory
DATABASE	pig_mt /pine/scr/k/f/kfurtado/RNASeq_202006/Fastq_screen/pig_mt/pig_mt
##
## Human Mitochondria: NCBI Accession NC_012920.1
## Downloaded locally and transferred to Fastq_screen directory
DATABASE human_mt /pine/scr/k/f/kfurtado/RNASeq_202006/Fastq_screen/human_mt/human_mt
```

After running Fastq_screen, files ending in _filter are to be used. The `--no_hits` option extracts genomes that do not align to any of the reference genomes.

## 5. Map merged and screened reads to *C.difficile* genome using Bowtie2

Using Bowtie2 aligner because this is a non-splicing aligner, which is more appropriate for prokaryotic samples. 

Will need to build *C. difficile* genome index in the working directory from which the aligner will be run.

For Bowtie2: http://gensoft.pasteur.fr/docs/bowtie2/2.2.6/#the-bowtie2-aligner

* Will need to build an index for *C. difficile* genome using ```bowtie2-build```. Store indexes in the same directory from which bowtie2 will run; it will look for index files by default here. Or, place them in a designated directory and specify said directory in the script.
* Local alignment or Global (end-to-end) alignment may be used. 
	* My understanding is that a local alignment is preferred with lower quality or low-signal data (e.g., metatranscriptomics samples from an environment, where the signal for any given microbe may be low). 
	* In our case, unclear whether local or global alignments would be better; global may be more accurate with low-complexity genomes, but we may lose more reads that are unable to align. Will test both.

**Test with local and global alignments (also a chance to observe output).**

Before running, move fasta file for the reference genome to the current directory and specify the correct filepath in the script. Or, run `bowtie2-build` separately before running the alignment.

*The script below is a bit redundant; it first builds the reference genome index in the current directory and then runs bowtie2, specifying a different index location. This is because the test script required a lot of troubleshooting to run; one troubleshoot was to move a previously-built index to a separate folder and specify the filepath there. That may be the best option moving forward, rather than building a new index every time.*

Chose to subset the data to compare the global and local methods a bit faster.

Creating the test dataset: 

```
#Went to Fastq_screen directory
#Created test file with the following
#Since each FASTQ file has each read split across 4 lines, pulling the first 400,000 lines results in 100,000 reads being pulled.

zcat 1A_filtered_1P.tagged_filter.fastq.gz | head -n 400000 > test.fastq.gz

#move the test file ot the mapping directory using `mv`
```

Test dataset containing 100,000 reads

```
#!/bin/bash
##

#SBATCH --ntasks=10
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/test.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/test.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Load Bowtie2 Module, v. 2.4.1"

module load bowtie2/2.4.1

echo "Build index for C. difficile genome, reference FN545816.1"

time bowtie2-build -f ./cdifficile_fn545816.1.fasta cdiff

echo "Testing Bowtie2 in Global mode."

bowtie2 -t -q --phred33 --very-sensitive --dovetail -x /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/cdiff/cdiff -S /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/globaltest.SAM -U "/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/test.fastq.gz" 

echo "Global Test Complete"

echo "Testing Bowtie2 in Local mode."

bowtie2 -t -q --phred33 --very-sensitive-local --dovetail -x /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/cdiff/cdiff -S /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/localtest.SAM -U "/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/test.fastq.gz" 

echo "Bowtie2 tests Completed"
echo "Check stderr files for the percentage of successful alignments with each method."
```


Results from test:

```
(base) [kfurtado@longleaf-login1 Mapping]$ cat test.err
     Genome indexes for use with BOWTIE2 are available
       in /proj/seq/data .
     bwa/bowtie/bowtie2 indexes are located in "Sequence" directory
       under the top level genome directory with name in CAPS, e.g. MM9_UCSC
Building a SMALL index

real	0m1.976s
user	0m1.801s
sys	0m0.069s
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 00:00:07
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    1641 (1.64%) aligned 0 times
    65105 (65.11%) aligned exactly 1 time
    33254 (33.25%) aligned >1 times
98.36% overall alignment rate
Time searching: 00:00:07
Overall time: 00:00:07
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 00:00:07
100000 reads; of these:
  100000 (100.00%) were unpaired; of these:
    1492 (1.49%) aligned 0 times
    63589 (63.59%) aligned exactly 1 time
    34919 (34.92%) aligned >1 times
98.51% overall alignment rate
Time searching: 00:00:07
Overall time: 00:00:07
```

Based on this, it seems performance is comparable between global and local, with a global test finding slightly fewer duplicate alignments than a local test.

Repeating the alignment in Global mode to avoid reads mapping more than once, this time using all reads:

```
#!/bin/bash
##

#SBATCH --ntasks=10
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Load Bowtie2 Module, v. 2.4.1"

module load bowtie2/2.4.1

echo "Using Bowtie2 in Global mode."

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C 
	do bowtie2 -t -q --phred33 --very-sensitive --dovetail -x /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/cdiff/cdiff -S /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/${basename}.SAM -U "/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/${basename}_filtered_1P.tagged_filter.fastq.gz","/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/${basename}_filtered_2P.tagged_filter.fastq.gz","/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/${basename}_filtered_1U.tagged_filter.fastq.gz","/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Fastq_screen/${basename}_filtered_2U.tagged_filter.fastq.gz"
done

echo "Bowtie2 Completed"
echo "Check stderr files for the percentage of successful alignments with each method."
```

With `--ntasks=10` specified, each iteration of the script for each basename takes ~1.5 hours and results in a .SAM file ~30 GB. Use `ls -lah` to observe progress and assess the size of each file as it's created.

**Results from the mapping step:**
 
```
(base) [kfurtado@longleaf-login5 Mapping]$ cat run.err
     Genome indexes for use with BOWTIE2 are available
       in /proj/seq/data .
     bwa/bowtie/bowtie2 indexes are located in "Sequence" directory
       under the top level genome directory with name in CAPS, e.g. MM9_UCSC
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:01
Multiseed full-index search: 01:32:45
85900599 reads; of these:
  85900599 (100.00%) were unpaired; of these:
    1573140 (1.83%) aligned 0 times
    58679982 (68.31%) aligned exactly 1 time
    25647477 (29.86%) aligned >1 times
98.17% overall alignment rate
Time searching: 01:32:46
Overall time: 01:32:50
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:25:42
81153652 reads; of these:
  81153652 (100.00%) were unpaired; of these:
    1980127 (2.44%) aligned 0 times
    59656644 (73.51%) aligned exactly 1 time
    19516881 (24.05%) aligned >1 times
97.56% overall alignment rate
Time searching: 01:25:42
Overall time: 01:25:43
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:29:30
84086631 reads; of these:
  84086631 (100.00%) were unpaired; of these:
    1823233 (2.17%) aligned 0 times
    62992115 (74.91%) aligned exactly 1 time
    19271283 (22.92%) aligned >1 times
97.83% overall alignment rate
Time searching: 01:29:30
Overall time: 01:29:30
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:30:28
82991912 reads; of these:
  82991912 (100.00%) were unpaired; of these:
    1296867 (1.56%) aligned 0 times
    56445849 (68.01%) aligned exactly 1 time
    25249196 (30.42%) aligned >1 times
98.44% overall alignment rate
Time searching: 01:30:28
Overall time: 01:30:28
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:36:05
85926759 reads; of these:
  85926759 (100.00%) were unpaired; of these:
    1575610 (1.83%) aligned 0 times
    61670000 (71.77%) aligned exactly 1 time
    22681149 (26.40%) aligned >1 times
98.17% overall alignment rate
Time searching: 01:36:05
Overall time: 01:36:05
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:22:16
73164793 reads; of these:
  73164793 (100.00%) were unpaired; of these:
    1373136 (1.88%) aligned 0 times
    55448201 (75.79%) aligned exactly 1 time
    16343456 (22.34%) aligned >1 times
98.12% overall alignment rate
Time searching: 01:22:16
Overall time: 01:22:16
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 01:44:46
97493607 reads; of these:
  97493607 (100.00%) were unpaired; of these:
    1828437 (1.88%) aligned 0 times
    66934232 (68.65%) aligned exactly 1 time
    28730938 (29.47%) aligned >1 times
98.12% overall alignment rate
Time searching: 01:44:46
Overall time: 01:44:46
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 02:05:40
107626101 reads; of these:
  107626101 (100.00%) were unpaired; of these:
    2253693 (2.09%) aligned 0 times
    81695309 (75.91%) aligned exactly 1 time
    23677099 (22.00%) aligned >1 times
97.91% overall alignment rate
Time searching: 02:05:40
Overall time: 02:05:40
Time loading reference: 00:00:00
Time loading forward index: 00:00:00
Time loading mirror index: 00:00:00
Multiseed full-index search: 02:03:25
96642214 reads; of these:
  96642214 (100.00%) were unpaired; of these:
    1977995 (2.05%) aligned 0 times
    73316194 (75.86%) aligned exactly 1 time
    21348025 (22.09%) aligned >1 times
97.95% overall alignment rate
Time searching: 02:03:25
Overall time: 02:03:25
```
##6. Convert .SAM to .BAM (necessary for running FADU in Step 7)

Samtools is available on the Longleaf Server.

For options regarding Samtools, see: http://www.htslib.org/doc/samtools.html

Script:

```
#!/bin/bash
##

#SBATCH --ntasks=10
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/convert.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/convert.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Loading SamTools module"

module load samtools/1.13

echo "Converting .SAM to .BAM and sorting the file by leftmost position."

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C 
	do samtools view -uS ${basename}.SAM | samtools sort -o ${basename}_sorted.bam
done

echo "Check directory for sorted .BAM files. Proceed to FADU."
```
##7. Running FADU to count genomic features (genes, operons, etc.).

More information on running FADU can be found here: https://igs.github.io/FADU/

Ultimately, I realized that the best way for me to run FADU in the server is to set up a Julia environment in the server, add the necessary packages within the environment, then call the environment within a script and run the program. A version of Julia for conda has been developed and is available on conda forge, but check whether any updates are available or needed for future use. 

[Julia on Conda Forge](https://anaconda.org/conda-forge/julia)

To set up environment: 

```conda create -n julia-fadu -c conda-forge julia```

Then, activate environment: 

```conda activate julia-fadu```

Once in the environment, opened a Julia interactive session and installed packages using: 

```
#Opens interactive session
julia

#Adding required packages
import Pkg; Pkg.add("ArgParse")
import Pkg; Pkg.add("BGZFStreams")
import Pkg; Pkg.add("GenomicFeatures")
import Pkg; Pkg.add("GFF3")
import Pkg; Pkg.add("StructArrays")
import Pkg; Pkg.add("Indexes")
import Pkg; Pkg.add("XAM")

#Check which packages and versions are installed 
import Pkg; Pkg.installed()

#exit interactive session
exit()
```

Then deactivate the environment for now with:

```conda deactivate```

Output on Package versions:

``` 
┌ Warning: Pkg.installed() is deprecated
└ @ Pkg /home/conda/feedstock_root/build_artifacts/julia_1632964725130/work/usr/share/julia/stdlib/v1.6/Pkg/src/Pkg.jl:570
Dict{String, VersionNumber} with 7 entries:
  "XAM"             => v"0.2.7"
  "Indexes"         => v"0.1.3"
  "StructArrays"    => v"0.6.0"
  "GFF3"            => v"0.2.1"
  "BGZFStreams"     => v"0.3.0"
  "GenomicFeatures" => v"2.0.4"
  "ArgParse"        => v"1.1.4"
```

Note, I did not have any issues running FADU with the most current versions of these packages. However, if future versions cause issues, you can specify a particular package version upon install with something like the following: 

```
#Copied from Julia forum:
Pkg.add(Pkg.PackageSpec(;name="xxx", version="1.2.3"))

#Example given our GenomicFeatures package:
import Pkg; Pkg.add(Pkg.PackageSpec(;name="GenomicFeatures", version="2.0.0"))
```

Set up environment with necessary packages. Received error that index files for bam files could not be found:

```
[ Info: Parsed args:
[ Info:   bam_file  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/3C_sorted.bam
[ Info:   gff3_file  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fn545816.1.gff3
[ Info:   stranded  =>  no
[ Info:   attribute_type  =>  ID
[ Info:   pp_only  =>  false
[ Info:   max_fragment_size  =>  1000
[ Info:   rm_multimap  =>  false
[ Info:   output_dir  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu_output/
[ Info:   em_iter  =>  1
[ Info:   feature_type  =>  gene
[ Info: Processing annotation features...
[ Info: Getting unique coordinates per contig and feature...
[ Info: Initializing dictionary of feature count information
ERROR: LoadError: Attempted to find index file for BAM file /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/3C_sorted.bam but could not find one
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:33
 [2] main()
   @ Main /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu.jl:147
 [3] top-level scope
   @ /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu.jl:206
in expression starting at /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu.jl:206

```

So, went back to index the bam files. In the future, this can be added to the sorting step.

```
#!/bin/bash
##

#SBATCH --ntasks=10
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/convert.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/convert.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Loading SamTools module"

module load samtools/1.13

echo "Indexing .bam files."

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C 
	do samtools index ${basename}_sorted.bam
done

echo "Check directory for sorted, indexed .BAM files. Proceed to FADU."
```

Then, re-ran FADU. Note that the lines starting with ##, excepting the header, can be deleted at this point. Also note that setting the variables as in this script can be used in previous scripts too, to make editing filepaths easier. 

```
(base) [kfurtado@longleaf-login5 FADU]$ cat fadu.sh
#!/bin/bash
##

#SBATCH --ntasks=10
#SBATCH -o /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/run.out
#SBATCH -e /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/run.err
#SBATCH --mem-per-cpu=10000
#SBATCH -t 8-00:00:00

echo "Setting Variables"

GFF="/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fn545816.1.gff3"
BAM="/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/"
OUTPUT="/pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu_output/"

##echo "Loading Julia & Python Modules"

##module load julia/1.6.2

echo "Open julia environment with required packages for FADU."

source activate julia-fadu

##echo "Activate Packages"

##julia import Pkg; Pkg.activate()

for basename in 1A 1B 1C 2A 2B 2C 3A 3B 3C 
	do julia --project=. fadu.jl -b $BAM${basename}_sorted.bam -g $GFF -o $OUTPUT
done

source deactivate julia-fadu

echo "Check output directory. Subset files and rearrange similar to htseq-count file for DESeq2."
```

**Output from FADU:** 

Text files containing a table of feature counts for each file are written to the ./fadu_output directory. First 50 lines from one table:

```
(base) [kfurtado@longleaf-login5 fadu_output]$ head -50 1A_sorted.counts.txt 
featureID	uniq_len	num_alignments	counts	tpm
gene-16S rRNA	1630	396743.0	395829.41	6241.39
gene-16S rRNA-2	1779	304514.5	303880.00	4390.23
gene-16S rRNA-3	1631	316205.0	315149.56	4966.19
gene-16S rRNA-4	1632	328029.0	319488.84	5031.49
gene-16S rRNA-5	1625	379795.5	379680.12	6005.17
gene-16S rRNA-6	1624	320660.0	320425.69	5071.10
gene-16S rRNA-7	1623	386315.0	385873.59	6110.65
gene-16S rRNA-8	1625	349669.0	349434.97	5526.80
gene-16S rRNA-9	1625	336510.0	336398.25	5320.61
gene-23S rRNA	2951	874654.0	874238.06	7614.15
gene-23S rRNA-10	2953	735283.0	734868.19	6395.98
gene-23S rRNA-2	2951	847814.5	847441.44	7380.76
gene-23S rRNA-3	2631	730358.5	727421.69	7106.02
gene-23S rRNA-4	2596	871081.0	764301.19	7566.95
gene-23S rRNA-5	2953	891101.0	890773.00	7752.90
gene-23S rRNA-6	2953	864049.5	863726.69	7517.50
gene-23S rRNA-7	2953	877452.5	876799.69	7631.29
gene-23S rRNA-8	2953	969697.0	969210.94	8435.59
gene-23S rRNA-9	2953	1022014.5	1021590.75	8891.49
gene-5S rRNA	117	67.5	26.54	5.83
gene-5S rRNA-2	117	377.5	57.39	12.61
gene-5S rRNA-3	117	22.0	16.52	3.63
gene-5S rRNA-4	117	86.0	21.24	4.67
gene-5S rRNA-5	117	22.0	15.66	3.44
gene-5S rRNA-6	117	51.5	12.26	2.69
gene-5S rRNA-7	117	29.5	26.28	5.77
gene-5S rRNA-8	117	50.5	38.71	8.50
gene-CDR20291_0001	510	2188.5	2031.96	102.40
gene-CDR20291_0002	882	2167.0	1976.02	57.58
gene-CDR20291_0003	1272	19710.5	19596.31	395.96
gene-CDR20291_0004	456	6018.5	5853.78	329.94
gene-CDR20291_0005	1638	11084.5	9906.72	155.45
gene-CDR20291_0006	351	2486.5	1690.51	123.79
gene-CDR20291_0007	600	12337.5	12159.07	520.85
gene-CDR20291_0008	1206	5282.5	5230.86	111.48
gene-CDR20291_0009	414	232.0	188.77	11.72
gene-CDR20291_0010	3432	140024.0	139906.91	1047.74
gene-CDR20291_0011	1929	1401.5	1370.18	18.26
gene-CDR20291_0012	465	2923.5	2238.52	123.73
gene-CDR20291_0013	495	5732.5	4877.74	253.26
gene-CDR20291_0014	1026	14499.5	13626.21	341.34
gene-CDR20291_0015	2250	51054.5	50618.39	578.21
gene-CDR20291_0016	1374	6022.0	5845.91	109.35
gene-CDR20291_0017	1071	1322.5	1186.62	28.48
gene-CDR20291_0018	1095	19384.5	19204.34	450.76
gene-CDR20291_0019	2178	7553.0	7297.50	86.11
gene-CDR20291_0020	846	109.0	108.27	3.29
gene-CDR20291_0021	1167	29.0	27.17	0.60
gene-CDR20291_0022	1647	64.0	60.74	0.95
```

Example Standard Error file for one file (same for all files):

```
[ Info: Parsed args:
[ Info:   bam_file  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/Mapping/1A_sorted.bam
[ Info:   gff3_file  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fn545816.1.gff3
[ Info:   stranded  =>  no
[ Info:   attribute_type  =>  ID
[ Info:   pp_only  =>  false
[ Info:   max_fragment_size  =>  1000
[ Info:   rm_multimap  =>  false
[ Info:   output_dir  =>  /pine/scr/k/f/kfurtado/RNASeq_Analysis_202006/FADU/fadu_output/
[ Info:   em_iter  =>  1
[ Info:   feature_type  =>  gene
[ Info: Processing annotation features...
[ Info: Getting unique coordinates per contig and feature...
[ Info: Initializing dictionary of feature count information
[ Info: Opening BAM alignment file...
[ Info: Now finding overlaps between alignment and annotation records...
[ Info: Counting and adjusting multimapped alignment feature counts via Expectation-Maximization algorithm...
[ Info: Writing counts output to file...
[ Info: FADU is complete!  Exiting!
```

Next steps: Extract the FeatureID and Counts columns from the output tables, then input to DESeq2!

### End of Unix/Bash Workflow. Next Steps to be done in R on local machine.
