# Contig scaffold comparisons

The purpose of this program is to allow direct comparison of contig and scaffold files from [SPAdes](http://cab.spbu.ru/software/spades/) outputs. If you want to compare contig and scaffold files now, it is very difficult since the NODE names don't match.

## Dependencies

The only dependencies are a couple of Perl modules I regularly use. I think these are actually available in the Perl download, so not sure you really need to install anything. The modules are: ```strict```, ```Getopt::Std```, and ```List::MoreUtils```.

## Installation

Download the repository from GitHub: [contig\_scaffold\_comparisons](https://github.com/seanmcallister/contig_scaffold_comparisons)

## Usage

After installation, pull up the help information:

```
perl contig_to_scaffold_comparitor.pl -h
```
This will give you:

```
Help called:
Options:
-c = contigs.paths file from spades
-f = contigs.fasta file from spades
-q = QC'ed contig fasta file
-s = scaffolds.fasta file from spades
-p = scaffolds.paths file from spades
-m = sample header prefix for simplified headers
-l = minimum sequence length to export (Note: will export two files: one with all sequences and one to this minimum length)
-h = This help message
```

Explanation of inputs:

1. contigs.paths, contigs.fasta, scaffolds.fasta, and scaffolds.paths are all automatically produced by SPAdes.
2. The QC'ed contig fasta file (```-q```), just a fasta file of the assembly with contigs passing your QC. For me, I QC contigs to remove contigs that don't have at least 1X read coverage over 90% of their length, this for reads mapping to their own assembly (i.e. the reads that went into the assembly in the first place).
3. Sample header prefix will export a file with the prefix + NODE + the node number (i.e. S1\_NODE\_1001). This is perfect for use with binning programs.
4. The minimum length is helpful if you are looking to put a smaller subset into a binning program.


### Example

Here is an example command:

```
contig_to_scaffold_comparitor.pl -c ../contigs.paths -f ../contigs.fasta -q ../contig_qc/S1_qccontigs_simpname.fasta -s ../scaffolds.fasta -p ../scaffolds.paths -m S1 -l 1000
```

While running, you should get print statments as you pass waypoints:

```
   Reading contig fasta file . . . 
   Reading scaffold fasta file . . . 
   Reading qc fasta file . . . 
   Reading contig paths file . . . 
   Reading scaffold paths file . . . 
   Tagging contigs dictionary with QC pass/fail . . . 
   Identifying multicontig scaffolds . . . 
   Adding new headers for single contig scaffolds . . . 
   Processing multicontig scaffolds . . . 
   Checking for problems . . . 
   Saving outfiles . . . 
```

The most problematic part of this comparison is when you get to multicontig scaffolds (scaffolds with more than one contig). Essentially, we are matching scaffold and contig names based on the EDGE names (which stay the same between the two paths files). Sometimes a single contig will consist of multiple EDGEs, so we'll get an error, such that the predicted number of contigs within a scaffold (based on the number of N's in the sequence) is different from the "actual" number (which would ideally be 1:1 EDGE to contig). I've made the program go through one pass to try to fix this, which works when EDGEs are exclusive to contigs, but sometimes two different contigs can share 1 or more EDGEs. This gets a bit complicated to program, so I've made it display the options and have the user pick (I know a bit archaic, but easy enough).

This is an example where the problem was fixed automatically (it will go to the next problematic multicontig scaffold):

```
WARNING: Scaffold does not match the predicted number of contigs: NODE_131_length_18414_cov_126.217     Predicted = 2   Actual = 11
Scaffold to contig match fixed
Scaffold edges: 15731,19763,37746,37747,41622,56629,68457
Contig match: NODE_3096_length_4310_cov_127.455 NODE_285_length_13566_cov_132.077
Contig edges: 15731,19763,37746,37747,41622,68457       56629   
```

This is an example where the user must decide which matches to keep:

```
WARNING: Scaffold does not match the predicted number of contigs: NODE_1456_length_6427_cov_12.7279     Predicted = 2   Actual = 3
NODE_1456_length_6427_cov_12.7279       THERE'S STILL A PROBLEM Predicted = 2   Actual = 3
Let's get some human input here:
Scaffold node edges: 7688,7689,17696
Possible matches are numbered:
   Option 1: NODE_7615_length_2482_cov_12.448   7688,17696
   Option 2: NODE_80098_length_132_cov_33.4     7688
   Option 3: NODE_3644_length_3935_cov_13.359   7689
Input contig matches to keep - comma delimited (i.e. 1,4,5)
$
```
At the ```$```, the program will wait for input (**it doesn't necessarily print a ```$```**). The proper input here is ```1,3```, which tells the program to keep options 1 and 3. Why these? There are only 2 contigs in the scaffold, and options 1 and 3 cover all three of the node's EDGEs, so is the proper match.

Here is another example:

```
WARNING: Scaffold does not match the predicted number of contigs: NODE_1413_length_14616_cov_3.59266    Predicted = 4   Actual = 6
NODE_1413_length_14616_cov_3.59266      THERE'S STILL A PROBLEM Predicted = 4   Actual = 6
Let's get some human input here:
Scaffold node edges: 21777,21778,106655,106657,243153,339709,339711
Possible matches are numbered:
   Option 1: NODE_27323_length_3001_cov_4.81663 21777,339709
   Option 2: NODE_363734_length_667_cov_6.11667 21778
   Option 3: NODE_17531_length_3865_cov_3.11851 106655,106657,243153
   Option 4: NODE_442491_length_136_cov_2.55556 106657
   Option 5: NODE_434220_length_221_cov_23.5957 243153
   Option 6: NODE_5816_length_7053_cov_3.35706  339711
Input contig matches to keep - comma delimited (i.e. 1,4,5)
$
```

And the answer: ```1,2,3,6```

This is a bit tedious if there are a lot of multicontig scaffolds with problems. I may fix it in future (but likely not).

#### Outfiles

Several outfiles are produced:

1. ```scaffolds_not_passing_QC.txt```
2. ```scaffolds_old_new_names.txt```
3. ```scaffolds_multicontig_stats.txt```
4. ```scaffolds_newnames.fasta```
5. ```scaffolds_newnames_simplified.fasta```
6. ```scaffolds_newnames_QCpass.fasta```
7. ```scaffolds_newnames_simplified_QCpass.fasta```
8. ```scaffolds_newnames_QCpass_atleastXbp.fasta```
9. ```scaffolds_newnames_simplified_QCpass_atleastXbp.fasta```

Their purpose:

1. List of scaffolds (old and new names, length) that don't contain contigs that pass the original QC run.
2. List of scaffold old and new names.
3. List of multicontig scaffolds (new names) along with all the contigs that are in the scaffolds (and their QC result).
4. FASTA with **ALL** scaffolds with their new names.
5. FASTA with **ALL** scaffolds with their new **simplified** names (i.e. S1\_NODE\_1001).
6. FASTA with **QC passing** scaffolds with their new names.
7. FASTA with **QC passing** scaffolds with their new **simplified** names.
8. Same as #6, except filtering out scaffolds that aren't greater or equal to the min length cutoff.
9. Same as #7, except filtering out scaffolds that aren't greater or equal to the min length cutoff.

Good luck!!

Please open an issue if you have any issues or questions. This program is provided without warranty or a guarantee of support. Thanks!!
