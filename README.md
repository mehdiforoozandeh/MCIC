# MCIC
Using a sequence similarity based annotation and an ensemble machine learning approach, MCIC (metagenome cellulase identification and classification) aims to identify and classify cellulolytic enzymes from a given metagenomic data as well as any other amino-acid sequence on the basis of optimum temperature and pH. 

***USAGE:
(arguments are separated with a space)

	-firts argument : MCIC 

	-second argument options: [-h] [--help] [csp] [CelScreenPred] [cs] [CelScreen] [fp] [FastaPred] [sp] [SinglePred]

	-third argument: [query input file]

	-forth+ argument options: [-bs] [--bitscore] [-out] [-noexport]


***FUNCTIONS: 

1. [csp] or [CelScreenPred]: Accepts nucleotide sequences and screens the cellulolytic enzymes and predicts pH and temperature dependence.

2. [cs] or [CelScreen]: Accepts nucleotide sequences and screens the cellulolytic enzymes. (without prediction)

3. [fp] or [FastaPred]: Accepts a fasta file containing protein sequence of cellulolytic enzymes and predicts pH and temperature dependence.

4. [sp] or [SinglePred]: Accepts a single protein sequence or entry or accession number of cellulolytic enzymes and predicts pH and temperature dependence.


***HOW TO USE EACH FUNCTION:

1. [csp] or [CelScreenPred]: 

	options: 
	- [-bs]/[--bitscore]: choose your bitscore limit for filteration during the screening process. ==> "default: 50"
	- [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
	- [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

	<< By default, all results will be written in "results" folder.>>

	usage examples: 
		- MCIC csp input.format
		- MCIC CelScreenPred input.format -bs 100 -out John 
		- MCIC csp input.format --bitscore 500 -noexport
		- MCIC CelScreenPred input.format -out home/John
		- MCIC csp input.format -noexport

2. [cs] or [CelScreen]: 
	
	options: 
	- [-bs]/[--bitscore]: choose your bitscore limit for filteration during the screening process. ==> "default: 50"
	- [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
	- [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

	<< By default, all results will be written in "results" folder.>>

	usage examples: 
		- MCIC cs input.format
		- MCIC CelScreen input.format -bs 300 -out John 
		- MCIC cs input.format --bitscore 50 -noexport
		- MCIC CelScreen input.format -out home/John
		- MCIC cs input.format -noexport

3. [fp] or [FastaPred]:
	
	options:
	- [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
	- [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

	<< By default, all results will be written in "results" folder.>>

	usage examples: 
		- MCIC fp input.fasta
		- MCIC FastaPred input.fasta -out John
		- MCIC fp input.fasta -noexport

4. [sp] or [SinglePred]:
	
	This function does not have any options. 
	Two types of inputs are possible (1) amin acid sequence (2) protein entry/accession

	usage examples:
		- MCIC sp MKSCAILAALGCLA....
		- MCIC SinglePred Q7Z9M7
