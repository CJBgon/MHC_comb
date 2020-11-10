# MHC_comb
Combine the output of several HLA haplotype calling tools to a concencus.
Takes .json files as input and prints every biologically possible combination of HLA class I, II or both alleles to be picked up in bash script 

## Dependencies:
* Itertools
* json
* argparse
* pandas
* numpy

## Usage:
For a short explanation per argument type `haplo_comb.py -h`

* `-a, --arcashla`: Path to the arcasHLA .json output.
* `-x, --xhla`: Path to the xHLA tsv output.
* `-s, --seq2hlai`: Path to MHC-I seq2HLA output: ClassI.HLAgenotype4digits.
* `-q, --seq2hlaii`: Path to MHC-II seq2HLA output: ClassII.HLAgenotype4digits.
* `-r, --resolution`: At what resolution should the HLA types be returned? Default = 3 -> 'HLA-A:02:01:01'.
* `-w, --weight`: Set weights to determine consensus. Default = 'arcasHLA Polysolver seq2HLA xHLA'.
* `-o, --object`: Which HLA types to return. Options are: 'I', 'II' or 'both'.
* `-d, --sep`: What separator to use between the unique haplotypes. pvactools requires ','.
* `-l, --log`: Boolean if a log file should be created. DEFAULT = False".
* `-f, --logfile`: "Path to directory of log output.

Primarily for use in pvactools bash scripts where comma separated alleles must be provided. To read the json into an bash array use:

```Bash
readarray arcashlatype <<< "$(python haplo_comb.py -j /path/to/file.json -r 2 -o both -s ',')"
```
`$arcashlatype` can then be used as an argument for pvactools

Input of MHC_comb are the native output of haplotype calling tools.

output of arcasHLA in json format:
```
{"A": ["A*01:02:01", "A*02:01:01"],
"B": ["B*18:01:01", "B*15:01:01"],
"C": ["C*05:01:01", "C*07:01:01"],
"DMA": ["DMA*01:01:01"],
"DMB": ["DMB*01:03:01", "DMB*01:01:01"],
"DOA": ["DOA*01:01:02"],
"DOB": ["DOB*01:03", "DOB*01:03:01"],
"DPA1": ["DPA1*01:02:01"],
"DPB1": ["DPB1*01:02:01","DPB1*04:01:01"],
"DQA1": ["DQA1*02:02:01"],
"DQB1": ["DQB1*02:02:01", "DQB1*05:01:01"],
"DRA": ["DRA*02:01:01"],
"DRB1": ["DRB1*01:01:01", "DRB1*03:01:01"],
"DRB3": ["DRB3*02:02:01"]}
```
