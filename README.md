# MHC_comb
Takes .json files as input and prints every possible combination of HLA class II combinations to be picked up in bash script 

## Dependencies:
* Itertools
* Combinations
* json

## Usage:
For a short explanation per argument type `haplo_comb.py -h`

* `-j, --json`: Path to the arasHLA .json output.
* `-r, --res`: The resolution of the HLA alleles that get returned.
* `-o, --object`: Which HLA types to return. Options are: 'I', 'II' or 'both'
* `-s, --sep`: What separator to use between the unique haplotypes. pvactools requires ','.

Primarily for use in pvactools bash scripts where comma separated alleles must be provided. To read the json into an bash array  use:

```Bash
readarray arcashlatype <<< "$(python haplo_comb.py -j /path/to/file.json -r 2 -o both -s ',')"
```
`$arcashlatype` can then be provided to pvactools

At the moment the input json file requires to be structred like the arcasHLA output:

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

the output will look like:

echo $arcashlatype
