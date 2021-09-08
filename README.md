# MHC_comb
Combine the output of one or several HLA haplotype calling tools to a concencus HLA typing.

Supports: 
* arcasHLA (https://github.com/RabadanLab/arcasHLA)
* xHLA  (https://github.com/humanlongevity/HLA)
* seq2HLA (https://github.com/TRON-Bioinformatics/seq2HLA)
* POLYSOLVER (https://software.broadinstitute.org/cancer/cga/polysolver)


MHC_comb takes the output from in silico Haplotype tools as input. It determines a concencus HLA type from these and returns prints every biologically possible combination of HLA class I, II or both alleles.

Finally an output file compatible with LOHHLA can be written.

## Dependencies:
* Itertools
* functools
* collections.abc
* collections
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
* `-l, --log`: Boolean if a log file should be created. DEFAULT = False.
* `-f, --logfile`: Path to directory of log output.
* `--lohhla`: Boolean if a LOHHLA compatible .txt file should be created. DEFAULT = False.

To acquire an array of HLA alleles, like required in pvactools bash scripts (comma separated alleles must be provided), use:

```Bash
readarray arcashlatype <<< "$(python haplo_comb.py -a /path/to/file.json -r 2 -o both -d ',')"
```
`$arcashlatype` can then be used as an argument for pvactools

You can combine outputs of various tools by the associated arguments. Concencsus on haplotypes will initially be based on the number of occurances accross the tools. In the case of a tie, a tie breaker score gets calculated based on the weights assigned to each tool. With the highest weight assigned to the first place.


```Bash
readarray arcashlatype <<< "$(python haplo_comb.py -a /path/to/file.json -x /path/to/xhla.tsv -w 'arcasHLA xHLA' -r 2 -o both -d ',')"
```



output of LOHHLA compatible format:
|<!-- -->|    <!-- -->    |   <!-- -->      |
| ------ | -------------- | --------------- |
| HLA-A  | hla_a_01_02_01 | hla_a_01_02_01  |
| HLA-B  | hla_b_18_01_01 | hla_b_15_01_01  |
| HLA-C  | hla_c_05_01_01 | hla_c_05_01_01  |




