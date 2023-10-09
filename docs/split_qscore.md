# Split qscore script


This script allows you to split a bulk fastq or sam file by mean quality score.

It will first try to read the `mean_qscore=` tag in the fastq read header or the `qs:i:` sam tag.
If this isn't present, it will calculate the mean_qscore from the read quality string.

Note, that the calculated values won't match the values output by basecallers and present in the tags, as there is some slight scaling that happens that is not available to this script to correct for, however the values are still very close.

The script uses the following method to calculate the mean_qscore of a read

```python
def calculate_qscore(qstring):
    '''
    calculate a qscore from a qstring
    '''
    qs = (np.array(qstring, 'c').view(np.uint8) - 33)
    mean_err = np.exp(qs * (-np.log(10) / 10.)).mean()
    score = -10 * np.log10(max(mean_err, 1e-4))
    return score
```

## Example

```
python3 scripts/split_qscore.py -q 9 input.fastq output/path/
```

## Usage

```
usage: split_qscore.py [-h] [-p PREFIX] [-q QSCORE] [-f {fastq,sam}] input output

split a fastq or sam file by qscore value into pass/fail files

positional arguments:
  input                 fastq, sam file or - for STDIN/pipe
  output                output path to write pass/fail files

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        filename prefix to give pass/fail output files, eg -p example would give example.pass.fastq & example.fail.fastq (default: reads)
  -q QSCORE, --qscore QSCORE
                        qscore to split on: pass=>qscore (default: 9)
  -f {fastq,sam}, --format {fastq,sam}
                        if using - for stdin, give format as fastq or sam (default: None)
```