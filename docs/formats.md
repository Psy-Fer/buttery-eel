# sequencing_summary

NOT finalised. This is just a draft subject to major changes.
Close to ONT's sequencing summary, but documented, with a spec that tells how changes are handled.
A tab-delimited text file with the first row being the header.

Things to think:
1. Strict order or arbitrary order?
2. how are new additions and removal/changes to existing columns are handled?

| Col | Type   | Name                     | Example value | Description |
| --- | ------ | ------------------------ | ----          | ----------- |
| 1   | string | filename_fastq           | reads.fastq | name of the fastq file this read was written to    |
| 2   | string | filename_slow5           | reads.blow5 | name of the slow5 file this read was read from     |
| 3   | string | parent_read_id           |  | the original read_id if the read was split.  If not split, will be the same as the read_id column below.             |
| 4   | string | read_id                  | <sub>f750bda8-0417-429f-b2ec-54c700543ed6
</sub> | the unique read_id. This is a Universally unique identifier (UUID) version 4 and should be unique for any read from any device.   |
| 5   | string | run_id                   | dc60b20f5078b3546ded810fb828b49c438fbd89322
| The unique run ID which will be different for each run                           |
| 6   | int    | channel                  | 1 | The channel number.  See Table 5 of [slow5 specification](https://hasindu2008.github.io/slow5specs/slow5-v1.0.0.pdf) for details. |
| 7   | int    | mux                      | 2 | The MUX setting for the channel when the read began. See Table 5 of [slow5 specification](https://hasindu2008.github.io/slow5specs/slow5-v1.0.0.pdf) for details                                    |
| 8   | int    | minknow_events           | | an internal value                                  |
| 9   | int    | start_time               | | start time of the read                             |
| 10  | int    | duration                 | | time it took from start time to sequence read      |
| 11  | string | passes_filtering         | PASS | TRUE/FALSE for passing the minimum qscore          |
| 12  | ?      | template_start           | | .                              |
| 13  | int    | num_events_template      | | number of events present in read                   |
| 14  | ?      | template_duration        | 16178 | .           |
| 15  | int    | sequence_length_template | 11.885449 | number of bases in the read                        |
| 16  | float  | mean_qscore_template     | 3.239253 | the mean qscore of the read                        |
| 17  | float  | strand_score_template    | 106.933281 | .                                       |
| 18  | float  | median_template          | 19.054329 | median of the template                             |
| 19  | float  | mad_template             | PGXX22394 | median absolute difference for the template        |
| 21  | string | experiment_id            | PGXX22394 | experiment ID if there is one present              |
| 22  | string | sample_id                | HG002 | the sample ID if there is one present              |
| 23  | string | end_reason               | signal_positive| the reason why the read ended                      |


# barcode_summary

| Col | Type   | Name                          | Example value | Description |
| --- | ------ | ----------------------------- |----  | ----------- |
| 1   | string | parent_read_id                | <sub>88be501f-f700-4d9a-8341-13661d579fac</sub> 	| The original read_id if the read was split. If not split, will be the same as the read_id column below.       |
| 2   | string | read_id                       | <sub>88be501f-f700-4d9a-8341-13661d579fac</sub>   | The unique read_id. This is a Universally unique identifier (UUID) version 4 and should be unique for any read from any device.                               |
| 3   | string | barcode_arrangement           | barcode02 								| The barcode name. Usually barcodexx or unclassifed                                |
| 4   | string | barcode_full_arrangement      | NB02_var2							    | ?           										|
| 5   | string | barcode_kit                   | NB									    | Barcode kit the barcode comes from. Likely to be NB or RB for native and rapid               |
| 6   | string | barcode_variant               | var2								    | Which variant of the barcode was detected        |
| 7   | float  | barcode_score                 | 62.16666794 							| the alignment score of the barcode               |
| 8   | string | barcode_front_id              | NB02_FWD							    | the barcode name of the front detected barcode   |
| 9   | float  | barcode_front_score           | 62.16666794						    | the alignment score of the front barcode         |
| 10  | string | barcode_front_refseq          | <sub>ATTGCTAAGGTTAAACAGAC... GACTACAAACGGAATCGACAGCACC</sub> | the front barcode reference sequence             |
| 11  | string | barcode_front_foundseq        | <sub>ACTTGTAGGTGAAACAGACA... TCTAGCAACACTAAGAATGCCC</sub>	   | the sequence found to represent the front barcode|
| 12  | int    | barcode_front_foundseq_length | 42 | the length of the barcode_front_foundseq         |
| 13  | int    | barcode_front_begin_index     | 76 | the start position of the barcode_front_foundseq |
| 14  | string | barcode_rear_id               | NB02_REV | the barcode name of the rear detected barcode    |
| 15  | float  | barcode_rear_score            | 15.33333302 | the alignment score of the rear barcode          |
| 16  | string | barcode_rear_refseq           | <sub>AGGTGCTGTCGATTCCGTTT... GTAGTCGTCTGTTTAACCT</sub> | the rear barcode reference sequence              |
| 17  | string | barcode_rear_foundseq         | <sub>AGTGTGTCATGAGGGCACAA... ACAAATGTCTGTGTAAACAC</sub> | the sequence found to represent the rear barcode |
| 18  | int    | barcode_rear_foundseq_length  | 40 | the length of the barcode_rear_foundseq          |
| 19  | int    | barcode_rear_end_index        | 76 | the start position of the barcode_rear_foundseq  |


# FASTQ header tags
| Type   | Name                          | Description |
| ------ | ----------------------------- | ----------- |
| string | parent_read_id                | the original read_id if the read was split           |
| string | model_version_id              | the basecalling model used           |
| int | mean_qscore                      | the mean qscore of the read           |
| string | barcode                       | the barcode name. only if demultiplexing is enabled            |

# SAM tags

| Type   | Name               | Description |
| ------ | ------------------ | ----------- |
| int (i) | qs                | the mean qscore of the read           |
| int (i) | ns                | number of samples           |
| int (i) | ts                | number of trimmed samples           |
