# sequencing_summary

NOT finalised. This is just a draft subject to major changes.
Close to ONT's sequencing summary, but documented, with a spec that tells how changes are handled.
A tab-delimited text file with the first row being the header. 

Things to think:
1. Strict order or arbitary order?
2. how are new additions and removal/changes to existing columns are handled?

| Col | Type   | Name                     | Description |
| --- | ------ | ------------------------ | ----------- |
| 1   | string | filename_fastq           | name of the fastq file this read was written to    |
| 2   | string | filename_slow5           | name of the slow5 file this read was read from     |
| 3   | string | parent_read_id           | the original read_id if the read was split         |
| 4   | string | read_id                  | the unique read_id                                 |
| 5   | string | run_id                   | the unique id of the run                           |
| 6   | int    | channel                  | the well number the molecule was sequenced through |
| 7   | int    | mux                      | the mux number                                     |
| 8   | int    | minknow_events           | an internal value                                  |
| 9   | int    | start_time               | start time of the read                             |
| 10  | int    | duration                 | time it took from start time to sequence read      |
| 11  | string | passes_filtering         | TRUE/FALSE for passing the minimum qscore          |
| 12  | ?      | template_start           | .                              |
| 13  | int    | num_events_template      | number of events present in read                   |
| 14  | ?      | template_duration        | .           |
| 15  | int    | sequence_length_template | number of bases in the read                        |
| 16  | float  | mean_qscore_template     | the mean qscore of the read                        |
| 17  | float  | strand_score_template    | .                                       |
| 18  | float  | median_template          | median of the template                             |
| 19  | float  | mad_template             | median absolute difference for the template        |
| 20  | string | experiment_id            | experiment ID if there is one present              |
| 21  | string | sample_id                | the sample ID if there is one present              |
| 22  | string | end_reason               | the reason why the read ended                      |


# barcode_summary

| Col | Type   | Name                          | Description |
| --- | ------ | ----------------------------- | ----------- |
| 1   | string | parent_read_id                | the original read_id if the read was split       |
| 2   | string | read_id                       | the unique read_id                               |
| 3   | string | barcode_arrangement           | the barcode name           |
| 4   | string | barcode_full_arrangement      |            |
| 5   | string | barcode_kit                   | barcode kit the barcode comes from               |
| 6   | string | barcode_variant               | which variant of the barcode was detected        |
| 7   | float  | barcode_score                 | the alignment score of the barcode               |
| 8   | string | barcode_front_id              | the barcode name of the front detected barcode   |
| 9   | float  | barcode_front_score           | the alignment score of the front barcode         |
| 10  | string | barcode_front_refseq          | the front barcode reference sequence             |
| 11  | string | barcode_front_foundseq        | the sequence found to represent the front barcode|
| 12  | int    | barcode_front_foundseq_length | the length of the barcode_front_foundseq         |
| 13  | int    | barcode_front_begin_index     | the start position of the barcode_front_foundseq |
| 14  | string | barcode_rear_id               | the barcode name of the rear detected barcode    |
| 15  | float  | barcode_rear_score            | the alignment score of the rear barcode          |
| 16  | string | barcode_rear_refseq           | the rear barcode reference sequence              |
| 17  | string | barcode_rear_foundseq         | the sequence found to represent the rear barcode |
| 18  | int    | barcode_rear_foundseq_length  | the length of the barcode_rear_foundseq          |
| 19  | int    | barcode_rear_end_index        | the start position of the barcode_rear_foundseq  |


# FASTQ header tags
| Type   | Name                          | Description |
| ------ | ----------------------------- | ----------- |
| string | parent_read_id                | the original read_id if the read was split           |
| string | model_version_id              | the basecalling model used           |
| int | mean_qscore                      | the mean qscore of the read           |
| string | barcode                       | the barcode name. only if demultiplexing is enabled            |  

# SAM tags

| Type       | Name          | Description |
| ---------- | ------------- | ----------- |
| int (i)    | qs            | the mean qscore of the read           |
| int (i)    | ns            | number of samples           |
| int (i)    | ts            | number of trimmed samples           |
| string (Z) | pi            | the parent id/original id if the read was split  |
