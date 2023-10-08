# sequencing_summary

NOT finalised. This is just a draft subject to major changes.
Close to ONT's sequencing summary, but documented, with a spec that tells how changes are handled.
A tab-delimited text file with the first row being the header. 

Things to think:
1. Strict order or arbitary order?
2. how are new additions and removal/changes to existing columns are handled?

| Col | Type   | Name                     | Description |
| --- | ------ | ------------------------ | ----------- |
| 1   | string | filename_fastq           | .           |
| 2   | string | filename_slow5           | .           |
| 3   | string | parent_read_id           | .           |
| 4   | string | read_id                  | .           |
| 5   | string | run_id                   | .           |
| 6   | int    | channel                  | .           |
| 7   | int    | mux                      | .           |
| 8   | int    | minknow_events           | .           |
| 9   | int    | start_time               | .           |
| 10  | int    | duration                 | .           |
| 11  | string | passes_filtering         | .           |
| 12  | ?      | template_start           | .           |
| 13  | int    | num_events_template      | .           |
| 14  | ?      | template_duration        | .           |
| 15  | int    | sequence_length_template | .           |
| 16  | float  | mean_qscore_template     | .           |
| 17  | float  | strand_score_template    | .           |
| 18  | float  | median_template          | .           |
| 19  | float  | mad_template             | .           |
| 20  | string | pore_type                | .           |
| 21  | string | experiment_id            | .           |
| 22  | string | sample_id                | .           |
| 23  | string | end_reason               | .           |


# barcode_summary

| Col | Type   | Name                          | Description |
| --- | ------ | ----------------------------- | ----------- |
| 1   | string | parent_read_id                | .           |
| 2   | string | read_id                       | .           |
| 3   | string | barcode_arrangement           | .           |
| 4   | string | barcode_full_arrangement      | .           |
| 5   | string | barcode_kit                   | .           |
| 6   | string | barcode_variant               | .           |
| 7   | float  | barcode_score                 | .           |
| 8   | string | barcode_front_id              | .           |
| 9   | float  | barcode_front_score           | .           |
| 10  | string | barcode_front_refseq          | .           |
| 11  | string | barcode_front_foundseq        | .           |
| 12  | int    | barcode_front_foundseq_length | .           |
| 13  | int    | barcode_front_begin_index     | .           |
| 14  | string | barcode_rear_id               | .           |
| 15  | float  | barcode_rear_score            | .           |
| 16  | string | barcode_rear_refseq           | .           |
| 17  | string | barcode_rear_foundseq         | .           |
| 18  | int    | barcode_rear_foundseq_length  | .           |
| 19  | int    | barcode_rear_end_index        | .           |


# FASTQ header tags
| Type   | Name                          | Description |
| ------ | ----------------------------- | ----------- |
| string | parent_read_id                | .           |
| string | model_version_id              | .           |
| int | mean_qscore                      | .           |
| string | barcode                       | the barcode name. only if demultiplexing is enabled            |  

# SAM tags

| Type   | Name               | Description |
| ------ | ------------------ | ----------- |
| int (i) | qs                | .           |
| int (i) | ns                | .           |
| int (i) | ts                | .           |
