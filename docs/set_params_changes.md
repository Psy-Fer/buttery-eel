 The following options were either removed:

 ```
 |      * post_out (bool):                      Flag indicating whether to return full posterior state data. Default is False.
 |      * lamp_kit (str):                       String naming LAMP barcode kit to use. Default is to not do LAMP barcoding.
 |      * detect_mid_strand_barcodes (bool):    Flag indicating that read will be marked as unclassified if barcodes appear
 |                                              within the strand itself. Default is False.
 |      * detect_mid_strand_adapter (bool):     Flag indicating that read will be marked as unclassified if the adapter
 |                                              sequence appears within the strand itself. Default is False.
 |      * min_score_barcode_front (float):      Minimum score for a front barcode to be classified. Default is 60.
 |      * min_score_barcode_rear (float):       Minimum score for a rear barcode to be classified. Default is to use the front
 |                                              minimum.
 |      * min_score_barcode_mid (float):        Minimum score for mid barcodes to be detected. Default is 50.
 |      * detect_adapter (bool):                Enable detection of adapters at the front and rear of the sequence.
 |      * detect_primer (bool):                 Enable detection of primers at the front and rear of the sequence.
 |      * min_score_adapter (float):            Minimum score for a front or rear adapter to be classified. Default is 60.
 |      * min_score_primer (float):             Minimum score for a front or rear primer to be classified. Default is 60.
 |      * align_type (str):                     Type of alignment requested. Valid values are "auto", "coarse", and "full".
 |                                              Default is "auto".
 |      * min_score_lamp (float):               Minimum score for a LAMP barcode to be detected. Default is 80.
 |      * min_score_lamp_mask (float):          Minimum score for a LAMP barcode mask to be detected. Default is 50.
 |      * min_score_lamp_target (float):        Minimum score for a LAMP target to be detected. Default is 50.
 |      * additional_lamp_context_bases (int):  Number of additional bases from context to include in lamp barcode comparison.
 |                                              Default is 2.
```

 The follow were changed:

old
```
  |      * move_and_trace_enabled (bool):        Flag indicating whether to return trace and move data. Default is True.
```
new
```
 |      * move_enabled (bool):                  Flag indicating whether to return move data. Default is False.
```