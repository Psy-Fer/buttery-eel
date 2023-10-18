# The mysterious default parameters in Guppy

##  qscore cutoff 

The qscore cutoff (that determine pass/fail) in Guppy keeps on changing and there is no easy way to get this. So this is an attempt to document based on what we find.


| Guppy version | Model                         | default qscore cut off |
|--------       |------                         | ---                    | 
| 6.5.7         | dna_r10.4.1_e8.2_400bps_5khz_hac_prom   |      9                   |
|               |                               |                      | 
| 6.4.2         | dna_r9.4.1_450bps_fast_prom   |    8                   |
| 6.4.2         | dna_r9.4.1_450bps_hac_prom    |    9                   |
| 6.4.2         | dna_r9.4.1_450bps_sup_prom    |    10                  |
| 6.4.2         | dna_r10.4.1_e8.2_400bps_fast_prom   |      8                   |
| 6.4.2         | dna_r10.4.1_e8.2_400bps_hac_prom   |      9                   |
| 6.4.2         | dna_r10.4.1_e8.2_400bps_sup   |      10                  |
|               |                               |                      | 
| 4.0.11        | dna_r9.4.1_450bps_hac_prom    | 7                      |

- Guppy 6.4.2 [dna_r9.4.1_450bps_fast_prom] (qscore cut off - default:8)
- Guppy 6.4.2 [dna_r9.4.1_450bps_hac_prom] (qscore cut off - default:9)
- Guppy 6.4.2 [dna_r9.4.1_450bps_sup_prom] (qscore cut off - default:10)
- Guppy 6.4.2 [dna_r10.4.1_e8.2_400bps_fast_prom] (qscore cut off - default 8) 
- Guppy 6.4.2 [dna_r10.4.1_e8.2_400bps_hac_prom] (qscore cut off - default 9)
- Guppy 6.4.2 [dna_r10.4.1_e8.2_400bps_sup] (qscore cut off - default 10)

