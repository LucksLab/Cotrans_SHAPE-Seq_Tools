# Cotrans_SHAPE-Seq_Tools v0.0.2
Tools for Processing Cotranscriptional SHAPE-Seq Datasets
Used in doi:10.1038/nsmb.3316

-------------------------------------------------------------------------------------------------------------

Cotrans_targets.py:
Generates a targets file for spats when analyzing intermediate lengths
of an RNA with an attached adapter for RT priming.

For an RNA of N nt, produces N-y (y=user chosen length) intermediate lengths in a
targets file with the supplied adapter sequence at the 3' end.

Usage:
   `Cotrans_targets.py [options] <RNA name> <RNA sequence> <targets filename(.fa recommended)>`

Options
```
-h, --help                  opens help message
-v, --version               displays version number
-a, --adapter <sequence>    Adapter sequence (5'->3') to add at the 3' end of the RNA w/ barcodes (if applicable)
                            as it appears on the RT primer (after handle) (ex for IDT 2: gtccttggtgcccgagtg; for IDT2_mod: gtccttggtgcccgagtcag)
-n, --adapter-name <string> Name for adapter to include in names in targets file (ex: IDT2)                     
-i, --min-len <N>           Minimum length of RNA from 5' end to make target with (default = 20)
-m, --max-len <N>           Maximum length of RNA from 5' end to make target with (default is entire RNA)
-x, --entire-only           Only produce the target for the max length
```

-------------------------------------------------------------------------------------------------------------

Cotrans_plot_Matlab.py:
Takes a directory containing cotranscriptional SHAPE-Seq reactivities and converts it to a nice matrix form

Usage:                                                                                                          
   `python Cotrans_plot_Matlab.py [options] <reactivities_dir>`

General Options:                                                                                                
```
-h, --help                  Opens help message
-v, --version               Displays version number
-l, --linker <string>       Linker sequence to exclude from RNAs (default is IDT mod: CTGACTCGGGCACCAAGGA)
-r, --recalc_end            Recalculate rhos. Specify the last index of the rhos to be considered. Ex. -1 will drop the last rho value and renormalize based off of reduced length.
-i, --min-len <N>           Minimum length to include (default is first length)
-m, --max-len <N>           Maximum length to include (default is last length)
```

-------------------------------------------------------------------------------------------------------------

Cotrans_matrix_rhos_processing_2D.m and Cotrans_matrix_rhos_processing_3D.m takes tabular formated SHAPE-Seq
 reactivities and alignment numbers and plots then in 2D heatmap or 3D bar chart, respectively.

-------------------------------------------------------------------------------------------------------------

Cotrans_matrix_rhos_processing_differences.m takes two tabular formatted SHAPE-Seq reactivities and plots
 their differences in a 2D heatmap colored by reactivity.
