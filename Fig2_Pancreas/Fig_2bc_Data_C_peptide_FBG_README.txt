Fig_2bc_Data_C_peptide_FBG.csv - FBG column removed (2026-09-06)

This file originally also carried an FBG column. That column was found
to be misaligned against animal identity: for CTRL it matched the
verified per-animal FBG record (Fig_2a_Data_FBG.csv), but for every
IL10-KO and DSS-IL10-KO animal it did not - the FBG value stored under
a given Replicate label belonged to a different animal than the
C_peptide value stored under that same label.

Verification: Fig_2a_Data_FBG.csv is confirmed correct because it
reproduces the published Fig_2a_Summary_FBG.csv exactly at every group
and week. The removed FBG column reproduced no published summary.

Root cause: traced to Experiments/HOMA/R/HOMA_data.csv, an earlier,
now-superseded HOMA-index analysis that predates this manuscript. Even
that original analysis's own FBG boxplot script read FBG from this same
Replicate-labelled, misaligned structure rather than from an
ID-verified source. Fig_2bc_Data_C_peptide_FBG.csv inherited the FBG
values from HOMA_data.csv without correction; the C_peptide column is
unaffected (it reproduces Fig_2b_Summary_C_peptide.csv exactly) and is
still the authoritative C-peptide source.

Anything in this project needing FBG for these animals should use
Fig_2a_Data_FBG.csv (ID-indexed) instead, joining to C_peptide/other
per-replicate data by Group + Replicate. Fig_2c_Boxplot_FCP_FBG_Ratio.r
and Supplementary_Table_3_Mixed_Effects_Models.R were both updated to
do this; both previously used the removed column and their outputs have
been regenerated accordingly.
