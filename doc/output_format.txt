
Output
------

For default/ranked/weighted analysis, the first six columns are always the same:

1. kmer
2. position: first nucleotide of the kmer in the sequence, starts with 0.
3. shift
4. statistics (positive for enrichment and negative for depletion)
5. p.value (-log10)
6. Bonferoni corrected p.value (-log10)

And the last column is always FDR. Columns in the middle varies.

Default:
7. observed fraction of input sequences have this kmer at this position allowing this amount of shift
8. expected fraction of background sequences have this kmer at this position allowing this amount of shift
9. ratio: column 7 divided by column 8
10. column 7 divided by the average of other positions 


Ranked:

No other information is provided. The statistics and p-value in column 4 and 5 are based on normal approximation. See: Carine A. Bellera et al, Normal Approximations to the Distributions of the Wilcoxon Statistics: Accurate to What N? Graphical Insights, Journal of Statistics Education, Volume 18, Number 2, (2010) http://www.amstat.org/publications/jse/v18n2/bellera.pdf  

Weighted:
Additional columns are:
7. sample size for sample 1, i.e. the number of sequences have this kmer at this position when allowing this amount of shift. 
8. mean of weight for sample 1 sequences.
9. standard deviation of weight for sample 1 sequences.
10. sample size for sample 2, i.e. the number of sequences do not have this kmer at this position when allowing this amount of shift. 
11. mean of weight for sample 2 sequences.
12. standard deviation of weight for sample 2 sequences.
13. FDR