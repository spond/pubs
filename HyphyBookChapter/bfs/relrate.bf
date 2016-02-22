SetDialogPrompt("Select a nucleotide data file:");
lnLA=paramValues[1][0];
dfA=paramValues[1][1];
df0=paramValues[1][1];

LRT=-2*(lnL0-lnLA);
Pvalue=1-CChi2(LRT,dfA-df0);
fprintf(stdout,"\n\nThe statistic ",LRT," has P-value ", Pvalue,"\n\n");
