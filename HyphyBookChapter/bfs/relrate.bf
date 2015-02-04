SetDialogPrompt("Select a nucleotide data file:");DataSet myData = ReadDataFile (PROMPT_FOR_FILE);DataSetFilter myFilter = CreateFilter (myData,1);HarvestFrequencies (obsFreqs, myFilter, 1, 1, 1);F81RateMatrix = {{* ,mu,mu,mu}{mu,* ,mu,mu}{mu,mu,* ,mu}{mu,mu,mu,* }};Model F81 = (F81RateMatrix, obsFreqs);Tree myTree = (a,b,og);fprintf(stdout,"\n Unconstrained analysis:\n\n");LikelihoodFunction theLikFun = (myFilter, myTree, obsFreqs);Optimize (paramValues, theLikFun);fprintf  (stdout, theLikFun);
lnLA=paramValues[1][0];
dfA=paramValues[1][1];fprintf(stdout,"\n\n\n Constrained analysis:\n\n");myTree.a.mu := myTree.b.mu;Optimize (paramValues, theLikFun);fprintf  (stdout, theLikFun);lnL0=paramValues[1][0];
df0=paramValues[1][1];

LRT=-2*(lnL0-lnLA);
Pvalue=1-CChi2(LRT,dfA-df0);
fprintf(stdout,"\n\nThe statistic ",LRT," has P-value ", Pvalue,"\n\n");
