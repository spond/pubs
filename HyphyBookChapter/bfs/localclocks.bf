SetDialogPrompt("Select a nucleotide data file:"); 

fprintf(stdout,"\n\n LOCAL MOLECULAR CLOCK ANALYSIS: \n");
LikelihoodFunction theLikFun = (myFilter, myTree);