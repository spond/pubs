REPLACE_TREE_STRUCTURE = 1;

/* this will display the category info, along with sample mean and variance */

SetDialogPrompt("Select a nucleotide data file:"); 

global TV=0.5;

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (16, EQUAL, MEAN, GammaDist(_x_,alpha,alpha), CGammaDist(_x_,alpha,alpha),0,1e25,CGammaDist(_x_,alpha+1,alpha));

HKY85VarRateMatrix = 


fprintf(stdout,"\n\n  Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. \n");


Model HKY85Var	 = (HKY85VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);
GetInformation (catInfo, rateCat);



/* Discrete Gamma with heterogeneous transition and transversion rates, each with 4 equiprobable rate classes. */

global alphaS=2;
alphaS:>0.01; alphaS:<100;

global alphaV=2;
alphaV:>0.01; alphaV:<100;

global beta=2;
beta:>0.01; beta:<100;

category catTS = (4, EQUAL, MEAN, GammaDist(_x_,alphaS,alphaS), CGammaDist(_x_,alphaS,alphaS),0,1e25,CGammaDist(_x_,alphaS+1,alphaS));


category catTV = (4, EQUAL, MEAN, GammaDist(_x_,alphaV,beta), CGammaDist(_x_,alphaV,beta),0,1e25,CGammaDist(_x_,alphaV+1,beta)*alphaV/beta);


HKY85TwoVarRateMatrix = 

fprintf(stdout,"\n\n Discrete Gamma with heterogeneous transition and transversion rates, each with 4 equiprobable rate classes.  \n");


Model HKY85Var	 = (HKY85TwoVarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);

GetInformation (catInfo, catTS);
GetInformation (catInfo, catTV);







