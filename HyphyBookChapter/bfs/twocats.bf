REPLACE_TREE_STRUCTURE = 1;/* This statement tells HyPhy to completely erase all parameters associated with a Tree when the Treeis being rewritten. Default behavior is to keep parameters from old analyses, including the onesno longer included in rate models.*/

/* this will display the category info, along with sample mean and variance */#include "displayFunction.bf";

SetDialogPrompt("Select a nucleotide data file:"); DataSet myData  = ReadDataFile(PROMPT_FOR_FILE);DataSetFilter myFilter = CreateFilter(myData,1);HarvestFrequencies(obsFreqs,myFilter,1,1,1);myTopology = "(((a,b),(c,(d,e))),f)"; 

global TV=0.5;/* Discrete Gamma with 16 equiprobable rate classes. */

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (16, EQUAL, MEAN, GammaDist(_x_,alpha,alpha), CGammaDist(_x_,alpha,alpha),0,1e25,CGammaDist(_x_,alpha+1,alpha));

HKY85VarRateMatrix = 		{{*,			t*rateCat*TV,		t*rateCat,		t*rateCat*TV}		 {t*rateCat*TV,		*,			t*rateCat*TV,		t*rateCat}		 {t*rateCat,		t*rateCat*TV,		*,			t*rateCat*TV}		 {t*rateCat*TV,		t*rateCat,		t*rateCat*TV,		*}};


fprintf(stdout,"\n\n  Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. \n");


Model HKY85Var	 = (HKY85VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);
GetInformation (catInfo, rateCat);catInfo = echoCatVar (catInfo);



/* Discrete Gamma with heterogeneous transition and transversion rates, each with 4 equiprobable rate classes. */

global alphaS=2;
alphaS:>0.01; alphaS:<100;

global alphaV=2;
alphaV:>0.01; alphaV:<100;

global beta=2;
beta:>0.01; beta:<100;

category catTS = (4, EQUAL, MEAN, GammaDist(_x_,alphaS,alphaS), CGammaDist(_x_,alphaS,alphaS),0,1e25,CGammaDist(_x_,alphaS+1,alphaS));


category catTV = (4, EQUAL, MEAN, GammaDist(_x_,alphaV,beta), CGammaDist(_x_,alphaV,beta),0,1e25,CGammaDist(_x_,alphaV+1,beta)*alphaV/beta);


HKY85TwoVarRateMatrix = 		{{*,			t*catTV*TV,		t*catTS,		t*catVarTV*TV}		 {t*catTV*TV,		*,			t*catTV*TV,		t*catTS}		 {t*catTS,		t*catTV*TV,		*,			t*catTV*TV}		 {t*catTV*TV,		t*catTS,		t*catTV*TV,		*}};

fprintf(stdout,"\n\n Discrete Gamma with heterogeneous transition and transversion rates, each with 4 equiprobable rate classes.  \n");


Model HKY85Var	 = (HKY85TwoVarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);

GetInformation (catInfo, catTS);catInfo = echoCatVar (catInfo);
GetInformation (catInfo, catTV);catInfo = echoCatVar (catInfo);








