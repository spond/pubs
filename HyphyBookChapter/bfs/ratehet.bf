REPLACE_TREE_STRUCTURE = 1;/* This statement tells HyPhy to completely erase all parameters associated with a Tree when the Treeis being rewritten. Default behavior is to keep parameters from old analyses, including the onesno longer included in rate models.*/

/* this will display the category info, along with sample mean and variance */#include "displayFunction.bf";

SetDialogPrompt("Select a nucleotide data file:"); DataSet myData  = ReadDataFile(PROMPT_FOR_FILE);DataSetFilter myFilter = CreateFilter(myData,1);HarvestFrequencies(obsFreqs,myFilter,1,1,1);myTopology = "(((a,b),(c,(d,e))),f)"; F81RateMatrix  = {{*,m,m,m}{m,*,m,m}{m,m,*,m}{m,m,m,*}};Model F81 = (F81RateMatrix, obsFreqs);

fprintf(stdout,"\n\n  Homogeneous Rates Analysis: \n");Tree myTree = myTopology;LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);


/* Variable Rates Model 1: 4 equiprobable rate classes, including an invariant class. */

fprintf(stdout,"\n\n  Variable Rates Model 1: 4 equiprobable rate classes, including an invariant class. \n");
category rateCat = (4, EQUAL, MEAN, , {{0}{1}{2}{4}}, 0, 4);

F81VarRateMatrix = {{*,rateCat*m,rateCat*m,rateCat*m}{rateCat*m,*,rateCat*m,rateCat*m}
		    {rateCat*m,rateCat*m,*,rateCat*m}{rateCat*m,rateCat*m,rateCat*m,*}};

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);
GetInformation (catInfo, rateCat);catInfo = echoCatVar (catInfo);


/* Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. */

fprintf(stdout,"\n\n  Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. \n");

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (4, EQUAL, MEAN, GammaDist(_x_,alpha,alpha), CGammaDist(_x_,alpha,alpha),0,1e25,CGammaDist(_x_,alpha+1,alpha));

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);
GetInformation (catInfo, rateCat);catInfo = echoCatVar (catInfo);



/* Variable Rates Model 3: Discrete Exponential with 4 equiprobable rate classes. */

fprintf(stdout,"\n\n  Variable Rates Model 3: Discrete Exponential with 4 equiprobable rate classes. \n");

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (4, EQUAL, MEAN, alpha*Exp(-alpha*_x_), 1-Exp(-alpha*_x_), 0, 1e25, -_x_*Exp(-alpha*_x) + (1-Exp(-alpha*_x_))/alpha);

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);Optimize(results,theLikFun);fprintf(stdout,theLikFun);
GetInformation (catInfo, rateCat);catInfo = echoCatVar (catInfo);









