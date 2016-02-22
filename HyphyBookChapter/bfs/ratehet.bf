REPLACE_TREE_STRUCTURE = 1;

/* this will display the category info, along with sample mean and variance */

SetDialogPrompt("Select a nucleotide data file:"); 

fprintf(stdout,"\n\n  Homogeneous Rates Analysis: \n");


/* Variable Rates Model 1: 4 equiprobable rate classes, including an invariant class. */

fprintf(stdout,"\n\n  Variable Rates Model 1: 4 equiprobable rate classes, including an invariant class. \n");
category rateCat = (4, EQUAL, MEAN, , {{0}{1}{2}{4}}, 0, 4);

F81VarRateMatrix = {{*,rateCat*m,rateCat*m,rateCat*m}{rateCat*m,*,rateCat*m,rateCat*m}
		    {rateCat*m,rateCat*m,*,rateCat*m}{rateCat*m,rateCat*m,rateCat*m,*}};

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);
GetInformation (catInfo, rateCat);


/* Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. */

fprintf(stdout,"\n\n  Variable Rates Model 2: Discrete Gamma with 4 equiprobable rate classes. \n");

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (4, EQUAL, MEAN, GammaDist(_x_,alpha,alpha), CGammaDist(_x_,alpha,alpha),0,1e25,CGammaDist(_x_,alpha+1,alpha));

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);
GetInformation (catInfo, rateCat);



/* Variable Rates Model 3: Discrete Exponential with 4 equiprobable rate classes. */

fprintf(stdout,"\n\n  Variable Rates Model 3: Discrete Exponential with 4 equiprobable rate classes. \n");

global alpha=2;
alpha:>0.01; alpha:<100;

category rateCat = (4, EQUAL, MEAN, alpha*Exp(-alpha*_x_), 1-Exp(-alpha*_x_), 0, 1e25, -_x_*Exp(-alpha*_x) + (1-Exp(-alpha*_x_))/alpha);

Model F81Var = (F81VarRateMatrix, obsFreqs);

Tree myTree = myTopology;
LikelihoodFunction theLikFun = (myFilter, myTree);
GetInformation (catInfo, rateCat);








