RequireVersion ("2.1");  
_simulation_replicates = 100; 	
   	
BASE_DIRECTORY = PATH_TO_CURRENT_BF;   	

range_lambda_k={};
iter_range=-1;
for (lam=0; lam<=5; lam+=1) {
    lambda=2^(-lam);
    for (K=0.25; K<=3.5; K+=0.25) {
        iter_range += 1;
        range_lambda_k+{};
        range_lambda_k[iter_range]+lambda;
        range_lambda_k[iter_range]+K;	
    }
        
}
for (lam=1; lam<=5; lam+=1) {	
    lambda=2^lam;
    for (K=0.25; K<=3.5; K+=0.25) {
        iter_rang += 1;
        range_lambda_k+{};
        range_lambda_k[iter_range]+lambda;
        range_lambda_k[iter_range]+K;	
    }
        
}

	
likelihood_file_name="`BASE_DIRECTORY`/Results/models.txt";
fprintf (likelihood_file_name,CLEAR_FILE,"Network,SettingsID,lambda,K,mu,P,logL\n");

run = 0;
while (1)	 {
	
	/*--------------------read inputs------*/
	
	run += 1;
	
	UseModel (USE_NO_MODEL);

    ACCEPT_ROOTED_TREES = 1;
    
    _tree_file = "`BASE_DIRECTORY`/Inputs/AN_"+run+".txt";
    
    if ((!_tree_file) == 0) {
        break;
    }
    
    fscanf (_tree_file,"Tree",tree1);
    _tree_file = "`BASE_DIRECTORY`/Inputs/PL_"+run+".txt";
    fscanf (_tree_file,"Tree",tree2);

    _interactions = "`BASE_DIRECTORY`/Inputs/interaction_"+run+".txt";
    fscanf(_interactions,"Lines",line_element);


    interactionsDict = {};	
    for (ind_line_element=0; ind_line_element<Columns (line_element); ind_line_element += 3) {
        interactionsDict[ ((line_element[ind_line_element])&&1) + "," + ((line_element[ind_line_element+1])&&1) ] = 0+(line_element[ind_line_element+2]);
    }

    /* convert both trees to post-order traversal annotated node lists */
    

    tree1AVL = (tree1)^0;
    tree2AVL = (tree2)^0;
    /*-------------------------------------------------------------------------*/



    /*-------------------------------------------------------------------------------------------
    /*first, compute the max depth of the existing nodes in the trees*/
    
    max1=0;
    for (i=1; i<Abs(tree1AVL); i=i+1)
    {
        if ((tree1AVL[i])["Depth"]>=max1)
        {
            max1=(tree1AVL[i])["Depth"];
        }
    }
    max2=0;
    for (i=1; i<Abs(tree2AVL); i=i+1)
    {
        if ((tree2AVL[i])["Depth"]>=max2)
        {
            max2=(tree2AVL[i])["Depth"];
        }
    }

    /*define which tree to scale (the one with smaller depth)*/
    
    if (max1>max2)
    {
        K2=K;
        K1=1;
    }
    else
    {
        K1=K;
        K2=1;
    }

    /*separate the leaves and the internal node of the 2 trees*/

    leaves1={};
    not_leaves1={};
    count1=0;
    count1_not_leaves=0;
    for (node=1; node<Abs(tree1AVL);node=node+1)
    {	if (Abs((tree1AVL[node])["Children"])==0)
        {	

            leaves1[count1]=node;
            count1=count1+1;
        }
        else
        {
            not_leaves1[count1_not_leaves]=node;
            count1_not_leaves=count1_not_leaves+1;

        }
        if ((tree1AVL[node])["Depth"]==0)
            {
                (tree1AVL[node])["Length"]=0;
            }

    }

    leaves2={};
    not_leaves2={};
    count2=0;
    count2_not_leaves=0;
    for (node=1; node<Abs(tree2AVL);node=node+1)
    {	if (Abs((tree2AVL[node])["Children"])==0)
        {	

            leaves2[count2]=node;
            count2=count2+1;
        }
        else
        {
            not_leaves2[count2_not_leaves]=node;
            count2_not_leaves=count2_not_leaves+1;
        }

        if ((tree2AVL[node])["Depth"]==0)
            {
                (tree2AVL[node])["Length"]=0;
            }

    }
	
	for (it_range=0; it_range<Abs(range_lambda_k);it_range=it_range+1) {
		lambda=(range_lambda_k[it_range])[0];
		K=(range_lambda_k[it_range])[1];
			
		fprintf (stdout, "\n");	
		reportProgress("NETWORK: " + run + ",lambda: " + lambda + ", K: " + K);


		/*----------------------------------------------------------*/

		/*--------assign the branch lengths to animal tree---------*/

		/*when the branch is an internal branch*/
		for (node=0; node<Abs(not_leaves1); node=node+1)
		{
			l=(tree1AVL[(not_leaves1[node])])["Depth"];
			L=max1;
			kk=K1*((((l/L)^lambda)*(L/l)));
			((tree1AVL[(not_leaves1[node])])["Length"])=kk;
		}
		/*When the branch is giving a leaf with max depth*/
		for (node=0; node<Abs(leaves1); node=node+1)
		{
			if ((tree1AVL[(leaves1[node])])["Depth"]==max1)
			{
				kk=K1*(((max1/L)^lambda));
				((tree1AVL[(leaves1[node])])["Length"])=kk;
				leaf=leaves1[node];
			}
		}
		/*Compute the length from the leaves to the root*/
		P_1=0;
		N_1=leaf;
		while((tree1AVL[N_1])["Depth"]!=0)
		{
			P_1=P_1+(tree1AVL[N_1])["Length"];
			N_1=(tree1AVL[N_1])["Parent"];
		}
		/*Assign bl for the remaining of the branches*/
		for (node=0; node<Abs(leaves1); node=node+1)
		{
			if ((tree1AVL[(leaves1[node])])["Depth"]!=max1)
			{
				p_1=0;
				n_1=(tree1AVL[(leaves1[node])])["Parent"];
				while((tree1AVL[n_1])["Depth"]!=0)
				{
					p_1=p_1+(tree1AVL[n_1])["Length"];
					n_1=(tree1AVL[n_1])["Parent"];
				}
				kk=P_1-p_1;
				((tree1AVL[(leaves1[node])])["Length"])=kk;

			} 
		}

		/*check that all the leaves are at the same point in time*/
		for (node=0; node<Abs(leaves1); node=node+1)
		{
			age=0;
			nod=leaves1[node];
			while((tree1AVL[nod])["Depth"]!=0)
			{
				age=age+(tree1AVL[nod])["Length"];
				nod=(tree1AVL[nod])["Parent"];
			}
			/*fprintf (stdout, age,"\n");*/


		}



		/*------------------------*/




		/*---------assign branch length to plant tree---------*/

		/*when the branch is an internal branch*/
		for (node=0; node<Abs(not_leaves2); node=node+1)
		{
			l=(tree2AVL[(not_leaves2[node])])["Depth"];
			L=max2;
			kk=K2*(((l/L)^lambda)*(L/l));
			((tree2AVL[(not_leaves2[node])])["Length"])=kk;
		}
		/*When the branch is giving a leaf with max depth*/
		for (node=0; node<Abs(leaves2); node=node+1)
		{
			if ((tree2AVL[(leaves2[node])])["Depth"]==max2)
			{
				kk=K2*((max2/L)^lambda);
				((tree2AVL[(leaves2[node])])["Length"])=kk;
				leaf=leaves2[node];
			}
		}
		/*Compute the length from the leaves to the root*/
		P_2=0;
		N_2=leaf;
		while((tree2AVL[N_2])["Depth"]!=0)
		{
			P_2=P_2+(tree2AVL[N_2])["Length"];
			N_2=(tree2AVL[N_2])["Parent"];
		}
		/*Assign bl for the remaining of the branches*/
		for (node=0; node<Abs(leaves2); node=node+1)
		{
			if ((tree2AVL[(leaves2[node])])["Depth"]!=max2)
			{
				p_2=0;
				n_2=(tree2AVL[(leaves2[node])])["Parent"];
				while((tree2AVL[n_2])["Depth"]!=0)
				{
					p_2=p_2+(tree2AVL[n_2])["Length"];
					n_2=(tree2AVL[n_2])["Parent"];
				}
				kk=P_2-p_2;
				((tree2AVL[(leaves2[node])])["Length"])=kk;

			} 
		}

		/*check that all the leaves are at the same point in time*/
		for (node=0; node<Abs(leaves2); node=node+1)
		{
			age=0;
			nod=leaves2[node];
			while((tree2AVL[nod])["Depth"]!=0)
			{
				age=age+(tree2AVL[nod])["Length"];
				nod=(tree2AVL[nod])["Parent"];
			}

			/*fprintf (stdout, "age",age,"\n");*/


		}

		/*-------------------------------------------------------------------------------------------*/
		
		/* compute the ages of internal nodes (backwards from the present) assumes molecular clock*/
		

		size1 = Abs(tree1AVL) - 1;
		size2 = Abs(tree2AVL) - 1;

		A={size1+size2,2};

		nodeAges = A["_MATRIX_ELEMENT_ROW_"];


		sequenceData = "";

		/*labeling each node*/
		for (k = 1; k <= size1; k = k+1)
		{
			childrenCount = Abs((tree1AVL[k])["Children"]);
			if (childrenCount)
			{
				(tree1AVL[k])["Label"] = "";
				nodeAges [k-1][0] = -(tree1AVL[(((tree1AVL[k])["Children"])[0])])["Length"];
				(tree1AVL[k])["Length"] = Max(0,(tree1AVL[k])["Length"]) + (tree1AVL[(((tree1AVL[k])["Children"])[0])])["Length"];
				for (k2 = 0; k2 < childrenCount; k2 = k2+1)
				{
					if (Abs((tree1AVL[k])["Label"]))
					{
						(tree1AVL[k])["Label"] = (tree1AVL[k])["Label"] + "|" + (tree1AVL[(((tree1AVL[k])["Children"])[k2])])["Label"];
					}
					else
					{
						(tree1AVL[k])["Label"] = (tree1AVL[(((tree1AVL[k])["Children"])[k2])])["Label"];
					}
				}
			}
			else
			{
				nodeAges [k-1][0] = 0;
				(tree1AVL[k])["Label"] = (tree1AVL[k])["Name"] && 1;
			}
		}




		for (k = 1; k <= size2; k = k+1)
		{
			childrenCount = Abs((tree2AVL[k])["Children"]);
			if (childrenCount)
			{
				(tree2AVL[k])["Label"] = "";
				nodeAges [k-1+size1][0] = -(tree2AVL[(((tree2AVL[k])["Children"])[0])])["Length"];
				(tree2AVL[k])["Length"] = Max(0,(tree2AVL[k])["Length"]) + (tree2AVL[(((tree2AVL[k])["Children"])[0])])["Length"];
				for (k2 = 0; k2 < childrenCount; k2 = k2+1)
				{
					if (Abs((tree2AVL[k])["Label"]))
					{
						(tree2AVL[k])["Label"] = (tree2AVL[k])["Label"] + "|" + (tree2AVL[(((tree2AVL[k])["Children"])[k2])])["Label"];
					}
					else
					{
						(tree2AVL[k])["Label"] = (tree2AVL[(((tree2AVL[k])["Children"])[k2])])["Label"];
					}
				}
			}
			else
			{
				nodeAges [k-1+size1][0] = 0;
				(tree2AVL[k])["Label"] = (tree2AVL[k])["Name"] && 1;
			}
		}


		/*sort the ages of internal nodes, from the youngest (0) to the oldest (roots)*/


		nodeAges 		= nodeAges % 0;
		

		/*--------------------------------------------------------------------------------*/
		/*sort the element of nodeAges according to their depth*/
		AGES={};
		AGES+nodeAges[0][0];
		for (nodeId=0; nodeId<(size1+size2); nodeId=nodeId+1)
		{
			IN=0;
			for(kk=0; kk<Abs(AGES); kk=kk+1)
			{
				if (Abs((nodeAges[nodeId][0])-(AGES[kk]))<1e-16)
				{
					IN=1;
					break;
				}

			} 
			if (IN==0)
			{
				AGES+nodeAges[nodeId][0];
			}
		}
		A={size1+size2,2};
		new_nodeAges=A["_MATRIX_ELEMENT_ROW_"];
		age_count=0;
		for (el=0;el<Abs(AGES);el=el+1)
		{
			nodes_by_age={};
			for (nodeId=0; nodeId<(size1+size2); nodeId=nodeId+1)
			{
				if (Abs((nodeAges[nodeId][0])-(AGES[el]))<1e-16)
				{
					nodes_by_age+nodeAges[nodeId][1];
				}
			}
			SIZE={Abs(nodes_by_age),2};
			nodes_by_age_matrix=SIZE["_MATRIX_ELEMENT_ROW_"];
			for (kk=0;kk<Abs(nodes_by_age);kk=kk+1)
			{
				nodes_by_age_matrix[kk][0]=kk;
				nodes_by_age_matrix[kk][1]=-(nodes_by_age[kk]);		
			}


			nodes_by_age_matrix=nodes_by_age_matrix%1;

			for (kk=0;kk<Abs(nodes_by_age);kk=kk+1)
			{
				nodeAges[age_count][1]=-(nodes_by_age_matrix[kk][1]);
				age_count=age_count+1;			
			}

		}

		/*-----------------------------------------------------------------------------------*/
		
		/*----------------construct the combined tree-----------------------------------------*/
		reportProgress ("Combine the two phylogenetic trees");
		interactionTree = {};


		nodes1			  = {};
		nodes2			  = {};
		idToLabel		  = {};

		insertInteractionPair ("0", -1, (tree1AVL [size1])["Label"] + "," + (tree2AVL [size2])["Label"], size1, size2, nodeAges[0][0]);

		nodesByLevel = {};

		lastNodeAge  = -1e100;


		for (splits = 0; splits<Rows(nodeAges); splits = splits+1)
		{
			nodeAge = nodeAges[splits][0];
			if (nodeAge > lastNodeAge) /*if the nodeage is equal to 0*/
			{

				if (Abs(nodesByLevel [lastNodeAge]) == 0)
				{
					nodesByLevel [nodeAge] = {};
	
				}
				else
				{
					nodesByLevel [nodeAge] = nodesByLevel [lastNodeAge];
	
				}
				lastNodeAge = nodeAge;
			}

			if (nodeAges[splits][1] < size1) /*if the node (here splits) belongs to tree1*/
			{
				idx 		  = nodeAges[splits][1]+1;
				childrenCount = Abs((tree1AVL[idx])["Children"]);
				keys		  = Rows (nodes1[idx]);
				for (dependents = 0; dependents < Abs (nodes1[idx]); dependents = dependents + 1)
				{
	
					interactionNode = 0+keys[dependents];
					nodesByLevel 	 [nodeAge] - interactionNode;
					if (childrenCount)
					{
						for (children = 0; children < childrenCount; children = children + 1)
						{
							childIdx		= ((tree1AVL[idx])["Children"])[children];
							secondIdx		= (interactionTree[interactionNode]) ["T2"];
							nn = insertInteractionPair ("" + Abs(interactionTree), 
													interactionNode, (tree1AVL [childIdx])["Label"] + "," + (tree2AVL [secondIdx])["Label"],
													childIdx, secondIdx, nodeAges[splits][0]
													);
							(nodesByLevel 	 [nodeAge])[nn] = 1;
						}	
					}
					else
					{
						(interactionTree[interactionNode])["Length"] = 0-(interactionTree[interactionNode])["Age"];
					}
				}

			}
			else
			{
				idx = nodeAges[splits][1]-size1+1;
				childrenCount = Abs((tree2AVL[idx])["Children"]);
				keys		  = Rows (nodes2[idx]);

				for (dependents = 0; dependents < Abs (nodes2[idx]); dependents = dependents + 1)
				{
					interactionNode = 0+keys[dependents];
					(nodesByLevel 	 [nodeAge]) - interactionNode;
					if (childrenCount)
					{
						for (children = 0; children < childrenCount; children = children + 1)
						{
							childIdx		= ((tree2AVL[idx])["Children"])[children];
							secondIdx		= (interactionTree[interactionNode]) ["T1"];
							nn = insertInteractionPair ("" + Abs(interactionTree),  
													interactionNode, (tree1AVL [secondIdx])["Label"] + "," + (tree2AVL [childIdx])["Label"],
													secondIdx, childIdx, nodeAges[splits][0]
													);
							(nodesByLevel 	 [nodeAge])[nn] = 1;
						}
					}
					else
					{
						(interactionTree[interactionNode])["Length"] = 0-(interactionTree[interactionNode])["Age"];
					}
				}
			}
		}



		interactionTreeString         = subtreeStringFunction (0,0,0);




		/*------------------------------------------------------------------------------*/

		function insertInteractionPair (name ,parent,label, intree1, intree2, age)
		{
			if (parent >= 0)
			{
				(interactionTree[parent])["Length"] = age - (interactionTree[parent])["Age"];
			}

			idToLabel[name] = label;
			newNode = {"Name":   name,
					   "Length": 0,
					   "Parent": parent,
					   "Label":  label,
					   "T1":     intree1,
					   "T2":     intree2,
					   "Children": 0,
					   "Age":	 age};
					   
					 
			interactionTree + newNode;
			insertionIndex = Abs (interactionTree)-1;
			if (parent >= 0)
			{
				if (Abs(	(interactionTree[parent])["Children"] ) == 0)
				{
						(interactionTree[parent])["Children"] = {};
				}
				((interactionTree[parent])["Children"] )+ insertionIndex;
			}

			if (Abs(nodes1[intree1]) == 0)
			{
				nodes1[intree1] = {}; 
			}
			(nodes1[intree1]) [insertionIndex] = 1;
			nodes1[intree1]  - parent;

			if (Abs(nodes2[intree2]) == 0)
			{
				nodes2[intree2] = {}; 
			}
			(nodes2[intree2]) [insertionIndex] = 1;
			nodes2[intree2]  - parent;

			return insertionIndex;
		}

		/*------------------------------------------------------------------------------*/

		function subtreeStringFunction (node,current,childrenCount)
		{
			childrenCount = Abs((interactionTree[node])["Children"]);
			if (childrenCount == 0)
			{
				sequenceID = (interactionTree[node])["Name"];
				sequenceData = sequenceData + "\n>" + sequenceID + "\n" + interactionsDict[idToLabel[sequenceID]];
				return sequenceID + ":" + (interactionTree[node])["Length"];
			}
			else
			{
				toExec = "\"(\"+subtreeStringFunction(" + ((interactionTree[node])["Children"])[0] + ",current,childrenCount)";
				for (current = 1; current < childrenCount; current = current+1)
				{
					toExec = toExec + "+\",\"+subtreeStringFunction(" + ((interactionTree[node])["Children"])[current] + ",current,childrenCount)";
				}
				ExecuteCommands ("current="+toExec + "+\")" + (interactionTree[node])["Name"] + ":" + (interactionTree[node])["Length"] + "\"");
				return current;
			}
		}

		/*------------------------------------------------------------------------------*/

		/*this is just to see if we got the right number of leaf nodes in the combined tree*/
		count3=0;
		leaves3={};
		for (node=0; node<Abs(interactionTree);node=node+1)
		{
			if (Abs((interactionTree[node])["Children"])==0)
			{
				leaves3[count3]=node;
				count3=count3+1;
			}
		}
		/*fprintf (stdout, "counttree \t", count3,"\n");*/
		

		/*-------------------------------ANALYSIS and model paramter inference----------------------------*/
		reportProgress("Parameter inference");
		/* 1.Read in the data and store the result in a dataset variable.*/
		
		DataSet myData1 = ReadFromString (sequenceData);

		/* 2. Filter the data, specifying that all of the data is to be used
			  and that it is to be treated as nucleotides.*/

		DataSetFilter myFilter1 = CreateFilter (myData1,1);

		/* 3. Collect observed frequencies from the filtered data. observedFreqs1 will
			  store the vector of frequencies. */

		HarvestFrequencies (obsFreqs1, myFilter1, 1, 1, 1);
		
		/* define paramters*/
		global P = obsFreqs1[0];
		P :< 1;
		mlFreqs = {{P,1-P}}; /*The same as obsFreqs1 but in a line */
		global mu=0;
		/* 4. Define the substitution matrix. '*' is defined to be -(sum of off-diag row elements) */
		RateMatrix =
					{{*,mu*t}
					 {mu*t,*}};	
		
		
		/*5.  Define the model, by combining the substitution matrix with the vector of observed (equilibrium)
			  frequencies. */

		Model 	mymodel = (RateMatrix, mlFreqs);

		/*6.  Now we can define the tree variable, using the tree string read from the data file*/

		Tree givenTree1 = interactionTreeString;

		ReplicateConstraint ("this1.?.t:=this2.?.t__", givenTree1, givenTree1); 

		/* construct the likelihood function. */

		LikelihoodFunction  theLnLik1 = (myFilter1, givenTree1);

		/*8.  Maximize the likelihood function, storing parameter values in the matrix paramValues */
            
		Optimize (paramValues1, theLnLik1);
		

		/*9.  Print the value of likelihood in a file*/
		fprintf (likelihood_file_name,Join(",", {{run, it_range,lambda,K,mu,P,paramValues1[1][0]}}), "\n");
		
		/*----------------------------------------------------------------------------------------------*/
				
		/*reconstruct the ancestors*/
		reportProgress("Ancestors reconstruction");
		DataSet mlAncestors= ReconstructAncestors(theLnLik1);
		DataSetFilter ancestralFilter=CreateFilter(mlAncestors,1);
		GetInformation(ancestralsequences,ancestralFilter);
        GetString       (internal_node_names,ancestralFilter,-1);
        _node_count         = Columns(internal_node_names);
        
		GetInformation(observedsequences,myFilter1);
        GetString       (tip_node_names,myFilter1,-1);
        
        _tip_count          = Columns (tip_node_names);
        _simulation_results = {_tip_count, _simulation_replicates+1};
		
		
        data_file_name="`BASE_DIRECTORY`./Results/interactions_network_"+run+"_lambda_"+lambda+"_K_"+K+".csv";
        fprintf (data_file_name, CLEAR_FILE, KEEP_OPEN, "Species1,Species2,Interaction");

        for (_node_index=0; _node_index<_node_count; _node_index += 1) {
            fprintf (data_file_name, "\n", idToLabel[internal_node_names[_node_index]], ",", ancestralsequences[_node_index]);
        }

        for (_node_index=0; _node_index<_tip_count; _node_index += 1) {
            _simulation_results[_node_index][0] = 0+observedsequences[_node_index];
        }

        fprintf (data_file_name, CLOSE_FILE);
		
	
	

		/*-----------------------------------------------------------------------------------------*/
				/*SIMULATIONS*/

		/*simulate under the Lkfunction*/
		reportProgress("Simulations");
        data_file_name="`BASE_DIRECTORY`./Results/simulations_network_"+run+"_lambda_"+lambda+"_K_"+K+".csv";
        fprintf (data_file_name, CLEAR_FILE, KEEP_OPEN);
        
        _columns = {_simulation_replicates + 1, 1};
        _columns[0] = "Observed";

		for (_i=1; _i<=_simulation_replicates; _i+=1) {
		    _columns [_i] = "Replicate" + _i;
            SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (_i) + " of " + _simulation_replicates + " generated", 0);
            DataSet sim    = SimulateDataSet (theLnLik1);
            DataSetFilter   myFilter2 = CreateFilter (sim,1);
            GetInformation  (simulated_states,myFilter2);
            

            for (_node_index=0; _node_index<_tip_count; _node_index += 1) {
                _simulation_results[_node_index][_i] = 0+simulated_states[_node_index];	
            }
        }

        fprintf (data_file_name, "Species1,Species2,", Join (",", _columns), "\n");
		for (_i=0; _i<_node_count; _i+=1) {
            fprintf (data_file_name, idToLabel[tip_node_names[_i]], ",", 
                            Join (",", _simulation_results [_i][-1]), "\n");
        }
        
        fprintf (data_file_name, CLOSE_FILE);

    }

}

function reportProgress (text) {
    fprintf (stdout, "[", text, "]\n");
    return 0;
}

