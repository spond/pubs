This code implements the probabilistic model of evolution for ecological networks
(reference to be inserted)

1. Download and install HyPhy (www.hyphy.org). You will need version 2.13 or later, which 
may only be available as a source download. You can always obtain the latest versions of 
HyPhy from https://github.com/veg/hyphy (the qtbranch)

2. Run Ecological_Network_Evolution.bf through HyPhy, e.g. ```$HYPHYMP ./Ecological_Network_Evolution.bf```

3. Inputs should also be in a folder called ```Inputs``` (case sensitive on some systems) 
in the same directory as ```Ecological_Network_Evolution.bf```. The script will over all 
the networks which are defined in the ```Inputs``` folder (batch processing). 
For each integer N = 1..K (K = number of networks to process), there should be:
 - two files contaning respectively the 2 phylogenetic trees: for the set of 'plants' and the set of 'animals'
 (or host/parasite, or predator/prey etc).
   They are given in their Newick representation (branch lengths need not be specified). 
    Files are respectively named: AN_N.txt and PL_N.txt , where N is a number denoting the specific network.  
 - an interaction file which gives the interactions at present time. It should be named: interaction_N.txt
   Each line of the file contains repeatedly: the animal name (line 1n), the plant name (line 2n) and the state of the interaction -0 or 1- (line 3n), 
   with n possible pairs of animal and plant. Animal names and plant names should be the same as those at the tips of the phylogenetic trees (NOT
   case sensitive).


4. It is  important to create a folder called ```Results``` 
(case sensitive on some systems) in the same directory as ```Ecological_Network_Evolution.bf```, 
before running the analysis. Output files (in the ```Results``` folder) are:
 - A single comma-separated file (```models.txt```), which, for each network/parameter setting set, 
 all model parameter estimates and the log likelihood. (**) below encode the network/parameter set,
 e.g. interactions_network_1_lambda_1_K_0.25
 - A interactions(**).csv file, listing, the inferred value for each pair of _unobserved_ interactions
    An entry of the form 
    ```MEGACHILE_SEMI,CYNOGLOSSUM_|PHACELIA_SECUNDA,1```
    implies that the model inferred an existing interaction between 'MEGACHILE_SEMI', and the most
    recent common ancestor of  CYNOGLOSSUM_ and PHACELIA_SECUNDA.
 - A simulations(**).csv file, listing, for each _observed_ interaction, the observed value, and 
    values based on '_simulation_replicates' (100 by default, edit ```Ecological_Network_Evolution.bf``` to change) 
    simulations from the likelihood function inferred for a given
    set of model parameters/global (lambda,K) values.
 
