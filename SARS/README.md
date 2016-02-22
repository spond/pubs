To generate simulated data as described in "A case for the ancient origin of coronaviruses" 
follow these steps

1. Download and install HyPhy (www.hyphy.org). You will need version 2.13 or later, which 
may only be available as a source download. You can always obtain the latest versions of 
HyPhy from https://github.com/veg/hyphy

2. Run BS-REL_simulator.bf through HyPhy, e.g. ```$HYPHYMP ./BS-REL_simulator.bf```

3. When propmted, 
   - supply the name of the file containing model fit information of the RdRp gene (rdrp.fit)
   - indicate the total desired tree length (expected substitutions per nucleotide/site)
   - choose the number of replicates to be generated 
   - provide the file system path to save the replicates to
