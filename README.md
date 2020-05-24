Code to infer model parameters in the county-level metapopulation model for COVID-19 spread in the US.

The transmission model is run in cpp. Need to mex the cpp function before use using "mex model_eakf.cpp" in MATLAB.

Inter-county commuting is stored using neighbor list: nl, part and C.
nl (neighbor list) and part (partition) both describe the network structure. 
For instance, the neighbors of location i are nl(part(i):part(i+1)-1). This neighbor set include i itself.
C has the same size with nl. The commuters from location i to location j=nl(part(i)+x) is C(part(i)+x). 
