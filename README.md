NOTES FOR INFERELATOR
====================
Overall comments:
- Inferelator.R is v1, which uses lars and is more or less finished
- Inf-v2.R is v2, which will use glmnet and is under active development

Requirements:
- A recent version (2.8?) of R,
- glmnet 
- lars (has some internal fxns used by cv.glmnet())
- inf-v2.R, which has cv.glmnet()

Features for v2:
- Get rid of correlated predictors, removes need for preclust() in v1
- Changeable alphas ("alpha" in glmnet(), 0 = ridge-regression, 1 = lasso)
- Changeable lambdas 
- Variable penalty ("penalty.factor" in glmnet())
- Clean up output, v1 is a little messy w/ separate .noa, .eda and .sif files

Usage:
There are 3 major parts to using Inferelator version 1:
I. learning a network from expression data
II. validating the network by comparing predicted profiles with observations
III. generating output files & plotting

I. Learning a network from expression data
------------------------------------------
1. Start up R
2. source("inferelator.R")
3. inf.result = preclust.and.inferelate(clusterStack=clusterStack.egrin, data=ratios.egrin, predictors=halo.tfs, col.map=col.map.egrin)
This will run with the default parameters. You may want to change some of them, such as 'preclust' and 'filter.by.aic'
To run with environmental factors set 'data=rbind(ratios.redox,env.data' and 'predictors=c(halo.tfs,env.factors)'
When the run finishes there should be tfGroups.noa and types.noa in the directory where R was launched.

II. Validating the network
--------------------------
1. Load the new ratios and its corresponding col.map (unnecessary if validating on training set)
2. pred = predictelator(coeffs=inf.result, clusterStack=clusterStack.egrin, oldratios=ratios.egrin, newratios=ratios.egrin, tfgroups=inf.result$tf.groups, col.map=col.map.redox) ##tf.file='tfgroups.noa'
The output pred is a m x n matrix of m mean bicluster profiles across n conditions
There are also some global vars:
tmp3.steady & tmp3.timeseries, one of which should be identical to your output depending on whether col.map was supplied
mean.old.profiles : mean bicluster profiles

III. Generating output files and plotting
-----------------------------------------
1. make.network.files(inf.result=inf.result, clusterStack=clusterStack.redox, sif.filename="myfile.sif")
This will make the remaining .noa and .eda files. IMPORTANT: saves to current directory, will overwrite existing files
2. count.predictors(inf_result)
This will return a vector of all unique single and combinatory predictors found in the result
3. plot.prediction(k)
This will plot observed bicluster profile, predicted profiles w/ and w/o timeseries correction in the same plot 