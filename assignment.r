# Assignment 9 #
# Martin Gubler #

# 1. Consider the state sequence object biofam.seq and the matrix dOM of pairwise OM
# dissimilarities based on the properties matrix considered in the previous assignments.

### CODE COPY-PASTE FROM PREVIOS ASSIGNMENT / SAMPLE SOLUTION

library(TraMineR)
data(biofam)
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
                     labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"),
                     right=FALSE)
bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child",
               "Left/Child", "Left/Married/Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
biofam.seq <- seqdef(biofam[,10:25], states=bf.shortlab,
                     labels=bf.states, weights=biofam$wp00tbgs)

properties <- matrix(c(# left, married, child, divorced
  0, 0, 0, 0,  # parent
  1, 0, 0, 0,  # left
  0, 1, .5, 0, # marr
  1, 1, 0, 0,  # left+marr
  0, 0, 1, 0,  # child
  1, 0, 1, 0,  # left+child
  1, 1, 1, 0,  # left+marr+child
  .5, 1, .5, 1 # divorced
), 8, 4, byrow=TRUE)
sm <- as.matrix(dist(properties))
indel <- .5*max(sm)
dOM <- seqdist(biofam.seq, method="OM", indel=indel, sm=sm,
               full.matrix = FALSE)
weight <- attr(biofam.seq, "weight")

# 2. Plot the single representative of the whole dataset for each of the `frequency' and
# `centrality' criteria.

par(mfrow = c(2,2))
seqrplot(biofam.seq, dist.matrix = dOM, criterion = "freq", nrep = 1)
seqrplot(biofam.seq, dist.matrix = dOM, criterion = "dist", nrep = 1)

# 3. Plot the densest-neighborhood representative of the whole dataset by setting 
# successively the neighborhood radius tsim as :1; :2; :3; :4.

seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", nrep = 1, tsim = 0.1)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", nrep = 1, tsim = 0.2)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", nrep = 1, tsim = 0.3)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", nrep = 1, tsim = 0.4)

# 4. Plot the set of representatives using a neighborhood radius of :15 and setting
# successively the minimum wanted coverage (trep) as :25; :4; :5; :75; 1

seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", tsim = 0.15, trep = 0.25)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", tsim = 0.15, trep = 0.4)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", tsim = 0.15, trep = 0.5)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", tsim = 0.15, trep = 0.75)
seqrplot(biofam.seq, dist.matrix=dOM, criterion = "density", tsim = 0.15, trep = 1)

# 5. Consider the 6-group solution of the PAM clustering of the set of biofam sequences
# and represent the clusters with representative plots using a neighborhood radius of
# .15 and a minimum coverage of 70%.

  # COPY/PASTE TO GET PAM SOLUTION

library(WeightedCluster)
set.seed(13425)
cluster.pam6 <- wcKMedoids(dOM, k = 6, weight = weight)
cl.pam6 <- factor(cluster.pam6$clustering)
pam.labels <- c("Staying with parents","Late Parenthood","Married without child","Solo","Staying married with parents","Early Parenthood")
cl.pam6.fact <- factor(cl.pam6, labels=pam.labels)

# Start plot

seqrplot(biofam.seq, dist.matrix = dOM, group = cl.pam6.fact, trep = 0.7, tsim = 0.15, border = NA)

#### QUESTION:  I do not provide a criterion but I seem to get the same solution as in the sample, where "density" is provided.
####            Does that mean that "density" is the default criterion used?


# 6. Examine the quality measures of the representatives. How do you explain that the
# overall Q is negative when there is a single representative?

#### I want to make sure that only the latest official version of TraMineR is installed on my computer (1.8-3).
#### I was not sure about the consequences when downloading a development version.
#### Therefore I will try this assignment with the old version, even if the values are false.

install.packages("TraMineRextras", repos="http://R-Forge.R-project.org")
library(TraMineRextras)
q.gr <- seqrep.grp(biofam.seq, group = cl.pam6.fact, mdis = dOM, trep = 0.7, tsim = 0.15)
q.gr[1:4]

#### QUESTION: This command only produces error messages 
#### (This error message re TraMineRextras 'combinat' occurs almost always when I try to use the package - not sure why):
# Loading required package: combinat
# Error: package ‘combinat’ could not be loaded
# In addition: Warning message:
#   In library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc) :
#   there is no package called ‘combinat’
# > q.gr <- seqrep.grp(biofam.seq, group = cl.pam6.fact, mdis = dOM, trep = 0.7, tsim = 0.15)
# Error: could not find function "seqrep.grp"
# > q.gr[1:4]
# Error: object 'q.gr' not found

# 7. Compute the discrepancy (pseudo-variance) of the whole set of data and the within
# dicrepancy in each cluster. (Tip: Currently, dissvar.grp accepts only a matrix as
# first argument. If you have computed the dissimilarity matrix with full.matrix=FALSE 
# pass the dissimilarity matrix as as.matrix(dOM))

#### Did not understand this - look forward to hearing the explanations in the webinar.