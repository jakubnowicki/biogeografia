###################################################################################################
# Diversity Analysis using Vegan and simple R functions
# written by M. Kowalewski, University of Florida
# kowalewski@ufl.edu
# November 17, 2014
###################################################################################################

library(vegan)			# you need to upload vegan for parts of this script
library(lattice)			# you need to upload lattice for parts of this script

# For details see: http://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf

###################################################################################################
# let's use the dataset we have used previously in ordination exercies ("dummyorddata.csv")
a<- read.csv("dummyorddata.csv")		# a dataset (species=columns and samples=rows). Values are counts.
							# Note that the first variable is sample label, so we need to
							# remove it a[,-1] to ensure that vegan functions work properly
a
###################################################################################################
# DIVERSITY INDICES
# we will use here functions "diversity", "fisher.alpha", and "specnumber".
# All these indices are pretty simple to compute so you can also develop your own code. This is also
# just a sample of indices (many other indices exist, although they tend to be quite redundant)
# Notice that some of those indices cannot be computed if species richness is 1
# I also added a one-line statement to compute Parker-Berger Index [PB] (the relative abundance of the most abundant taxon)

n<- apply(a[,-1],1,sum)				# compute sample sizes
PB<- apply(a[,-1]/n,1,max)			# compute Parker-Berger index
richness<- specnumber(a[,-1])			# per sample species richness (=number of species in a sample)
S_a<- diversity(a[,-1], index="simpson")	# simpson diversity index
H_a<- diversity(a[,-1], index="shannon")	# shannon diversity index (H)
fisher<- fisher.alpha(a[,-1],se=FALSE)	# Fisher's alpha with standard errors
J<- H_a/log(richness)				# Evenness index J (based on H)
final<- cbind(a[,1],n,richness,S_a,H_a,J,fisher,PB)	# assemble all results into a single matrix
colnames(final)<- c("sample","no.speciemns","no.species","Simpson","Shannon",
			  "Shannon-J", "Fisher's alpha","Parker-Berger")	# label columns
final							# print final output
cor(final[-which(final[,3]<2),-c(1,2,3,7)])	# check correlation of indices (remove samples with richness<2)

###################################################################################################
# RAREFACTION

# For multiple samples, we can reduce all sample sizes to the same number of individuals:
a_rar <- rarefy(a[,-1], min(rowSums(a[,-1])))
a_rar
# We can use rarecurve function to plot rarefaction curves of all compared samples
rarecurve(a[,-1], step=1, xlab = "Sample Size", ylab = "Species")

###################################################################################################
# RAD MODELS
# This set of functions in vegan provides maximum likelihood estimates for some common RAD models

rad3<- radfit(a[,3])	# we can run it for a given sample
rad3
radall<- radfit(apply(a[,-1],2,sum))		# we can also run it for pooled data
radall
plot(radall)			# you can plot all models against the actual data
radlattice(radall)		# or you can plot it as a series of plots
str(radall)				# you can examine the content of radall and extract values you need. For example:
plot(1:ncol(a[,-1]),log2(radall$y),xlab="species rank order",ylab="species abundance", pch=16, col="black")
points(1:ncol(a[,-1]),log2(radall$models$Null$fitted.values),type="l", col="red")
plot(1:ncol(a[,-1]),radall$y,xlab="species rank order",ylab="species abundance", pch=16, col="black",log="y")
points(1:ncol(a[,-1]),radall$models$Null$fitted.values,type="l", col="red")

