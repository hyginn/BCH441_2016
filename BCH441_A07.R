# BCH441_A07.R
#
# Purpose:  Supporting scripts for BCH441 (Bioinformatics) at the University of
# Toronto, Fall 2016 - Assignment 07
#
# Version: 1.0
#
# Date:    2016  11
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First code
#
# TODO:
#
#
# ==============================================================================

# DO NOT SIMPLY  source()  THIS FILE!

# ==============================================================================

# Tidbit: Here is the formula to calculate the number of trees one can create
# from n OTUs, as an R function. It's the number of unrooted binary trees
# with n labeled leaves, and unlabeled internal nodes.
#
nTrees <- function(nOTU) {
    if (nOTU < 3)  { return(1) }
    if (nOTU > 87) { return(Inf) }
    return(factorial((2 * nOTU) - 4) / ((2 ^ (nOTU - 2)) * factorial(nOTU - 2)))
}
nTrees(5)  # 15
nTrees(10) # 2,027,025
nTrees(22) # approximately Loschmidt's number

# ==============================================================================
#        PART ONE: Choosing sequences
# ==============================================================================

# Start by loading libraries. You already have the packages installed.
library(Biostrings)
library(msa)
library(stringr)

# What is the latest version of myDB that you have saved?
list.files(pattern = "myDB.*")

# ... load it (probably myDB.05.RData - if not, change the code below).
load("myDB.05.RData")

# The database contains the ten Mbp1 orthologues from the reference species
# and the Mbp1 RBM for YFO.
#
# We will construct a phylogenetic tree from the proteins' APSES domains.
# You have annotated their ranges as a feature.

# Collect APSES domain sequences from your database. The function
# dbGetFeatureSequence() retrieves the sequence that is annotated for a feature
# from its start and end coordinates. Try:

dbGetFeatureSequence(myDB, "MBP1_SACCE", "APSES fold")

# Lets put all APSES sequences into a vector:
APSESnames <- myDB$protein$name[grep("^MBP1_", myDB$protein$name)]
APSES <- character(length(APSESnames))

for (i in 1:length(APSESnames)) {
    APSES[i] <- dbGetFeatureSequence(myDB, APSESnames[i], "APSES fold")
}

# Let's name the rows of our vector with the BiCode part of the protein name.
# This is important so we can keep track of which sequence is which. We use the
# gsub() funcion to substitute "" for "MBP1_", thereby deleting this prefix.
names(APSES) <- gsub("^MBP1_", "", APSESnames)

# inspect the result: what do you expect? Is this what you expect?
head(APSES)

# Let's add the E.coli Kila-N domain sequence as an outgroup, for rooting our
# phylogegetic tree (see the Assignment Course Wiki page for details on the
# sequence).

APSES[length(APSES) + 1] <-
"IDGEIIHLRAKDGYINATSMCRTAGKLLSDYTRLKTTQEFFDELSRDMGIPISELIQSFKGGRPENQGTWVHPDIAINLAQ"
names(APSES)[length(APSES)] <- "ESCCO"


# ==============================================================================
#        PART TWO: Multiple sequence alignment
# ==============================================================================

# This vector of sequences with named elements fulfills the requirements to be
# imported as a Biostrings object - an AAStringSet - which we need as input for
# the MSA algorithms in Biostrings.
#

APSESSeqSet <- AAStringSet(APSES)

APSESMsaSet <- msaMuscle(APSESSeqSet, order = "aligned")

# inspect the alignment.
writeSeqSet(APSESMsaSet, format = "ali")


# What do you think? Is this a good alignment for phylogenetic inference?

# ==============================================================================
#        PART THREE: reviewing and editing alignments
# ==============================================================================

# Head back to the assignment 7 course wiki page and read up on the background
# first.
#



# Let's mask out all columns that have observations for
# less than 1/3 of the sequences in the dataset. This
# means they have more than round(nrow(msaSet) * (2/3))
# hyphens in a column.
#
# We take all sequences, split them into single
# characters, and put them into a matrix. Then we
# go through the matrix, column by column and decide
# whether we want to include that column.

# Step 1. Go through this by hand...

# get the length of the alignment
lenAli <- APSESMsaSet@unmasked@ranges@width[1]

# initialize a matrix that can hold all characters
# individually
msaMatrix <- matrix(character(nrow(APSESMsaSet) * lenAli),
                    ncol = lenAli)

# assign the correct rownames
rownames(msaMatrix) <- APSESMsaSet@unmasked@ranges@NAMES
for (i in 1:nrow(APSESMsaSet)) {
    seq <- as.character(APSESMsaSet@unmasked[i])
    msaMatrix[i, ] <- unlist(strsplit(seq, ""))
}

# inspect the result
msaMatrix[1:5, ]

# Now let's make a logical vector with an element
# for each column that selects which columns should
# be masked out.

# To count the number of elements in a vector, R has
# the table() function. For example ...
table(msaMatrix[ , 1])
table(msaMatrix[ , 10])
table(msaMatrix[ , 20])
table(msaMatrix[ , 30])


# Since the return value of table() is a named vector, where
# the name is the element that was counted in each slot,
# we can simply get the counts for hyphens from the
# return value of table(). We don't even need to assign
# the result to an intermediate variable, but we
# can attach the selection via square brackets,
# i.e.: ["-"],  directly to the function call:
table(msaMatrix[ , 1])["-"]

# ... to get the number of hyphens. And we can compare
# whether it is eg. > 4.
table(msaMatrix[ , 1])["-"] > 4

# Thus filling our logical vector is really simple:

# initialize the mask
colMask <- logical(lenAli)

# define the threshold for rejecting a column
limit <- round(nrow(APSESMsaSet) * (2/3))

# iterate over all columns, and write TRUE if there are less-or-equal to "limit"
# hyphens, FALSE if there are more.
for (i in 1:lenAli) {
    count <- table(msaMatrix[ , i])["-"]
    if (is.na(count)) { # No hyphen
        count <- 0
    }
    colMask[i] <- count <= limit
}

# inspect the mask
colMask

# How many positions were masked? R has a simple trick
# to count the number of TRUE and FALSE in a logical
# vector. If a logical TRUE or FALSE is converted into
# a number, it becomes 1 or 0 respectively. If we use
# the sum() function on the vector, the conversion is
# done implicitly. Thus ...
sum(colMask)

# ... gives the number of TRUE elements.

cat(sprintf("We are masking %4.2f %% of alignment columns.\n",
            100 * (1 - (sum(colMask) / length(colMask)))))


# Next, we use colMask to remove the masked columns from the matrix
# in one step:
maskedMatrix <- msaMatrix[ , colMask]

# check:
ncol(maskedMatrix)


# ... then collapse each row back into a sequence ...

apsMaskedSeq <- character()
for (i in 1:nrow(maskedMatrix)) {
    apsMaskedSeq[i] <- paste(maskedMatrix[i, ], collapse="")
}
names(apsMaskedSeq) <- rownames(maskedMatrix)

# ... and read it back into an AAStringSet object

apsMaskedSet <- AAStringSet(apsMaskedSeq)

# inspect ...
writeSeqSet(apsMaskedSet, format = "ali")



# Step 2. Turn this code into a function...

# Even though the procedure is simple, doing this more than once is tedious and
# prone to errors. I have assembled the steps we just went through into a
# function maskSet() and put it into the utilities.R file, from where it has
# been loaded when you started this sesssion.

maskSet

# Check that the function gives identical results
# to what we did before by hand:
identical(apsMaskedSet, maskSet(APSESMsaSet))

# The result must be TRUE. If it's not TRUE you have
# an error somewhere.

# We save the aligned, masked domains to a file in multi-FASTA format.
writeSeqSet(maskSet(APSESMsaSet), file = "APSES.mfa",   format = "mfa")

# ==============================================================================
#        PART FOUR: Calculating trees
# ==============================================================================

# After you have installed Phylip on your computer, install the R package that
# provides an interface to the Phylip functions.

if (!require(Rphylip, quietly=TRUE)) {
    install.packages("Rphylip")
    library(Rphylip)
}

# This will install RPhylip, as well as its dependency, the package "ape".

# The next part may be tricky. You will need to figure out where
# on your computer Phylip has been installed and define the path
# to the proml program that calculates a maximum-likelihood tree.
# I give you instructions for the Mac below.
# You'll need to figure out the equivalent Windows commands and
# please post instructions on the mailing list once you have got
# this to work on Windows.

# On the Mac, the standard installation places a phylip folder
# in the /Applications directory. That folder contains all the
# individual phylip programs as <name>.app files. These are not
# the actual executables, but "app" files are actually directories
# that contain the required resources for a program to run.

# The executable is in a subdirectory and you can point Rphylip
# directly to that subdirectory to find the program it needs:
PROMLPATH <- "/Applications/phylip-3.695/exe/proml.app/Contents/MacOS"

# Now read the mfa file you have saved, as a "proseq" object with the
# read.protein() function of the RPhylip package:

apsIn <- read.protein("APSES.mfa")

# ... and you are ready to build a tree.

# Building maximum-likelihood trees can eat as much computer time
# as you can throw at it. Calculating a tree of 48 APSES domains
# with default parameters of Rproml() runs for more than half a day
# on my computer. But we have only twelve sequences here, so the
# process will take us about 5 to 10 minutes.

apsTree <- Rproml(apsIn, path=PROMLPATH)


# A quick first look:

plot(apsTree)

# Back to the Assignment 7 course Wiki page...


# ==============================================================================
#        PART FIVE: Tree analysis
# ==============================================================================

# A Entrez restriction command
cat(paste(paste(c(myDB$taxonomy$ID, "83333"), "[taxid]", sep=""), collapse=" OR "))

# The Common Tree from NCBI
# Download the EDITED phyliptree.phy
commonTree <- read.tree("phyliptree.phy")

# Plot the tree
plot(commonTree, cex=1.0, root.edge=TRUE, no.margin=TRUE)
nodelabels(text=commonTree$node.label, cex=0.6, adj=0.2, bg="#D4F2DA")


# === Visualizing your tree ====================================================

# The trees that are produced by Rphylip are stored as an object of class
# "phylo". This is a class for phylogenetic trees that is widely used in the
# community, practically all R phylogenetics packages will options to read and
# manipulate such trees. Outside of R, a popular interchange format is the
# Newick_format that you have seen above. It's easy to output your calculated
# trees in Newick format and visualize them elsewhere.

# The "phylo" class object is one of R's "S3" objects and methods to plot and
# print it have been added to the system. You can simply call plot(<your-tree>)
# and R knows what to do with <your-tree> and how to plot it. The underlying
# function is plot.phylo(), and documentation for its many options can by found
# by typing:

?plot.phylo

plot(apsTree) # default type is "phylogram"
plot(apsTree, type="unrooted")
plot(apsTree, type="fan", no.margin = TRUE)

# rescale to show all of the labels:
# record the current plot parameters ...
tmp <- plot(apsTree, type="fan", no.margin = TRUE, plot=FALSE)
# ... and adjust the plot limits for a new plot
plot(apsTree,
     type="fan",
     x.lim = tmp$x.lim * 1.8,
     y.lim = tmp$y.lim * 1.8,
     cex = 0.8,
     no.margin = TRUE)

# Inspect the tree object
str(apsTree)
apsTree$tip.label
apsTree$edge
apsTree$edge.length

# show the node / edge and tip labels on a plot
plot(apsTree)
nodelabels()
edgelabels()
tiplabels()

# show the number of nodes, edges and tips
Nnode(apsTree)
Nedge(apsTree)
Ntip(apsTree)


# Finally, write the tree to console in Newick format
write.tree(apsTree)

# === Rooting Trees ============================================================

# In order to analyse the tree, it is helpful to root it first and reorder its
# clades. Contrary to documentation, Rproml() returns an unrooted tree.

is.rooted(apsTree)

# You can root the tree with the command root() from the "ape" package. ape is
# automatically installed and loaded with Rphylip.

plot(apsTree)

# add labels for internal nodes and tips
nodelabels(cex=0.5, frame="circle")
tiplabels(cex=0.5, frame="rect")

# The outgroup of the tree is tip "8" in my sample tree, it may be a different
# number in yours. Substitute the correct node number below for "outgroup".
apsTree <- root(apsTree, outgroup = 8, resolve.root = TRUE)
plot(apsTree)
is.rooted(apsTree)

# This tree _looks_ unchanged, beacuse when the root trifurcation was resolved,
# an edge of length zero was added to connect the MRCA (Most Recent Common
# Ancestor) of the ingroup.

# The edge lengths are stored in the phylo object:
apsTree$edge.length

# ... and you can assign a small arbitrary value to the edge
# to show how it connects to the tree without having an
# overlap.
apsTree$edge.length[1] <- 0.1
plot(apsTree, cex=0.7)
nodelabels(text="MRCA", node=12, cex=0.5, adj=0.1, bg="#ff8866")


# This procedure does however not assign an actual length to a root edge, and
# therefore no root edge is visible on the plot. Why? , you might ask. I ask
# myself that too. We'll just add a length by hand.

apsTree$root.edge <- mean(apsTree$edge.length) * 1.5
plot(apsTree, cex=0.7, root.edge=TRUE)
nodelabels(text="MRCA", node=12, cex=0.5, adj=0.8, bg="#ff8866")


# === Rotating Clades ==========================================================

# To interpret the tree, it is useful to rotate the clades so that they appear
# in the order expected from the cladogram of species.

# We can either rotate around individual internal nodes:
layout(matrix(1:2, 1, 2))
plot(apsTree, no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#ff8866")
plot(rotate(apsTree, node=17), no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#88ff66")
layout(matrix(1), widths=1.0, heights=1.0)

# ... or we can plot the tree so it corresponds as well as possible to a
# predefined tip ordering. Here we use the ordering that NCBI Global Tree
# returns for the reference species - we have used it above to make the vector
# apsMbp1Names. You inserted your YFO name into that vector - but you should
# move it to its correct position in the cladogram.

# (Nb. we need to reverse the ordering for the plot. This is why we use the
# expression [nOrg:1] below instead of using the vector directly.)

nOrg <- length(apsTree$tip.label)

layout(matrix(1:2, 1, 2))
plot(commonTree,
     no.margin=TRUE, root.edge=TRUE)
nodelabels(text=commonTree$node.label, cex=0.5, adj=0.2, bg="#D4F2DA")

plot(rotateConstr(apsTree, apsTree$tip.label[nOrg:1]),
     no.margin=TRUE, root.edge=TRUE)
add.scale.bar(length=0.5)
layout(matrix(1), widths=1.0, heights=1.0)

# Study the two trees and consider their similarities and differences. What do
# you expect? What do you find?
#

# Print the two trees on one sheet of paper, write your name and student number,
# and bring it to class as your deliverable for this assignment.




# [END]
