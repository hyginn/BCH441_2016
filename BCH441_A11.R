# BCH441_A11.R
#
# Purpose:  Supporting scripts for BCH441 (Bioinformatics) at the University of
# Toronto, Fall 2016 - Assignment 11
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
#        PART ONE: ...
# ==============================================================================


# This tutorial covers basic concepts of graph theory and analysis in R. You
# should have typed init() to configure some utilities in the background.


# ==============================================================================
#        PART ONE: REVIEW
# ==============================================================================

# I assume you'll have read the Pavlopoulos review of graph theory concepts.
# Let's explore some of the ideas by starting with a small random graph."

set.seed(112358)

N <- 20
# create a vector of N random node names
#
makeRandomGenenames <- function(N) {
    nam <- character()
    while (length(nam) < N) {
        a <- paste(c(sample(LETTERS, 1), sample(letters, 2)),
                   sep="", collapse="") # three letters
        n <- sample(1:9, 1)             # one number
        nam[length(nam) + 1] <- paste(a, n, sep="") # store in vector
        nam <- unique(nam)   # delete if this was a duplicate
    }
    return(nam)
}
Nnames <- makeRandomGenenames(N)
Nnames

# One way to represent graphs in a computer is as an "adjacency matrix". In this
# matrix, each row and each column represents a node, and the cell at the
# intersection of a row and column contains a value/TRUE if there is an edge,
# 0/FALSE otherwise. It's easy to see that an undirected graph has a symmetric
# adjacency matrix (i, j) == (j, i); and we can put values other than {1, 0}
# into a cell if we want to represent a weighted edge.

# At first, lets create a random graph: let's say a pair of nodes has
# probability p <- 0.12 to have an edge, and our graph is symmetric and has no
# self-edges. We use our nNames as node labels, but I've written the function so
# that we could also just ask for a number of un-named nodes, for later.

makeRandomGraph <- function(nam, p = 0.1) {
    # nam: either a character vector of unique names, or a single
    #        number that will be converted into a vector of integers.
    # p:   probability that a random pair of nodes will have an edge.
    #
    # Value: an adjacency matrix
    #
    if (is.numeric(nam) && length(nam) == 1) {
        nam <- as.character(1:nam)
    }
    N <- length(nam)
    G <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
    rownames(G) <- nam
    colnames(G) <- nam
    for (iRow in 1:(N-1)) { # Note how we make sure iRow != iCol
        for (iCol in (iRow+1):N) {
            if (runif(1) < p) {  # runif() creates uniform random numbers
                # between 0 and 1
                G[iRow, iCol] <- 1   # symmetric
                G[iCol, iRow] <- 1
            }
        }
    }
    return(G)
}

set.seed(112358)
G <- makeRandomGraph(Nnames, p = 0.09)
G


# That's not very informative - we should plot this graph. We'll go into more details of the iGraph package a bit later, for now we just use it to plot:

if (!require(igraph)) {
    install.packages("igraph")
    library(igraph)
}

iG <- graph_from_adjacency_matrix(G)
iGxy <- layout_with_graphopt(iG, charge=0.001)   # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 800 + (150 * degree(iG)),
     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)


# The simplest descriptor of a graph are the number of nodes, edges, and the
# degree-distribution. In our example, the number of nodes was given: N; the
# number of edges can easily be calculated from the adjacency matrix. In our
# matrix, we have entered 1 for every edge. Thus we simply sum over the matrix:
sum(G)

# Is that correct? Is that what you see in the plot?

# Yes and no: we entered every edge twice: once for a node [i,j], and again for
# the combination [j, i]. Whether that is correct depends on what exactly we
# want to do with the matrix. If these were directed edges, we would need to
# keep track of them separately. Since we didn't intend them to be directed,
# we'll just divide the number of edges by 2. Why didn't we simply use an
# upper-triangular matrix? Because then we need to keep track of the ordering of
# edges if we want to know whether a particular edge exists or not. For example
# we could sort the nodes alphabetically, and make sure we always query a pair
# in alphabetical order. Then a triangular matrix would be efficient.

# What about the degree distribution? We can get that simply by summing over the
# rows (or the columns):"

rowSums(G)

rs <- rowSums(G)
brk <- seq(min(rs)-0.5, max(rs)+0.5, by=1)
hist(rs, breaks=brk, col="#A5CCD5")

# The degree distribution is actually quite an important descriptor of graphs,
# since it is very sensitive to the generating mechanism. For biological
# networks, that is one of the key questions we are interested in: how was the
# network formed?

# ==============================================================================
#        PART TWO: DEGREE DISTRIBUTIONS
# ==============================================================================

# Let's simulate a few graphs that are a bit bigger to get a better sense of their degree distributions:
#

set.seed(31415927)
G200 <- makeRandomGraph(200, p = 0.015)
iG200 <- graph_from_adjacency_matrix(G200)
iGxy <- layout_with_graphopt(iG200, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG200,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG200)+1))[degree(iG200)+1],
     vertex.size = 200 + (30 * degree(iG200)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

dg <- degree(iG200)/2   # we use the iGraph function, not rowsums()
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5CCD5")

# Note the characteristic peak of this distribution: this is not "scale-free". Here is a log-log plot of frequency vs. degree-rank:

freqRank <- table(dg)
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b")

# the iGraph package has a function to make random graphs according to the Barabasi-Albert model of scale-free graphs.


set.seed(31415927)
GBA <- sample_pa(200, power = 0.8)

iGxy <- layout_with_graphopt(GBA, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(GBA,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(GBA)+1))[degree(GBA)+1],
     vertex.size = 200 + (30 * degree(GBA)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# This is a very obviously different graph!
#
dg <- degree(GBA)/2
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5D5CC")

# Most nodes have a degree of 1, but one node has a degree of 14.

freqRank <- table(dg)
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b")

# sort- of linear, many of the higher ranked nodes have a frequency of only one, that behaviour smooths out in larger graphs.
#
X <- sample_pa(100000, power = 0.8)  # 100,000 nodes
freqRank <- table(degree(X)/2)
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b")
rm(X)

# Finally, let's simulate a random geometric graph and look at the degree distribution. We'll randomly place our nodes in a box. Then we'll define the
# probability for two nodes to have an edge to be a function of their distance.

makeRandomGeometricGraph <- function(nam, B = 25, Q = 0.001, t = 0.6) {
    # nam: either a character vector of unique names, or a single
    #        number that will be converted into a vector of integers.
    # B, Q, t:   probability that a random pair (i, j) of nodes gets an
    #              edge determined by a generalized logistic function
    #              p <- 1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / 0.9)))
    #
    # Value: a list with the following components:
    #        G$mat : an adjacency matrix
    #        G$nam : labels for the nodes
    #        G$x   : x-coordinates for the nodes
    #        G$y   : y-coordinates for the nodes
    #
    nu <- 1  # probably not useful to change
    G <- list()

    if (is.numeric(nam) && length(nam) == 1) {
        nam <- as.character(1:nam)
    }
    G$nam <- nam
    N <- length(G$nam)
    G$mat <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
    rownames(G$mat) <- G$nam
    colnames(G$mat) <- G$nam
    G$x <- runif(N)
    G$y <- runif(N)
    for (iRow in 1:(N-1)) { # Same principles as in makeRandomGraph()
        for (iCol in (iRow+1):N) {
            # geometric distance ...
            d <- sqrt((G$x[iRow] - G$x[iCol])^2 +
                      (G$y[iRow] - G$y[iCol])^2)  # Pythagoras
            # distance dependent probability
            p <- 1 - 1/((1 + (Q * (exp(-B * (d-t)))))^(1 / nu))
            if (runif(1) < p) {
                G$mat[iRow, iCol] <- 1
                G$mat[iCol, iRow] <- 1
            }
        }
    }
    return(G)
}

# getting the parameters of a generalized logistic right takes a bit of
# experimenting. If you are interested, you can try a few variations. Or you can
# look up the function at
# https://en.wikipedia.org/wiki/Generalised_logistic_function

# genLog <- function(x, B = 25, Q = 0.001, t = 0.5) {
#     # generalized logistic (sigmoid)
#     nu <- 1
#     return(1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / nu)))
# }
#
# x <- seq(0, sqrt(2), length.out = 50)
# plot(x, genLog(x), type="l", col="#AA0000", ylim = c(0, 1),
#      xlab = "d", ylab = "p(edge)")


set.seed(112358)
GRG <- makeRandomGeometricGraph(200, t=0.4)


iGRG <- graph_from_adjacency_matrix(GRG$mat)
iGRGxy <- cbind(GRG$x, GRG$y) # use our node coordinates for layout

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])) * 1.1,
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iGRG)+1))[degree(iGRG)+1],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

dg <- degree(iGRG)/2
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#CCA5D5")

freqRank <- table(dg)
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b")

# ====================================================================
#        PART THREE: A CLOSER LOOK AT THE igraph PACKAGE
# ====================================================================


# == BASICS ==========================================================

# The basic object of the igraph package is a graph object. Let's explore the
# first graph some more, the one we built with our random gene names:
summary(iG)

# This output means: this is an IGRAPH graph, with D = directed edges and N =
# named nodes, that has 20 nodes and 40 edges. For details, see
?print.igraph

mode(iG)
class(iG)

# This means an igraph graph object is a special list object; it is opaque in
# the sense that a user is never expected to modify its components directly, but
# through a variety of helper functions which the package provides. There are
# many ways to construct graphs - from adjacency matrices, as we have just done,
# from edge lists, or by producing random graphs according to a variety of
# recipes, called _games_ in this package.

# As with many R objects, loading the package provides special functions that
# can be accessed via the same name as the basic R functions, for example:

print(iG)
plot(iG)

# ... where plot() allows the usual flexibility of fine-tuning the plot. We
# first layout the node coordinates with the Fruchtermann-Reingold algorithm - a
# force-directed layout that applies an ettractive potential along edges (which
# pulls nodes together) and a repulsive potential to nodes (so they don't
# overlap). Note the use of the degree() function to color and scale nodes and
# labels by degree and the use of the V() function to retrieve the vertex names.
# See ?plot.igraph for details."

iGxy <- layout_with_fr(iG)   # calculate layout coordinates

# Plot with some customizing parameters
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 9 + (2 * degree(iG)),
     vertex.label.cex = 0.5 + (0.05 * degree(iG)),
     edge.arrow.size = 0,
     edge.width = 2,
     vertex.label = toupper(V(iG)$name))
par(oPar)


# == Components

# The igraph function components() tells us whether there are components of the
# graph in which there is no path to other components.
components(iG)

# In the _membership_ vector, nodes are annotetd with the index of the component
# they are part of. Sui7 is the only node of component 2, Cyj1 is in the third
# component etc. This is perhaps more clear if we sort by component index
sort(components(iG)$membership)

# Retrieving e.g. the members of the first component from the list can be done by subsetting:

                                components(iG)$membership == 1  # logical ..
      components(iG)$membership[components(iG)$membership == 1]
names(components(iG)$membership)[components(iG)$membership == 1]



# == RANDOM GRAPHS AND GRAPH METRICS =================================


# Let's explore some of the more interesting, topological graph measures from
# the Pavlopoulos paper. We start by building a somewhat bigger graph. We aren't
# quite sure whether biological graphs are small-world, or random-geometric, or
# preferential-attachment ... but igraph has ways to simulate the basic ones
# (and we could easily simulate our own). Look at the following help pages:

?sample_gnm                      # see also sample_gnp for the Erdös-Rényi models
?sample_smallworld               # for the Watts & Strogatz model
?sample_pa                       # for the Barabasi-Albert model

# But note that there are many more sample_ functions. Check out the docs!

# sample_pa() is the stochastic algorithm for a scale-free that we already used
# above. It is based on the Barabasi-Albert preferential attachment model. Let's
# plot a moderately large graph, with a force-directed layout, colored by
# node-degree."


set.seed(27172)
gBA <- sample_pa(1000, power=1.05, out.dist = c(0, 1, 0.3, 0.09, 0.027, 0.0081))
gBAxy <- layout_with_graphopt(gBA, charge=0.001)   # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(gBA,
     layout = gBAxy,
     rescale = FALSE,
     xlim = c(min(gBAxy[,1]), max(gBAxy[,1])),
     ylim = c(min(gBAxy[,2]), max(gBAxy[,2])),
     vertex.color=heat.colors(max(degree(gBA)+1))[degree(gBA)+1],
     vertex.size = 800 + (80 * degree(gBA)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# The large nodes arise from the fact that a node of high degree is more likely
# to attract even more additional edges. "


# == Diameter

diameter(gBA)  # The maximum length shortest path
# let's plot this path though:
lines(gBAxy[get_diameter(gBA),], lwd=3, col="#00AA00")

# ... a good reminder that as nice as the network plots are, they don't
# necessarily show us the true topological structure of a graph well."


# == Centralization scores

?centralize
bCgBA <- centr_betw(gBA)  # calculate betweenness centrality

# replot or graph, and color by log_betweenness"

nodeBetw <- bCgBA$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(gBA,
     layout = gBAxy,
     rescale = FALSE,
     xlim = c(min(gBAxy[,1]), max(gBAxy[,1])),
     ylim = c(min(gBAxy[,2]), max(gBAxy[,2])),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 800 + (80 * degree(gBA)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)

# Note that the betweenness - the number of shortest paths that pass through a
# node, is in general higher for high-degree nodes - but not always: one of the
# larger nodes has a very low betweenness and is colored red: this measure
# really depends on the detailed local topology of the graph."

# Lets plot the same thing for our random geometric graph:

bCiGRG <- centr_betw(iGRG)  # calculate betweenness centrality

nodeBetw <- bCiGRG$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)




# == CLUSTERING ======================================================

# Clustering finds "communities" in graphs - and depending what the edges
# represent, these could be complexes, pathways, biological systems or similar.
# There are many graph-clustering algorithms. One approach with many attractive
# properties is the Map Equation, developed by Martin Rosvall. See:
# http://www.ncbi.nlm.nih.gov/pubmed/18216267 and htttp://www.mapequation.org


iGRGclusters <- cluster_infomap(iGRG)
modularity(iGRGclusters)
membership(iGRGclusters)
table(membership(iGRGclusters))


gBAclusters <- cluster_infomap(gBA)
modularity(gBAclusters)   # ... measures how separated the different membership
                          # types are from each other
membership(gBAclusters)
table(membership(gBAclusters))  # The sizes of the clusters

# Lets plot our graph again, coloring the nodes of the first five communities by
# their cluster membership:

# first, make a vector with as many grey colors as we have communities ...
commColors <- rep("#f1eef6", max(membership(iGRGclusters)))
# ... the overwrite the first five with "real colors" - something like rust,
# lilac, pink, and mauve or so.
commColors[1:5] <- c("#980043", "#dd1c77", "#df65b0", "#c994c7", "#d4b9da")


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=commColors[membership(iGRGclusters)],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)




# ==============================================================================
#        PART FOUR: PLAY AND EXPLORE
# ==============================================================================

# In order for you to explore some real, biological networks, I give you a
# dataframe of functional relationships of human proteins that I have downloaded from the
# STRING database. The full table has 8.5 million records, here is a subset of
# records with combined confidence scores > 980

# The selected set of edges with a confidence of > 980 is a dataframe with about
# 50,000 edges and 6,500 unique proteins. You can load the saved dataframe here
# (and also read more about what the numbers mean at
# http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

load("STRINGedges.RData")

head(STRINGedges)


# make a graph from this dataframe
?graph_from_data_frame

gSTR <- graph_from_data_frame(STRINGedges)

compSTR <- components(gSTR)
summary(compSTR) # our graph is fully connected!


hist(log(degree(gSTR)))
dDegSTR <- degree_distribution(gSTR, cumulative=TRUE)[-1]
plot(log(degree(gSTR)), log(dDegSTR[degree(gSTR)]))

# This is very cool! What does it mean?

# Now explore some more:

# === CLIQUES   ========
# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.
cliques <- largest_cliques(gSTR)

cliques[[1]]

# Pick one of the proteins and find out what it is (you can simply Google for
# it). Is this expected?

# === BETWEENNESS CENTRALITY   =======================================

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- centr_betw(gSTR)

# remember: bC$res contains the results
head(BC$res)

BC$res[1]   # betweeness centrality of node 1 in the graph ...
            # ... which one is node 1?
V(gSTR)[1]

# to get the ten-highest nodes, we simply we define the element indices of BC to be the names ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)
head(sBC)

# This ordered vector means: node 3,862 has the highest betweeness centrailty, node 1,720 has the second highest.. etc.

BCsel <- as.numeric(names(sBC)[1:10])
BCsel

# Now we can get the IDs...
ENSPsel <- names(V(gSTR)[BCsel])

# We are going to use these IDs to produce some output for you to print out and bring to class, so I need you to personalize ENSPsel part with the following two lines of code:
set.seed(myStudentNumber)
ENSPsel <- sample(ENSPsel)

# Here are the IDs:
ENSPsel

#  Next, to find what these proteins are...

# We could now Google for all of them to learn more. But really, that would be
# lame. Let's instead use the very, very useful biomaRt package to translate
# these Ensemble IDs into gene symbols.

# == biomaRt =========================================================

# IDs are just labels, but for _bio_informatics we need to learn more about the
# biological function of the genes or proteins that we retrieve via graph data
# mining. biomaRt is the tool of choice. It's a package distributed by the
# bioconductor project. This here is not a biomaRt tutorial (that's for another
# day), simply a few lines of sample code to get you started on the specific use
# case of retrieving descriptions for ensembl protein IDs.


if (!require(biomaRt)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    library("biomaRt")
}

# define which dataset to use ...
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# what filters are defined?
filters <- listFilters(myMart)
filters

# and what attributes can we filter for?
attributes <- listAttributes(myMart)
attributes

# Now - lets look for the correct name of filters that are useful for ENSP IDs ...
filters[grep("ENSP", filters$description), ]

# ... and the correct attribute names for gene symbols and descriptions ...
attributes[grep("symbol", attributes$description, ignore.case=TRUE), ]
attributes[grep("description", attributes$description, ignore.case=TRUE), ]


# ... so we can put this together: here is a syntax example:
getBM(filters = "ensembl_peptide_id",
      attributes = c("hgnc_symbol",
                     "wikigene_description",
                     "interpro_description",
                     "phenotype_description"),
      values = "ENSP00000000442",
      mart = myMart)

# A simple loop will now get us the information for our 10 most central genes
# from the human subset of STRING.

# But first, we need to remove the string "9606." from the ID:
ENSPsel <- gsub("9606\\.", "", ENSPsel)

CPdefs <- list()  # since we don't know how many matches one of our queries
                  # will return, we'll put the result dataframes into a list.

for (ID in ENSPsel) {
    CPdefs[[ID]] <- getBM(filters = "ensembl_peptide_id",
                          attributes = c("hgnc_symbol",
                                         "wikigene_description",
                                         "interpro_description",
                                         "phenotype_description"),
                          values = ID,
                          mart = myMart)
}

# So what are the proteins with the ten highest betweenness centralities?
#  ... are you surprised? (I am! Really.)

# Final task: Write a loop that will go through your list and print the first
# column's symbol and wikigene description. (Hint, you can structure your loop
# in the same way as the loop that created CPdefs in the frist place. )

# Print the R code for your loop and its output onto a sheet of paper, write
# your student number and name on it, and bring this to class.


# [END]
