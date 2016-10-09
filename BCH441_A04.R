# BCH441_A04.R
#
# Purpose:  Supporting scripts for BCH441 (Bioinformatics) at the University of
# Toronto, Fall 2016 - Assignment 04
#
# Version: 1.2
#
# Date:    2016  10 Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.2    Sections up to APSES domain annotation complete
# V 1.0    First code
#
# TODO:
#
#
# ==============================================================================

# DO NOT SIMPLY  source()  THIS FILE!

# ==============================================================================
#        PART ONE: Database maintenance
# ==============================================================================

# === Merging databases ========================================================

# From time to time, we will want to update the systems database - either
# because the schema has changed, e.g. we may include additional tables, or
# because I have added reference data: sequences, annotations and the like. The
# goal here is to make this as simple as possible, keeping all information you
# may have entered intact, and avoiding to duplicate entries. This is not quite
# as easy as one would hope. First of all, since we need to merge records from
# different databases - your own, and the reference - we can't assume that the
# primary keys are unique. You may have introduced the same proteinID that I
# used later for the reference database. Thus we need to update keys. But if we
# update primary keys, we need to make sure that every record that references a
# key as a foreign key is updated as well. And ensuring database integrity for
# "cascading updates" takes a _lot_ of care to guarantee correctness. You don't
# want to read through that code. I don't want to write it. And we should not
# end up pretending to build a real and proper database engine here!
#
# We'll do something simple instead. The sanest way is to create
# primary keys with "namespaces", ie. the keys will reflect the source of the
# information. This way, even though the databases are merged, the information
# sources can be identified. Still, it can be accessed through a common
# mechanism. A database merge with separate namespaces can simply proceed as
# follows:

#
#    for each table in the database
#        merge the table from both source databases
#        throw out any row with a duplicate ID
#
# In this way, if we merge your working database with the reference database
# all reference database records will be updated, and all your working records
# will remain.

# And while we are adding semnatics to keys, we should also identify in which
# table they are being used. You may remember that I told you: keys should not
# try to express semantics? Here we have a slightly different situation: since
# we are not actually constructing an underlying database engine, we will often
# encounter keys that we use "manually". And in that case it is very helpful to
# make the key's domain explicit.

# The code I wrote to implement this was loaded when you typed init().

# Let's have a brief look at the new dbAutoincrement() function. You'll see the
# code as usual when you execute the function name without its parentheses:

dbAutoincrement

# The default namespace is "my" - that's the namespace you use. I will use
# "ref" for all entries of the reference database. Thus you normally won't
# need to enter the namespace identifier - only I do that.

refDB$protein$ID
dbAutoincrement(refDB$protein$ID)
dbAutoincrement(refDB$protein$ID, ns = "ref")

# Once the namespaces are thus kept well apart, we won't overwrite each other's
# work. Merging two table becomes a breeze:

tableUnion

# Here is an example
this <- data.frame(ID = "ref_1", type = "one crow", stringsAsFactors = FALSE)
this <- rbind(this, data.frame(ID = dbAutoincrement(this$ID, ns="ref"),
                               type = "for",
                               stringsAsFactors = FALSE))
this <- rbind(this, data.frame(ID = dbAutoincrement(this$ID, ns="ref"),
                               type = "sorrow",
                               stringsAsFactors = FALSE))
that <- data.frame(ID = "my_1", type = "two crows", stringsAsFactors = FALSE)
that <- rbind(that, this[2:3, ])  # repeat some rows
that <- rbind(that, data.frame(ID = dbAutoincrement(that$ID),
                               type = "for ",
                               stringsAsFactors = FALSE))
that <- rbind(that, data.frame(ID = dbAutoincrement(that$ID),
                               type = "joy ...",
                               stringsAsFactors = FALSE))

that
this

rhyme <- tableUnion(this, that)
rhyme$type

# Finally, we need to do this for all tables in the datamodel. Well, almost all:
# we introduce a field $version, through which we can ensure that the datamodels
# have the same schema. If they would have different schemas the merging would
# break. But we don't merge this single field. (Incidentally: a schema version
# is not "data", rather we call this "metadata": information _about_ the data.)
# So we exclude the "version" from the names of elements to merge.
#
# But how to iterate over all tables in a list by name? You might be tempted to
# try something like
#
# n <- names(db)
# myDB$n[1], myDB$n[2] etc.
# ... but NO! That's not how the $ operator works.
#
# The R community keeps particularly poignant comments from the R-help mailing
# list in a package called "fortunes", and fortune(312) reads:
#
# " The problem here is that the $ notation is a magical shortcut and like any
#   other magic if used incorrectly is likely to do the programmatic equivalent
#   of turning yourself into a toad.
#
# -- Greg Snow (in response to a user that wanted
#    to access a column whose name is stored in y
#    via x$y rather than x[[y]])
#    R-help (February 2012)

# So in order to avoid turning into toads, we create a vector "tables", iterate
# over "table" elements from "tables" and use them as ref[[table]] ... etc.

# === Putting data into the new schema =========================================

# Entering the YFO data into the new schema takes almost exactly
# the same syntax as the code you wrote for the last assignment.
#

# === TASK: ===
# Initialize myDB with the YFO homologue of yeast Mbp1

# First: initialize "myDB", this is the name of the R object that
#        will contain data you collect.

myDB <- dbInit()   # This creates an empty database with the latest schema

# Next: Add your protein and the YFO to the database. Copy the code-template
#       below to your myCode.R file, edit it to replace the placeholder
#       items with your data, and execute the code. Then continue below
#       the code template.


# ===== begin code template: add a protein and an organism to the database =====

# == edit placeholder items
myBinomial <- "BINOMIAL NAME"
myTaxonomyId <- as.integer(TAX_ID)
myProteinName <- "PROTEIN NAME"
myProteinRefSeqID <- "REFSEQID"
myProteinUniProtID <- "UNIPROTID"
mySeq <- "
SEQUENCE
"

# == create the protein entry
proteinRow <- data.frame(ID = dbAutoincrement(myDB$protein$ID),
                    name = myProteinName,
                    RefSeqID = myProteinRefSeqID,
                    UniProtID = myProteinUniProtID,
                    taxonomy.ID = myTaxonomyId,
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
myDB$protein <- rbind(myDB$protein, proteinRow)

# == create the taxonomy entry
taxonomyRow <- data.frame(ID = myTaxonomyId,
                    species = myBinomial,
                    stringsAsFactors = FALSE)
myDB$taxonomy <- rbind(myDB$taxonomy, taxonomyRow)
# ===== end code template ===========================================

# ... continue here.

# myDB now contains one record each in two tables. The remaining tables exist
# but they are empty.

# Now let's merge myDB with the data from refDB. refDB should already have been
# loaded from .utilities.R ; you can also explore the original script with which
# refDB was created, for your reference: it is create_refDB.R  The database
# refDB is the default argument for dbMerge(), so you don't need to
# specify it. By assigning the result of dbMerge() back to myDB we overwrite the
# previous version.

myDB <- dbMerge(myDB)

str(myDB)

# check the protein table
View(myDB$protein[ , c("ID", "name", "RefSeqID")])

# Let's compute sequence lengths on the fly (with the function nchar() ) and
# add them to our view. Then we'll open this with the table viewer function
# View()

View(cbind(myDB$protein[ , c("ID", "name", "RefSeqID")],
           length = nchar(myDB$protein$sequence)))

# Where does your protein's length fall relative to the reference proteins?
# About the same? Much shorter? Much longer?

# Is that the right sequence?
myDB$protein$sequence[myDB$protein$ID == "my_pro_1"]

# If not, don't continue! Fix the problem first.
# Let me repeat: If this does not give you the right sequence of the YFO
#                Mbp1 homologue, DO NOT CONTINUE. Fix the problem.

# Is that the right taxonomy ID and binomial name for YFO?
myDB$taxonomy[myDB$taxonomy$species == YFO, ]

# If not, or if the result was "<0 rows> ... " then DO NOT CONTINUE.
# Fix the problem first.
#
# === Saving and loading data ==================================================
#
# There are many ways to save data to a file on disk and read it back in. One of
# the most convenient is the function pair save() and load().
#
# save() saves an R object to a file. Its signature is

#    save(<object-name/s>, file = <file-name>)

# The object, or objects that are saved like this are identically recreated when
# load() is executed. The object name is not changed - you don't assign the
# result of load(). You use load() for its "side-effect": re-creating the saved
# object, not for using the return-value. Therefore the signature for load() is
# simply

#     load(<file-name>)

# All R objects in <file-name> are created by load(), if they don't yet exist.
# If they already exist, the will be overwritten. The only effect you see is
# that the object appears in the Environment pane of RStudio (if it wasn't there
# already), or you may notice that its content has changed if it did exist.

# == TASK: ==
# Save myDB so you can recreate it easily when you next open RStudio.

save(myDB, file = "myDB.01.RData")  # Note that I give this a version number!

# Let's confirm:
rm(myDB)
nrow(myDB$protein)  # must give an error
load("myDB.01.RData")
nrow(myDB$protein)  # must be 11. If not, don't continue etc.

# Now you can return to the Assignment on the Course Wiki.



#
# ==============================================================================
#        PART TWO: Dotplot and MDM
# ==============================================================================

# If you come here after restarting RStudio, you need to reload the database.
# load("myDB.01.RData")

# First, we install and load the Biostrings package.
if (!require(Biostrings, quietly=TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    library(Biostrings)
}


# This is a large collection of tools ...
help(package = "Biostrings")

# ... with mutation matrices and other useful datasets
data(package = "Biostrings")

# Let's load BLOSUM62
data(BLOSUM62)

# ... and see what it contains. (You've seen this before, right?)
BLOSUM62

# We can simply access values via the row/column names to look at the data
# for the questions I asked in the Assignment on the Wiki:
BLOSUM62["H", "H"]
BLOSUM62["S", "S"]

BLOSUM62["L", "K"]
BLOSUM62["L", "I"]


BLOSUM62["R", "W"]
BLOSUM62["W", "R"]   # the matrix is symmetric!

# Now let's craft code for a dotplot. That's surprisingly simple. We build a
# matrix that has as many rows as one sequence, as many columns as another. Then
# we go through every cell of the matrix and enter the pairscore we encounter
# for the amino acid pair whose position corresponds to the row and column
# index. Finally we visualize the matrix in a plot.
#

# First we fetch our sequences and split them into single characters.
MBP1_SACCE <- s2c(myDB$protein$sequence[myDB$protein$name == "MBP1_SACCE"])
my <- paste("MBP1_", biCode(YFO), sep = "")
MBP1_YFO <- s2c(myDB$protein$sequence[myDB$protein$name == my])

# Check that we have two character vectors of the expected length.
str(MBP1_SACCE)
str(MBP1_YFO)

# How do we get the pairscore values? Consider: a single pair of amino acids can
# be obtained from sequence SACCE and YFO eg. from position 13 and 21 ...
MBP1_SACCE[13]
MBP1_YFO[21]

# ... using these as subsetting expressions, we can pull the pairscore
# from the MDM
BLOSUM62[MBP1_SACCE[13], MBP1_YFO[21]]

# First we build an empty matrix that will hold all pairscores ...
dotMat <- matrix(numeric(length(MBP1_SACCE) * length(MBP1_YFO)),
                 nrow = length(MBP1_SACCE), ncol = length(MBP1_YFO))

# ... then we loop over the sequences and store the scores in the matrix.
#
for (i in 1:length(MBP1_SACCE)) {
    for (j in 1:length(MBP1_YFO)) {
        dotMat[i, j] <- BLOSUM62[MBP1_SACCE[i], MBP1_YFO[j]]
    }
}

# Even though this is a large matrix, this does not take much time ...
# Let's have a look at a small block of the values:

dotMat[1:10, 1:10]

# Rows in this matrix correspond to an amino acid from MBP1_SACCE, columns in
# the matrix correspond to an amino acid from MBP1_YFO.

# To plot this, we use the image() function. Here, with default parameters.

image(dotMat)

# Be patient, this takes a few moments to render: more than 500,000 values.
# Nice.
# What do you expect?
# What would similar sequences look like?
# What do you see?

#You migh notice a thin line of yellow along the diagonal, moving approximately
# from bottom left to top right, fading in and out of existence. This is the
# signature of extended sequence similarity.

# Let's magnify this a bit by looking at only the first 200 amino acids ...
image(dotMat[1:200, 1:200])

# ... and, according to our normal writing convention, we would like the
# diagonal to run from top-left to bottom-right since we write from left to
# right and from top to bottom...
image(dotMat[1:200, 1:200], ylim = 1.0:0.0)

# ... and we would like the range of the x- and y- axis to correspond to the
# sequence position ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1))

# ... and labels! Axis labels would be nice ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1),
      xlab = "MBP1_YFO", ylab = "MBP1_SACCE" )

# ... and why don't we have axis-numbers on all four sides? Go, make that right
# too ...
len <- 200
image(x = 1:len, y = 1:len,  dotMat[1:len, 1:len], ylim=c(len,1),
      xlab = "MBP1_YFO", ylab = "MBP1_SACCE", axes = FALSE)
box()
axis(1, at = c(1, seq(10, len, by=10)))
axis(2, at = c(1, seq(10, len, by=10)))
axis(3, at = c(1, seq(10, len, by=10)))
axis(4, at = c(1, seq(10, len, by=10)))

# ... you get the idea, we can infinitely customize our plot. However a good way
# to do this is to develop a particular view for, say, a report or publication
# in a script and then put it into a function. I have put a function into the
# utilities file and called it dotPlot2(). Why not dotPlot() ... that's because
# there already is a dotplot function in the seqinr package:

dotPlot(MBP1_SACCE, MBP1_YFO)                                 # seqinr
dotPlot2(MBP1_SACCE, MBP1_YFO, xlab = "SACCE", ylab = "YFO")  # Our's

# Which one do you prefer? You can probably see the block patterns that arise
# from segments of repetitive, low complexity sequence. But you probably have to
# look very closely to discern the faint diagonals that correspond to similar
# sequence.


# Let's see if we can enhance the contrast between distributed noise and the
# actual alignment of conserved residues. We can filter the dot matrix with a
# pattern that enhances diagonally repeated values. Every value in the matrix
# will be replaced by a weighted average of its neighborhood. Here is  a
# diagonal-filter:

myFilter <- matrix(numeric(25), nrow = 5)
myFilter[1, ] <- c( 1, 0, 0, 0, 0)
myFilter[2, ] <- c( 0, 1, 0, 0, 0)
myFilter[3, ] <- c( 0, 0, 1, 0, 0)
myFilter[4, ] <- c( 0, 0, 0, 1, 0)
myFilter[5, ] <- c( 0, 0, 0, 0, 1)

# I have added the option to read such filters (or others that you could define on your own) as a parameter of the function.

dotPlot2(MBP1_SACCE, MBP1_YFO, xlab = "SACCE", ylab = "YFO", f = myFilter)

# I think the result shows quite nicely how the two sequences are globally
# related and where the regions of sequence similarity are. Play with this a bit
# ...  Can you come up with a better filter? If so, eMail us.

# Back to the Course Wiki.



# ==============================================================================
#        PART THREE: Biostrings Pairwise Alignment
# ==============================================================================

# Biostrings is one of the basic packages that the Bioconductor software
# landscape builds on. It stores sequences in "AAstring" objects and these are
# complex software structures that are designed to be able to handle
# genome-scale sequences. Biostrings functions - such as the alignment functions
# - expect their input to be Biostrings objects.
?AAString

AAString("ACDE")
s <- AAString("ACDE")
str(s)
# See: it's complicated. This is an "S4" object. Bioconductor uses these objects
# almost exclusively, but we will not be poking around in their internals. Just
# this: how do we get the sequence back out of an AAString object? The help page
# for XString - the parent "class" of AAStrings - mentions the  alternatives:

as.character(s)  # the base R version
toString(s)      # using the Biostrings function toString()

# While we need to remember to convert our sequences from the character vectors
# that we store in our database, to AAStrings that we can align, the alignment
# itself is really straightforward. The pairwiseAlignment() function was written
# to behave exactly like the functions you encountered on the EMBOSS server.

# First: make AAString objects ...
aaMBP1_SACCE <- AAString(myDB$protein$sequence[myDB$protein$name ==
                                                   "MBP1_SACCE"])
my <- paste("MBP1_", biCode(YFO), sep = "")
aaMBP1_YFO <-   AAString(myDB$protein$sequence[myDB$protein$name ==
                                                   my])

?pairwiseAlignment

# ... and align.
# Global optimal alignment with end-gap penalties is default. (like EMBOSS needle)
ali1 <-  pairwiseAlignment(
    aaMBP1_SACCE,
    aaMBP1_YFO,
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5)

str(ali1)  # Did you think the AAString object was complicated ?

# This is a Biostrings alignment object. But we can use Biostrings functions to
# tame it:
ali1
writePairwiseAlignments(ali1)   # That should look familiar

# And we can make the internal structure work for us  (@ is for classes as
# $ is for lists ...)
str(ali1@pattern)
ali1@pattern
ali1@pattern@range
ali1@pattern@indel
ali1@pattern@mismatch

# or work with "normal" R functions
# the alignment length
nchar(ali1@pattern)

# the number of identities
sum(s2c(as.character(ali1@pattern)) ==
    s2c(as.character(ali1@subject)))

# ... e.g. to calculate the percentage of identities
    100 *
    sum(s2c(as.character(ali1@pattern)) ==
          s2c(as.character(ali1@subject))) /
    nchar(ali1@pattern)
# ... which should be the same as reported in the writePairwiseAlignments()
# output. Awkward to type? Then it calls for a function:
#
percentID <- function(al) {
    # returns the percent-identity of a Biostrings alignment object
    return(100 *
           sum(s2c(as.character(al@pattern)) ==
               s2c(as.character(al@subject))) /
           nchar(al@pattern))
}

percentID(ali1)

# Compare with local optimal alignment (like EMBOSS Water)
ali2 <-  pairwiseAlignment(
    aaMBP1_SACCE,
    aaMBP1_YFO,
    type = "local",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 50,
    gapExtension = 10)

writePairwiseAlignments(ali2)   # This has probably only aligned the N-terminal
                                # DNA binding domain - but that one has quite
                                # high sequence identity:
percentID(ali2)

# == TASK: ==

# Compare the two alignments. I have weighted the local alignment heavily
# towards an ungapped alignment by setting very high gap penalties. Try changing
# the gap penalties and see what happens: how does the number of indels change,
# how does the length of indels change...

# Fine. Please return to the Wiki to study BLAST alignment...


# ==============================================================================
#        PART FOUR: APSES Domain annotation by alignment
# ==============================================================================

# In this section we define the YFO APSES sequence by performing a global,
# optimal sequence alignment of the yeast domain with the full length protein
# sequence of the protein that was the most similar to the yeast APSES domain.
#

# I have annotated the yeast APSES domain as a proteinAnnotation in the
# database. To view the annotation, we can retrieve it via the proteinID and
# featureID. Here is the yeast protein ID:
myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]
proID <- myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]

# ... and if you look at the feature table, you can identify the feature ID
myDB$feature[ , c("ID", "name", "description")]
myDB$feature$ID[myDB$feature$name == "APSES fold"]
ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"]

# ... and with the two annotations we can pull the entry from the protein
# annotation table
myDB$proteinAnnotation[myDB$proteinAnnotation$protein.ID == proID &
                           myDB$proteinAnnotation$feature.ID == ftrID, ]

myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                           myDB$proteinAnnotation$feature.ID == ftrID]

fanID <- myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                                myDB$proteinAnnotation$feature.ID == ftrID]

# The annotation record contains the start and end coordinates which we can use
# to define the APSES domain sequence with a substr() expression.
substr(myDB$protein$sequence[myDB$protein$ID == proID],
       myDB$proteinAnnotation$start[myDB$proteinAnnotation$ID == fanID],
       myDB$proteinAnnotation$end[myDB$proteinAnnotation$ID == fanID])

# Lots of code. But don't get lost. Let's recapitulate what we have done: we have selected from the sequence column of the protein table the sequence whose name is "MBP1_SACCE", and selected from the proteinAnnotation table the start and end coordinates of the annotation that joins an "APSES fold" feature with the sequence, and used the start and end coordinates to extract a substring. The expressions get lengthy, but it's not hard to wrap all of this into a function so that we only need to define name and feature.

dbGetFeatureSequence
dbGetFeatureSequence(myDB, "MBP1_SACCE", "APSES fold")


# Let's convert this to an AAstring and assign it:
aaMB1_SACCE_APSES <- AAString(dbGetFeatureSequence(myDB,
                                                   "MBP1_SACCE",
                                                   "APSES fold"))

# To align, we need the YFO sequence. Here is it's definition again, just
# in case ...

my <- paste("MBP1_", biCode(YFO), sep = "")
aaMBP1_YFO <-   AAString(myDB$protein$sequence[myDB$protein$name ==
                                                   my])

# Now let's align these two sequences of very different length without end-gap
# penalties using the "overlap" type. "overlap" turns the
# end-gap penalties off and that is crucially important since
# the sequences have very different length.

aliApses <-  pairwiseAlignment(
    aaMB1_SACCE_APSES,
    aaMBP1_YFO,
    type = "overlap",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5)

# Inspect the result. The aligned sequences should be clearly
# homologous, and have (almost) no indels. The entire "pattern"
# sequence from QIYSAR ... to ... KPLFDF  should be matched
# with the "query". Is this correct?
writePairwiseAlignments(aliApses)

# If this is correct, you can extract the matched sequence from
# the alignment object. The syntax is a bit different from what
# you have seen before: this is an "S4 object", not a list. No
# worries: as.character() returns a normal string.
as.character(aliApses@subject)

# Now, what are the aligned start and end coordinates? You can read them from
# the output of writePairwiseAlignments(), or you can get them from the range of
# the match.

str(aliApses@subject@range)

# start is:
aliApses@subject@range@start

# ... and end is:
aliApses@subject@range@start + aliApses@subject@range@width - 1

# Since we have this section defined now, we can create a feature annotation
# right away and store it in myDB.  Copy the code-template below to your
# myCode.R file, edit it to replace the placeholder items with your data:
#
#  - The "PROTEIN ID" is to be replaced with the ID of MBP1_YFO
#  - The "FEATURE ID" is to be replaced with the ID of "APSES fold"
#  - START and END are to be replaced with the coordinates you geot above
#
# Then execute the code and continue below the code template. If you make an
# error, there are instructions on how to recover, below.
#
# ===== begin code template: add a proteinAnnotation to the database =====

# == edit placeholder items
myProteinID <- "PROTEIN ID"
myFeatureID <- "FEATURE ID"
myStart <- START
myEnd   <- END

# == create the proteinAnnotation entry
panRow <- data.frame(ID = dbAutoincrement(myDB$proteinAnnotation$ID),
                         name = myProteinName,
                         protein.ID = myProteinID,
                         feature.ID = myFeatureID,
                         start = myStart,
                         end = myEnd,
                         stringsAsFactors = FALSE)
myDB$proteinAnnotation <- rbind(myDB$proteinAnnotation, panRow)

# == check that this was successful and has the right data
myDB$proteinAnnotation[nrow(myDB$proteinAnnotation), ]

# ===== end code template ===========================================

# ... continue here.
# I expect that a correct result would look something like
#          ID protein.ID feature.ID start end
# 63 my_fan_1   my_pro_1  ref_ftr_1     6 104

# If you made a mistake, simply overwrite the current version of myDB by loading
# your saved, good version:  load("myDB.01.RData") and correct your mistake

# If this is correct, save it
save(myDB, file = "myDB.02.RData")  # Note that it gets a new version number!

# Done with this part. Copy the sequence of the APSES domain of MBP!_YFO - you
# need it for the reverse BLAST search, and return to the course Wiki.


# ==============================================================================
#        PART FIVE: SMART Domain annotations
# ==============================================================================


# ... TBC

# ==============================================================================
#        PART SIX: Multiple sequence alignment
# ==============================================================================

# ... TBC




# [END]
