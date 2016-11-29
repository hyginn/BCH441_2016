# BCH441_A03.R
#
# Purpose:  Supporting scripts for BCH441 (Bioinformatics)
#           at the University of Toronto, Fall 2016
#           Assignment 03
#
# Version: 1.2
#
# Date:    2016  09 Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.2    Fixed more typos
# V 1.1    Fixed typo
# V 1.0    First code
#
# TODO:
#
#
# ==============================================================================

# DO NOT SIMPLY  source()  THIS FILE!

# ==============================================================================
#        PART ONE: YFO Species
# ==============================================================================

# A detailed description of the process of compiling the YFO list of genome
# sequenced fungi with protein annotations and Mbp1 homologues is in the file
# BCH441_A03_makeYFOlist.R I encourage you to study it. Here we load the
# resulting vector of species name, then pick one of them in a random but
# reproducible way, determined by your student number.

load("data/YFOspecies.RData")  # load the species names
set.seed(myStudentNumber)      # seed the random number generator
YFO <- sample(YFOspecies, 1)   # pick a species at random
# write the result to your personalized profile data so we can use the result in
# other functions
cat(sprintf("YFO <- \"%s\"\n", YFO), file = ".myProfile.R", append = TRUE)

YFO         # so, which species is it ... ?
biCode(YFO) # and what is it's "BiCode" ... ?
            # (note these on your Student Wiki Page).


# ==============================================================================
#        PART TWO: Implementing the (protein part of) the datamodel
# ==============================================================================

# Every entity in our datmodel can be implemented by a dataframe. To keep things
# organized, we create a list, that contains the enity tables. Here is a
# definition for a dataframe that holds the protein data for yeast Mbp1 as
# described in the data model schema:
#
# === Creating two tables ===
db <- list()      # we'll call our database "db" and start with an empty list

# Then add a dataframe with one row to the list
db$protein <- data.frame(ID = as.integer(1),
                         name = "Mbp1",
                         RefSeqID = "NP_010227",
                         UniProtID = "P39678",
                         taxonomy.ID = as.integer(4932),
                         sequence = "...",              # just a placeholder
                         stringsAsFactors = FALSE)

str(db)

# Next, we create the taxonomy table and add it into db
db$taxonomy <- data.frame(ID = as.integer(4932),
                          species = "Saccharomyces cerevisiae",
                          stringsAsFactors = FALSE)

# Let's add a second protein: We create a one-row dataframe with the data, then we rbind() it to the existing data frame:
myRow <- data.frame(ID = as.integer(2),
                    name = "Res2",
                    RefSeqID = "NP_593032",
                    UniProtID = "P41412",
                    taxonomy.ID = as.integer(4896),
                    sequence = "...",              # again, just a placeholder
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(4896),
                    species = "Schizosaccharomyces pombe",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

str(db)

# you can look at the contents of the tables in the usual way we would access
# elements from lists and dataframes. Here are some examples:
db$protein
db$protein$RefSeqID
db$protein[,"name"]
db$taxonomy
db$taxonomy$species
biCode(db$taxonomy$species)

# Here is an example to look up information in one table,
# based on a condition in another table:
# what is the species name for the protein
# whose name is "Mbp1"?

# First, get the taxonomy.ID for the Mbp1 protein. This is
# the key we need for the taxonomy table. We find it in a cell in the
# table: db$protein[<row>, <column>]
# <row> is that row for which the value in
# the "name" column is Mbp1:

db$protein$name == "Mbp1" # TRUE FALSE

# The <column> is called "taxonomy.ID". Simply
# insert these two expressions in the square
# brackets.

db$protein[db$protein$name == "Mbp1", "taxonomy.ID"]

# Assign the taxonomy.ID value ...
x <- db$protein[db$protein$name == "Mbp1", "taxonomy.ID"]

# ... and fetch the species_name value from db$taxonomy

db$taxonomy[db$taxonomy$ID == x, "species"]


# === An excursion into regular expressions ====================================

# The Mbp1 sequence as copied from the NCBI Website
mySeq <- "
  1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
 61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
181 klgvnlqrsq sdmgfprpai pnssisttql psirstmgpq sptlgileee rhdsrqqqpq
241 qnnsaqfkei dledglssdv epsqqlqqvf nqntgfvpqq qssliqtqqt esmatsvsss
301 pslptspgdf adsnpfeerf pgggtspiis miprypvtsr pqtsdindkv nkylsklvdy
361 fisnemksnk slpqvllhpp phsapyidap idpelhtafh wacsmgnlpi aealyeagts
421 irstnsqgqt plmrsslfhn sytrrtfpri fqllhetvfd idsqsqtvih hivkrksttp
481 savyyldvvl skikdfspqy rielllntqd kngdtalhia skngdvvffn tlvkmgaltt
541 isnkegltan eimnqqyeqm miqngtnqhv nssntdlnih vntnnietkn dvnsmvimsp
601 vspsdyityp sqiatnisrn ipnvvnsmkq masiyndlhe qhdneikslq ktlksisktk
661 iqvslktlev lkesskdeng eaqtnddfei lsrlqeqntk klrkrliryk rlikqkleyr
721 qtvllnklie detqattnnt vekdnntler lelaqeltml qlqrknklss lvkkfednak
781 ihkyrriire gtemnieevd ssldvilqtl iannnknkga eqiitisnan sha
//
"

mySeq     # "\n" means: line-break

mySeq <- gsub("[^a-zA-Z]", "", mySeq) # replace all non-letters with ""

mySeq

# Now to change the sequence to upper-case. R has toupper()
# and tolower().

toupper(mySeq)

# I have combined this into a function in the .utilities.R script. Execute:

dbSanitizeSequence

# (Executing a function name without parentheses displays the function
# code.)

# try the function:

test <- "f123  a$^&&*)m. i	l马马虎虎 é yßv+w"
dbSanitizeSequence(test)
dbSanitizeSequence(mySeq)


# === Updating the database ====================================================

# Modifying a field

# Here is code that modifies the sequence field in the protein table of the
# database:
#
mySeq <- "
  1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
 61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
181 klgvnlqrsq sdmgfprpai pnssisttql psirstmgpq sptlgileee rhdsrqqqpq
241 qnnsaqfkei dledglssdv epsqqlqqvf nqntgfvpqq qssliqtqqt esmatsvsss
301 pslptspgdf adsnpfeerf pgggtspiis miprypvtsr pqtsdindkv nkylsklvdy
361 fisnemksnk slpqvllhpp phsapyidap idpelhtafh wacsmgnlpi aealyeagts
421 irstnsqgqt plmrsslfhn sytrrtfpri fqllhetvfd idsqsqtvih hivkrksttp
481 savyyldvvl skikdfspqy rielllntqd kngdtalhia skngdvvffn tlvkmgaltt
541 isnkegltan eimnqqyeqm miqngtnqhv nssntdlnih vntnnietkn dvnsmvimsp
601 vspsdyityp sqiatnisrn ipnvvnsmkq masiyndlhe qhdneikslq ktlksisktk
661 iqvslktlev lkesskdeng eaqtnddfei lsrlqeqntk klrkrliryk rlikqkleyr
721 qtvllnklie detqattnnt vekdnntler lelaqeltml qlqrknklss lvkkfednak
781 ihkyrriire gtemnieevd ssldvilqtl iannnknkga eqiitisnan sha
//
"


str(db$protein) # before
db$protein$sequence[db$protein$name == "Mbp1"] <- dbSanitizeSequence(mySeq)
str(db$protein) # after

# Analyze the expression ! Note how we specify an element of a vector (column)
# in a data frame in a list using a logical expression. And note how we assign
# the output (return value) of a function. As far as database coding goes this
# is pretty minimal - there is no error checking done at all. In particular: can
# we really guarantee that the name "Mbp1" is unique in the protein table? No!
# We never required it to be unique. This is a check we need to perform so
# frequently that we will encapsulate it in a function:

dbConfirmUnique

# try this:
dbConfirmUnique(c("TRUE", "FALSE"))
dbConfirmUnique(c(TRUE, FALSE))
dbConfirmUnique(c(TRUE, TRUE))
dbConfirmUnique(c(FALSE, FALSE))



#
#
# Here is the update to the sequence field of Res2 but using our
# confirmUnique() function


mySeq <- "
  1 maprssavhv avysgvevye cfikgvsvmr rrrdswlnat qilkvadfdk pqrtrvlerq
 61 vqigahekvq ggygkyqgtw vpfqrgvdla tkykvdgims pilsldideg kaiapkkkqt
121 kqkkpsvrgr rgrkpsslss stlhsvnekq pnssisptie ssmnkvnlpg aeeqvsatpl
181 paspnallsp ndntikpvee lgmleapldk yeeslldffl hpeegripsf lyspppdfqv
241 nsvidddght slhwacsmgh iemiklllra nadigvcnrl sqtplmrsvi ftnnydcqtf
301 gqvlellqst iyavdtngqs ifhhivqsts tpskvaaaky yldcilekli siqpfenvvr
361 lvnlqdsngd tslliaarng amdcvnslls ynanpsipnr qrrtaseyll eadkkphsll
421 qsnsnashsa fsfsgispai ispscsshaf vkaipsissk fsqlaeeyes qlrekeedli
481 ranrlkqdtl neisrtyqel tflqknnpty sqsmenlire aqetyqqlsk rlliwlearq
541 ifdlerslkp htslsisfps dflkkedgls lnndfkkpac nnvtnsdeye qlinkltslq
601 asrkkdtlyi rklyeelgid dtvnsyrrli amscginped lsleildave ealtrek
"

db$protein$sequence[dbConfirmUnique(db$protein$name == "Res2")] <- dbSanitizeSequence(mySeq)
str(db$protein)

# These expressions can get rather long, it's easier to read if we write:

select <- dbConfirmUnique(db$protein$name == "Res2")
value <- dbSanitizeSequence(mySeq)
db$protein$sequence[select] <- value

# ... and that's the code paradigm that we will adopt to update
# database fields (for now).

# Adding a record

# Adding a record is easy in principle, simply defining the values we get
# from the NCBI or EBI databases ... except for the ID field. That is a field
# we need to define internally, and once again, we'll use a small function
# to get this right.

dbAutoincrement

# After reading the code, predict the results
dbAutoincrement(1)      # Hm.
dbAutoincrement(1L)
dbAutoincrement(1:4)
dbAutoincrement(c(1, 2, 3, 4))
dbAutoincrement(c(1L, 2L, 3L, 4L))
dbAutoincrement(c(1, "4", 9))
dbAutoincrement(TRUE)

# Therefore, here is sample code to add one entry to the protein table.
#
mySeq <- "
  1 msgdktifka tysgvpvyec iinnvavmrr rsddwlnatq ilkvvgldkp qrtrvlerei
 61 qkgihekvqg gygkyqgtwi pldvaielae ryniqgllqp itsyvpsaad spppapkhti
121 stsnrskkii padpgalgrs rratsietes evigaapnnv segsmspsps dissssrtps
181 plpadrahpl hanhalagyn grdannhary adiildyfvt enttvpslli npppdfnpdm
241 sidddehtal hwacamgrir vvklllsaga difrvnsnqq talmratmfs nnydlrkfpe
301 lfellhrsil nidrndrtvf hhvvdlalsr gkphaaryym etminrlady gdqladilnf
361 qddegetplt maararskrl vrlllehgad pkirnkegkn aedyiieder frsspsrtgp
421 agielgadgl pvlptsslht seagqrtagr avtlmsnllh sladsydsei ntaekkltqa
481 hgllkqiqte iedsakvaea lhheaqgvde erkrvdslql alkhainkra rddlerrwse
541 gkqaikrarl qaglepgals tsnatnapat gdqkskddak sliealpagt nvktaiaelr
601 kqlsqvqank telvdkfvar areqgtgrtm aayrrliaag cggiapdevd avvgvlcell
661 qeshtgarag aggerddrar dvammlkgag aaalaanaga p
"

myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "UMAG_1122",
                    RefSeqID = "XP_011392621",
                    UniProtID = "A0A0D1DP35",
                    taxonomy.ID = as.integer(5270),
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(5270),
                    species = "Ustilago maydis",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

# Here is a bit of code template for the same ... it won't execute of course but
# you can copy it into your own code script file to modify when you need to add
# your own protein.

# ===== begin code template =====================================
mySeq <- "
RAW SEQUENCE FROM NCBI RECORD
"

myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "NAME",
                    RefSeqID = "REFSEQ",
                    UniProtID = "UNIPROT",
                    taxonomy.ID = as.integer(TAX_ID),
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(TAX_ID),
                    species = "BINOMIAL NAME",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

# ===== end code template ===========================================




# deleting a record

# This is simple code without error checking and of course we can make mistakes.
# Often we can just overwrite the offending field with correct data. But
# sometimes it will be easier (and more robust) to delete the erroneous entry
# and add the correct one. For example, if your code is in a script, and you
# realize the entry had an error, I would not "patch" the error in the script
# but delete the entry on the console command line, correct the script and
# execute the correct block. That way we fix the mistake at its source.
#
# Removing a row from a datframe is trivial: just overwrite the dataframe with
# a selection statement in which the unique selection of the offending row is
# inverted:

# create an erroneous entry
myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "nonesuch",
                    RefSeqID = "NP_000000",
                    UniProtID = "A123456",
                    taxonomy.ID = as.integer(999999),
                    sequence = dbSanitizeSequence("utter nonsense"),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

# Make a logical vector that identifies it
select <- dbConfirmUnique(db$protein$name == "nonesuch")
select     # the selection
!select    # its logical inversion

str(db)    # before
db$protein <- db$protein[ ! select, ]  # overwrite the table with a copy
                                       # without the selected record
str(db)    # after

# Note: if you delete records "by hand" you need to be careful that you do
# not remove keys that are used as foreign keys in another table - if there
# are such dependecies, you need to update the other table(s) too. "Real"
# database systems include such dependencies in the creation instructions
# of the table schema: "on delete cascade ..."

# ==============================================================================
#        PART THREE: Sequence Analysis
# ==============================================================================


if (!require(seqinr, quietly=TRUE)) {
    install.packages("seqinr")
    library(seqinr)
}

help(package = seqinr) # shows the available functions

# Let's try a simple function
?computePI

# This takes as input a vector of upper-case AA codes
# Let's retrieve the YFO sequence from our datamodel
# (assuming it is the last one that was added):

db$protein[nrow(db$protein), "sequence"]

# We can use the function strsplit() to split the string
# into single characters

s <- db$protein[nrow(db$protein), "sequence"]
s <- strsplit(s, "") # splitting on the empty spring
                     # splits into single characters
s <- unlist(s)       # strsplit() returns a list! Why?
                     # (But we don't need a list now...)

# Alternatively, seqinr provides
# the function s2c() to convert strings into
# character vectors (and c2s to convert them back).

s <- s2c(db$protein[nrow(db$protein), "sequence"])
s

computePI(s)  # isoelectric point
pmw(s)        # molecular weight
AAstat(s)     # This also plots the distribution of
# values along the sequence

# A true Labor of Love has gone into the
# compilation of the "aaindex" data:

?aaindex
data(aaindex)  # "attach" the dataset - i.e. make it accessible as an
               # R object

length(aaindex)

# Here are all the index descriptions
for (i in 1:length(aaindex)) {
    cat(paste(i, ": ", aaindex[[i]]$D, "\n", sep=""))
}


# Lets use one of the indices to calculate and plot amino-acid
# composition enrichment:
aaindex[[459]]

# === Sequence Composition Enrichment
#
# Let's construct an enrichment plot to compare one of the amino acid indices
# with the situation in our sequence.

refData <- aaindex[[459]]$I   # reference frequencies in %
names(refData) <- a(names(refData))  # change names to single-letter
# code using seqinr's "a()" function
refData


# tabulate our sequence of interest and normalize
obsData <- table(s)                # count occurrences
obsData = 100 * (obsData / sum(obsData))   # Normalize
obsData

len <- length(refData)

logRatio <- numeric() # create an empty vector

# loop over all elements of the reference, calculate log-ratios
# and store them in the vector
for (i in 1:len) {
    aa <- names(refData)[i] # get the name of that amino acid
    fObs <- obsData[aa]  # retrieve the frequency for that name
    fRef <- refData[aa]
    logRatio[aa] <- log(fObs / fRef) / log(2)  # remember log Ratio from
                                               # the lecture?
}

barplot(logRatio)

# Sort by frequency, descending
logRatio <- sort(logRatio, decreasing = TRUE)

barplot(logRatio)  # If you can't see all of the amino acid letters in the
                   # x-axis legend, make the plot wider by dragging the
                   # vertical pane-separator to the left



# label the y-axis
# (see text() for details)
label <- expression(paste(log[2],"( f(obs) / f(ref) )", sep = ""))

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        cex.names=0.9)



# color the bars by type.
# define colors
chargePlus  <- "#404580"
chargeMinus <- "#ab3853"
hydrophilic <- "#9986bf"
hydrophobic <- "#d5eeb1"
plain       <- "#f2f7f7"

# Assign the colors to the different amino acid names
barColors <- character(len)

for (i in 1:length(refData)) {
    AA <- names(logRatio[i])
    if (grepl("[HKR]",      AA)) {barColors[i] <- chargePlus }
    else if (grepl("[DE]",       AA)) {barColors[i] <- chargeMinus}
    else if (grepl("[NQST]",     AA)) {barColors[i] <- hydrophilic}
    else if (grepl("[FAMILYVW]", AA)) {barColors[i] <- hydrophobic}
    else                               barColors[i] <- plain
}

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        col = barColors,
        cex.names=0.9)


# draw a horizontal line at y = 0
abline(h=0)

# add a legend that indicates what the colours mean
legend (x = 1,
        y = -1,
        legend = c("charged (+)",
                   "charged (-)",
                   "hydrophilic",
                   "hydrophobic",
                   "plain"),
        bty = "n",
        fill = c(chargePlus,
                 chargeMinus,
                 hydrophilic,
                 hydrophobic,
                 plain)
)

# == TASK ==
# Interpret this plot. (Can you?)




# [END]
