# .utilities.R
# sundry utility code for the BCH441 course
#

init <- function() {
    source(".init.R")
}

# ===== PACKAGES

if (!require(seqinr, quietly=TRUE)) {
    install.packages("seqinr")
    library(seqinr)
}

# ====== Utilities

biCode <- function(s) {
    # make a 5 character code from a binomial name by concatening
    # the first three letter of the first word and the first
    # two letters of the second word.
    #
    # There are many ways to do this, here we assemble the two parts
    # in a loop, this way the function is vectorized and can
    # work on a large vector of names.
    b <- character()
    for (i in 1:length(s)) {
      b[i] <- sprintf("%s%s",
                      toupper(unlist(substr(s[i], 1, 3))),
                      toupper(unlist(substr(strsplit(s[i], "\\s+")[[1]][2],
                                            1, 2))))
    }
    return(b)
}

dotPlot2 <- function(A, B,        # sequences
                     MDM = BLOSUM62,
                     f,           # filter
                     palette,     # a function that returns color values
                     xlab = "",
                     ylab = ""
                     ) {
    # Purpose:
    #     Create a dotplot to measure sequence similarity between
    #     two amino acid sequences
    # Version:  1.0
    # Date:     2016-09
    # Author:   Boris Steipe
    #
    # Parameters:
    #     A, B: strings that contain no letters that are not found in MDM
    #     MDM: Mutation Data Matrix. Deafaults to BLOSUM62 which must be
    #          defined. Load package biostrings from BioConductor and load
    #          the matrix with data(BLOSUM62).
    #     f: filter matrix to weight an average around the neighborhood of
    #        an amino acid pair. Default to the identity matrix if missing.
    #        Average over a window of length f if length(f) is 1.
    #     palette: rainbow(), cm.colors(), or another function that returns
    #              a palette of color hexcodes. If missing, make our own
    #              palette.
    # Value:
    #     none. creates a dotplot.


    if (missing(f)) {
        f <- matrix(1) # default
    } else if (length(f) == 1) {
        if (! f %% 2) {stop("Sorry: f must be odd.")}
        w <- f
        f <- matrix(numeric(w * w), nrow = w)
        for (i in 1:w) { f[i, i] <- 1 }  # identity matrix
    }

    if (missing(palette)) {
        palette <- colorRampPalette(c("#000000",  #black
                                      "#111111",
                                      "#222222",
                                      "#332222",  # grey
                                      "#CC7755",
                                      "#EE8844",
                                      "#FF4400"), # red
                                      bias = 0.8)
    }
    A <- unlist(strsplit(A, ""))
    V <- unlist(strsplit(B, ""))
    lA <- length(A)
    lB <- length(B)

    m <- matrix(numeric(lA * lB), nrow = lA, ncol = lB)
    for (i in 1:lA) {
        for (j in 1:lB) {
            m[i, j] <- MDM[A[i], B[j]]
        }
    }
    m2 <- m
    wr <- floor((dim(f)[1] - 1) / 2)  # half-window size for rows
    wc <- floor((dim(f)[2] - 1) / 2)  # half-window size for columns

    for (i in (wr + 1):(lA - wr)) {
        for (j in (wc + 1):(lB - wc)) {
            # apply the filter to each value in m by weighting and summing
            # over its wr x wc neighborhood. Put the new value in m2
            m2[i, j] <- sum(f * m[(i-wr):(i+wr), (j-wc):(j+wc)])
        }
    }
    image(1:lA, 1:lB, m2,
          col = palette(24),
          ylim=c(lB,1), xlim=c(1,lA),
          xlab = xlab,
          ylab = ylab,
          axes = FALSE)
    box()

    # find good values for axis ticks and gridlines
    steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500,
               1000, 2000, 5000, 10000, 20000, 50000)
    gridStep <- sum(steps < max(lA, lB))

    # draw axes
    axis(1, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
    axis(2, at = c(1, seq(steps[gridStep - 3], lB, by=steps[gridStep - 3])))
    axis(3, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
    axis(4, at = c(1, seq(steps[gridStep - 2], lB, by=steps[gridStep - 3])))

    # draw grid with thin, transparent lines
    for (pos in seq(steps[gridStep - 2], lA, by = steps[gridStep - 2])) {
        abline(v=pos, col = "#FFFFFF44", lwd = 0.5)
    }
    for (pos in seq(steps[gridStep - 2], lB, by = steps[gridStep - 2])) {
        abline(h=pos, col = "#FFFFFF44", lwd = 0.5)
    }
}



# ====== DATABASE FUNCTIONS ============================================

dbInit <- function() {
    # Return an empty instance of the systems database
    # For the schema, see:
    # https://docs.google.com/presentation/d/1_nYWiwQc-9Z4anUAwOeVqWXgXIvM1Zl_AffMTM_No2Q

    db <- list()

    db$version <- "1.0"
    attributes(db$version)$code <- "ver"

    db$system <- data.frame(
        ID = character(),
        name = character(),
        notes = character(),
        stringsAsFactors = FALSE)
    attributes(db$system$ID)$code <- "sys"

    db$component <- data.frame(
        ID = character(),
        protein.ID = character(),
        system.ID = character(),
        status = character(),     # include / exclude / provisional
        notes = character(),
        stringsAsFactors = FALSE)
    attributes(db$component$ID)$code <- "cmp"

    db$protein <- data.frame(
        ID = character(),
        name = character(),
        RefSeqID = character(),
        UniProtID = character(),
        taxonomy.ID = numeric(),
        sequence = character(),
        stringsAsFactors = FALSE)
    attributes(db$protein$ID)$code <- "pro"

    db$taxonomy <- data.frame(
        ID = numeric(),
        species = character(),
        stringsAsFactors = FALSE)
    attributes(db$taxonomy$ID)$code <- "tax"

    db$systemAnnotation <- data.frame(
        ID = character(),
        system.ID = character(),
        feature.ID = character(),
        stringsAsFactors = FALSE)
    attributes(db$systemAnnotation$ID)$code <- "san"

    db$componentAnnotation <- data.frame(
        ID = character(),
        component.ID = character(),
        feature.ID = character(),
        stringsAsFactors = FALSE)
    attributes(db$componentAnnotation$ID)$code <- "can"

    db$proteinAnnotation <- data.frame(
        ID = character(),
        protein.ID = character(),
        feature.ID = character(),
        start = numeric(),
        end = numeric(),
        stringsAsFactors = FALSE)
    attributes(db$proteinAnnotation$ID)$code <- "pan"

    db$feature <- data.frame(
        ID = character(),
        name = character(),
        type.ID = character(),
        description = character(),
        sourceDB = character(),
        accession = character(),
        stringsAsFactors = FALSE)
    attributes(db$feature$ID)$code <- "fea"

    # type is one of:
    #     reference e.g. literature reference
    #     source (of information) e.g. literature reference
    #     3D domain e.g
    #     ptm       e.g. phosphorylation
    #     feature   e.g. low-complexity rtegion

    db$type <- data.frame(
        ID = character(),
        name = character(),
        description = character(),
        stringsAsFactors = FALSE)
    attributes(db$type$ID)$code <- "typ"

    return(db)
}


dbMerge <- function(my, ref = refDB) {
    # my and ref must be databases with the same schema
    if (my$version != ref$version) {
        stop("PANIC: Database version mismatch. Update schema before merging.")
    }

    tables <- names(ref)
    tables <- tables[tables != "version"]  # exclude "version" from the
                                           # list of tables to merge

    for (table in tables) {
        ref[[table]] <- tableUnion(ref[[table]], my[[table]])
    }

    return(ref)
}

tableUnion <- function(a, b) {
    # rbind() two tables, then keep only the rows with the first
    # occurrence of any duplicated ID
    a <- rbind(a, b)
    return(a[! duplicated(a$ID), ])
}


dbSanitizeSequence <- function(s, unambiguous = TRUE) {
    # Flatten any structure that s may have.
    # Remove all non-letters.
    # Convert to uppercase.
    # If "unambiguous" stop if any remaining letter
    #    matches an ambiguity code
    # else, return as-is.
    s <- paste(unlist(s), collapse="")
    s <- toupper(gsub("[^a-zA-Z]", "", s))
    if (unambiguous) {
        amb <- "([bjouxzBJOUXZ])"  # parentheses capture the match
        ambChar <- unlist(regmatches(s, regexec(amb, s)))[1]
        if (! is.na(ambChar)) {
            stop(paste("Input contains ambiguous letter(s): \"",
                       ambChar, "\".", sep=""))
        }
    }
    return(s)
}


dbConfirmUnique <- function(x) {
    # x is a vector of logicals.
    # returns x if x has exactly one TRUE element.
    # stop() otherwise.

    if (any(!is.logical(x))) {
        stop("PANIC: Input fails is.logical() test.")
    } else if (sum(x) == 0) {
        stop("PANIC: No TRUE value in input.")
    } else  if (sum(x) > 1) {
        stop("PANIC: More than one TRUE value in input.")
    } else {
        return(x)
    }
}


dbAutoincrement <- function(x, ns = "my") {
    # input: x: a character vector of keys structured as <namespace>_<table>_<integer>
    # output: a unique key of the same structure, in which the integer part
    # increments the largest integer in that namespace in x by one

    # get the table code either from its original attribute, or (since it gets
    # lost after merging or rbinding) from the first element in the table

    if (length(grep("code", names(attributes(x)))) == 1) {
        code <- attributes(x)$code
    } else {
        code <- strsplit(x[1], "_")[[1]][2]
    }

    iKeys <- grep(sprintf("^%s_", ns), x)  # get index of existing keys in ns
                                           # namespace in x
    if (length(iKeys) == 0) { # no key yet, return first key
        return(sprintf("%s_%s_1", ns, code))
    } else {
        nums <- as.integer(gsub("^.+_", "", x[iKeys])) # remove prefix
        if (any(is.na(nums))) { # non-digits in remaining string
            stop("PANIC: Input contains malformed key(s).")
        } else {
            return(sprintf("%s_%s_%d",
                           ns,
                           code,
                           as.integer(max(nums) + 1)))
        }
    }
}




# ====== SUPPORT FUNCTIONS =====================================================

waitTimer <- function(t, nIntervals = 50) {
    # pause and wait for t seconds and display a progress bar as
    # you are waiting
    t <- as.numeric(t)

    if (t < 0.1) {return(invisible())}

    increment <- t / nIntervals

    bar <- "----:----|"  # One module for the progress bar:
    bar <- rep(bar, ceiling(nIntervals / 10))  # repeat,
    bar <- unlist(strsplit(bar, "")) # split into single characters,
    bar <- bar[1:nIntervals]  # truncate,
    bar <- paste(bar, collapse="") # and collapse.

    cat(sprintf("\nWaiting: |%s\n         |", bar))
    for (i in 1:(nIntervals - 1)) {
        Sys.sleep(increment)
        cat("=")
    }
    Sys.sleep(increment)
    cat("|\n\n")

    return(invisible())
}


# ====== DATA ==================================================================


# 10 species of fungi for reference analysis.
# http://steipe.biochemistry.utoronto.ca/abc/index.php/Reference_species_for_fungi
REFspecies <- c("Aspergillus nidulans",
                "Bipolaris oryzae",
                "Coprinopsis cinerea",
                "Cryptococcus neoformans",
                "Neurospora crassa",
                "Puccinia graminis",
                "Saccharomyces cerevisiae",
                "Schizosaccharomyces pombe",
                "Ustilago maydis",
                "Wallemia mellicola"
)

# load the reference database refDB  (cf. create_refDB.R)
load("data/refDB.RData")

