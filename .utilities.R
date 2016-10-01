# .utilities.R
# sundry utility code for the BCH441 course
#

init <- function() {
    source(".init.R")
}


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


dbSanitizeSequence <- function(s) {
    # drop all characters that are not one of the codes for an amino acid
    # from the input s, convert to uppercase and return
    return(toupper(gsub("[^a-zA-Z]", "", s)))
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


dbAutoincrement <- function(x) {
    # input: x: a vector of integers
    # output: an integer one larger than the largest integer in x
    if (any(!is.integer(x))) {
        stop("PANIC: Input contains non-integer(s).")
    } else {
        return(as.integer(max(x) + 1))
    }
}


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


