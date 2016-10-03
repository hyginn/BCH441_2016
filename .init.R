# .init.R
# Functions to initialize this Exercise session
# Boris Steipe
# ====================================================================

# if it doesn't exist yet, set up a profile:
if (! file.exists(".myProfile.R")) {
    # setup profile data
    cat("\nPlease enter the requested values correctly, no spaces, and\n")
    cat("press <enter>.\n")
    e <- readline("Please enter your UofT eMail address: ")
    n <- readline("Please enter your Student Number: ")

    conn <- file(".myProfile.R")
    writeLines(c(sprintf("myEMail <- \"%s\"", e),
                 sprintf("myStudentNumber <- %d", as.numeric(n))),
               conn)
    close(conn)
    rm(e, n, conn)
}

source(".myProfile.R")

file.edit("BCH441_2016.R")

# [End]
