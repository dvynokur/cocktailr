# 1. Update documentation (from roxygen comments in R/ files)
devtools::document()

# 2. Load package into current R session (without reinstall)
devtools::load_all()

# 3. Run tests in tests/testthat/
devtools::test()

# 4. Run CRAN-like checks
devtools::check()

# 5. Install into your library
devtools::install()

# 6. Load installed package
library(cocktailr)


# Diagnostics
# Check for non-ASCII characters
tools::showNonASCIIfile("R/plot.R")

# Or check all R scripts
lapply(list.files("R", "\\.[rR]$", full.names = TRUE), tools::showNonASCIIfile)


# Git Reminders
git status     # check changes
git add .      # stage all changes
git commit -m "Message"  # commit

# Delete the broken installed copy
remove.packages("cocktailr", lib = .libPaths()[1])
unlink(file.path(.libPaths()[1], "cocktailr"), recursive = TRUE, force = TRUE)
unlink(file.path(.libPaths()[1], "00LOCK*"), recursive = TRUE, force = TRUE)

