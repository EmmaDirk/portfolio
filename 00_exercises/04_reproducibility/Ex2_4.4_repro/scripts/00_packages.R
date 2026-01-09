# this script loads the required packages for this simulation study
# and installs them if not already installed
# ---------------------------------------------------------------------

# required packages
pkgs <- c("Hmisc", "mice", "tidyverse", "here")

# install packages if not already installed
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install) > 0) {
  install.packages(to_install)
}

# load packages
lapply(pkgs, library, character.only = TRUE)

