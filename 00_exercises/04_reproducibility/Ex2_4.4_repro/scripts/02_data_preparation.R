# this script prepares the data for the simulation study
# ---------------------------------------------------------------------

# rename variables:
# RIAGENDR - Gender
d$age <- d$ridageyr

# remove individuals <18 years
d$age[d$age<18] <- NA

# RIDAGEYR - Age in years at screening
d$sex <- d$riagendr

# BPXSY1 - Systolic: Blood pres (1st rdg) mm Hg
d$bp <- d$bpxsy1

# BMXBMI - Body Mass Index (kg/m**2)
d$bmi <- d$bmxbmi

# LBDTCSI - Total Cholesterol (mmol/L)
d$HbA1C <- d$lbxgh

# LBXGH - Glycohemoglobin (%)
d$chol <- d$lbdtcsi

# select complete cases:
dc <- cc(subset(d,select=c("age","sex","bmi","HbA1C","bp")))
