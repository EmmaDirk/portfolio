# this script loads the data for the simulation study
# ---------------------------------------------------------------------

# set the here package root to the current directory, in case it is not already set
# file.create("../.here")

# The data can be dowloaded in xpt form from https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2015

# read data:
d1 <- sasxport.get(here("data","DEMO_I.xpt"))
d2 <- sasxport.get(here("data","BPX_I.xpt"))
d3 <- sasxport.get(here("data","BMX_I.xpt"))
d4 <- sasxport.get(here("data","GHB_I.xpt"))
d5 <- sasxport.get(here("data","TCHOL_I.xpt"))

# subset and merge data:
d1.t <- subset(d1,select=c("seqn","riagendr","ridageyr"))
d2.t <- subset(d2,select=c("seqn","bpxsy1"))
d3.t <- subset(d3,select=c("seqn","bmxbmi"))
d4.t <- subset(d4,select=c("seqn","lbxgh"))
d5.t <- subset(d5,select=c("seqn","lbdtcsi"))

d <- merge(d1.t,d2.t)
d <- merge(d,d3.t)
d <- merge(d,d4.t)
d <- merge(d,d5.t)