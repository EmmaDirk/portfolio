# this script provides the analysis for the simulation study
# ---------------------------------------------------------------------

summary(lm(bp ~ HbA1C + age + as.factor(sex), data=dc))
confint(lm(bp ~ HbA1C + age + as.factor(sex), data=dc))
summary(lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc))
confint(lm(bp ~ HbA1C + bmi + age + as.factor(sex), data=dc))