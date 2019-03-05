#################################################
## Population dynamics of loggerhead sea turtles (= Caretta caretta)
## from Crowder et al. 1994 Ecological Applications
## DOI: 10.2307/1941948
#################################################
## 0. Preliminaries
#################################################
# install.packages("popbio")

library(popbio)

#################################################
## 1. Enter the turtle matrix data
#################################################
turtle_U <- matrix(c(0, 0, 0, 0, 0,
                     0.675, 0.703, 0, 0, 0,
                     0, 0.047, 0.657, 0, 0,
                     0, 0, 0.019, 0.682, 0,
                     0, 0, 0, 0.061, 0.8091),
                   ncol = 5,byrow = TRUE)
turtle_F <- matrix(c(0, 0, 0, 4.665, 61.896,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0),
                   ncol = 5,byrow = TRUE)
turtle_A <- turtle_U+turtle_F

#################################################
## 2. Conservation interventions
#################################################
## 2.1 Baseline lambda
lambda(turtle_A)

## 2.2 Improving hatchling survival
turtle_A_hatch<-turtle_A
turtle_A_hatch[2, 1] <- 1  # set hatchling survival to 100%

lambda(turtle_A_hatch)

## 2.3 Turtle exclusion devices
ctrl_surv <- colSums(turtle_U)  # calculate stage-specific survival

# "Although all of the US loggerhead populations have been subject to trawling pressure,
# a well-studied loggerhead population in Australia, that is relatively isolated from trawling
# activity, may provide reasonable survivorship estimates for a relatively unexploited population."
aust_surv <- c(NA, NA, 0.830, 0.885, 0.910)

# "To predict the population response for US turtles with TEDs in place, we took the difference
# between Australian and US survivorships for large juveniles, subadults, and adults, then multiplied
# those differences by 0.67."
ted_surv <- ctrl_surv+(aust_surv-ctrl_surv)*0.67
ted_surv[1:2] <- ctrl_surv[1:2]  # keep hatchling and small juvenile survival the same

# The result for each stage was then added to the current survivorship
# estimate. This assumes that sources of mortality are additive and that upon reducing trawling-related
# mortality, mortality from other causes would remain unchanged."
turtle_U_ted <- (turtle_U*matrix(rep(ctrl_surv^-1, times = 5), nrow = 5, byrow = TRUE))*
  matrix(rep(ted_surv, times = 5), nrow = 5, byrow = TRUE)

# recalculate the new A matrix for the TED scenario 
turtle_A_ted <- turtle_U_ted + turtle_F

lambda(turtle_A_ted)

#################################################
## 3. Sensitivity and elasticity analysis
#################################################
(turtle_sens <- sensitivity(turtle_A, zero = TRUE))
image2(turtle_sens)

(turtle_elas <- elasticity(turtle_A))
image2(turtle_elas)

#################################################
## 4. COMADRE Database
#################################################
download.file("http://www.compadre-db.org/Data/ComadreDownload", destfile = "comadre.RData")
load("comadre.RData")
View(comadre$metadata[comadre$metadata$SpeciesAccepted == "Caretta_caretta", ])
comadre$mat[comadre$metadata$SpeciesAuthor == "Caretta_caretta_2"]
