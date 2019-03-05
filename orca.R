#################################################
## Population dynamics of orca (= killer whales; Orcinus orca)
## from Brault & Caswell 1993 Ecology
## DOI: 10.2307/1940073
#################################################
## 0. Preliminaries
install.packages(c("popbio", "ggplot2", "dplyr", "tibble", "tidyr"))

library(popbio)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

# define function to generate random population composition
rdirichlet <- function(a) {
  y <- rgamma(length(a), a, 1)
  return(y / sum(y))
}

#################################################
## 1. Enter the orca matrix data
#################################################
orca_U <- matrix(c(0, 0, 0, 0,
                   0.9775, 0.9111, 0, 0,
                   0, 0.0736, 0.9534, 0,
                   0, 0, 0.0452, 0.9804),
                 nrow = 4, byrow = TRUE)
orca_F <- matrix(c(0, 0.0043, 0.1132, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0),
                 nrow = 4, byrow = TRUE)
orca_C <- matrix(c(0, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0),
                 nrow = 4, byrow = TRUE)

orca_A <- orca_U + orca_F + orca_C
orca_A  # this is the projection matrix

#################################################
## 2. Project the population forward in time
#################################################
# specify a starting population with 10 small adults
pop0 <- matrix(c(0,
                 10,
                 0,
                 0),
               nrow = 4, byrow = TRUE)

# recursive matrix multiplication
(pop1 <- orca_A %*% pop0)  # population composition change from t0->t1
sum(pop1)/sum(pop0)  # population growth from t0->t1
(pop2 <- orca_A %*% pop1)
sum(pop2)/sum(pop1)
(pop3 <- orca_A %*% pop2)
sum(pop3)/sum(pop2)
(pop4 <- orca_A %*% pop3)
sum(pop4)/sum(pop3)

# automate population projection
(orca_short <- pop.projection(orca_A, pop0, iterations = 5))  # to t4
(orca_long <- pop.projection(orca_A, pop0, iterations = 51))  # to t50

# make plots
plot(orca_long$pop.sizes, type = "l", ylab = "N", xlab = "time")  # population size
plot(orca_long$pop.changes, type = "l", ylab = "N", xlab = "time")  # population growth rate

orca_project_df <- data.frame(t(orca_long$stage.vectors)) %>%  # arrange data to plot stage structure
  rownames_to_column(var="time") %>%
  rename(calf = X1, sm_adult = X2, la_adult = X3, pr_adult = X4) %>%
  mutate(time = as.numeric(time),
         calf = calf/rowSums(.[2:5]), sm_adult = sm_adult/rowSums(.[2:5]),
         la_adult = la_adult/rowSums(.[2:5]), pr_adult = pr_adult/rowSums(.[2:5])) %>%
  gather("stage","proportion",2:5) %>%
  mutate(stage = ordered(stage, levels = c("calf", "sm_adult", "la_adult", "pr_adult")))
ggplot(orca_project_df, aes(x = time, y = proportion, group = stage, color = stage))+
  geom_line(size=3)

# equilibrium values
lambda(orca_A)
stable.stage(orca_A)

#################################################
## 3. Unstructured orca population dynamics
#################################################
# Population at stable stage structure (w)
orca_w <- stable.stage(orca_A)
orca_w_death <- 1-weighted.mean(colSums(orca_U), w = orca_w)
orca_w_birth <- weighted.mean(colSums(orca_F), w = orca_w)
(orca_w_lambda <- 1+(orca_w_birth-orca_w_death))
lambda(orca_A)

# Population not at stable stage structure
orca_notw <- rdirichlet(c(1, 1, 1, 1))  # generate a random population composition
orca_notw_death <- 1-weighted.mean(colSums(orca_U), w = orca_notw)
orca_notw_birth <- weighted.mean(colSums(orca_F), w = orca_notw)
(orca_notw_lambda <- 1+(orca_notw_birth-orca_notw_death))
lambda(orca_A)

