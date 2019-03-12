#################################################-
## Exercise key: Teasel ----
## W.K. Petry
#################################################-
## Preliminaries ----
#################################################-
library(popbio)
library(ggplot2)
library(dplyr)
library(tidyr)

#################################################-
## Question A ----
#################################################-
# survival & growth
Umat <- matrix(c(0,     0,     0,     0,     0,     0, 
                 0.966, 0,     0,     0,     0,     0, 
                 0.013, 0.010, 0.125, 0,     0,     0, 
                 0.007, 0,     0.125, 0.238, 0,     0, 
                 0.008, 0,     0.036, 0.245, 0.167, 0, 
                 0,     0,     0,     0.023, 0.750, 0),
               nrow = 6, byrow = TRUE)
colSums(Umat)

Fmat <- matrix(c(0, 0, 0, 0, 0, 322.38, 
                 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 3.448, 
                 0, 0, 0, 0, 0, 30.170, 
                 0, 0, 0, 0, 0, 0.862, 
                 0, 0, 0, 0, 0, 0),
               nrow = 6, byrow = TRUE)

Cmat <- matrix(c(0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0, 0),
               nrow = 6, byrow = TRUE)

teasel <- Umat + Fmat + Cmat
dimnames(teasel) <- list(c("DS1", "DS2", "RosS", "RosM", "RosL", "Flo"),
                         c("DS1", "DS2", "RosS", "RosM", "RosL", "Flo"))

image2(teasel)

#################################################-
## Question B ----
#################################################-
pop0 <- matrix(c(50,
                 50,
                 0,
                 0,
                 0,
                 0),
               nrow = 6, byrow = TRUE)
pop15 <- pop.projection(A = teasel, n = pop0)$stage.vectors[, "15"]
sum(pop15)

proj_df <- as.data.frame(t(pop.projection(A = teasel,
                                          n = pop0,
                                          iterations = 16)$stage.vectors)) %>%
  mutate(time = 0:15) %>%
  gather("stage", "number", -time)

ggplot(proj_df, aes(x = time, y = number, color = stage))+
  geom_line(size = 2)+
  # scale_y_log10()+
  scale_color_discrete()+
  theme_bw()

#################################################-
## Question C ----
#################################################-
lambda(teasel)

#################################################-
## Question D ----
#################################################-
stable.stage(teasel) * 100

# notice that the starting population is close to the equalibrium
# stage structure
pop15/sum(pop15) * 100
