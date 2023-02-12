########################
#    RAS attractors    #
########################

# Getting the attractors in wild type and labeling
# Getting each attractor after the mutation of each node (KO and Gain function)


# Libraries ----
library(BoolNet)
library(BoolNetPerturb)
library(tidyverse)
library(erer)

# Getting Attractors with labels ----

# Functions
source("mutants_function.R", local = TRUE)

## Getting network with attractor labels ----
netRAS <- data_net("meta/rules_ras.txt", "meta/labels_ras.csv")

## Getting attractors and labeling ----
rasAttr <-  attractors(netRAS$network, netRAS$labels, 
                         states = 10000,
                         geneON = c("Renin", "CPN"))

# Robustness, mutation of each node ----

## Getting mutants (KO and gain function) for each node ----
mutRasKO <- mutants(netRAS$network, netRAS$labels, DOWN_REG = TRUE,
                      set_gene = c("Renin", "CPN"), gene_value = c(1, 1), states = 10000)

mutRasG <- mutants(netRAS$network, netRAS$labels, UP_REG = TRUE,
                     set_gene = c("Renin", "CPN"), gene_value = c(1, 1), states = 10000)

## Saving results for each mutant ----
write.list(z = mutRasKO, file = "off_rasNet.csv", row.names = TRUE)
write.list(z = mutRasG, file = "on_rasNet.csv", row.names = TRUE)

