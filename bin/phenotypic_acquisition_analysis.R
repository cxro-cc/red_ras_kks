# =========================================================================
# Copyright
# Jorge Arturo Arciniega González (arturo dot arciniegago at gmail dot com)
# Instituto de Ecologia, UNAM y Centro de Ciencias de la Complejidad, UNAM.
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# ========================================================================

#### Libraries ####
library(BoolNet)
library(BoolNetPerturb)      
library(parallel)
library(erer)  
library(tidyverse)
library(foreach)
library(doParallel)

# We define a function called data_net to load the data that will be # needed for the analysis. 
# necessary for the analysis. 1.- Logical functions, 2.

data_net <- function(network, labels) {
  bool_net <- loadNetwork(file = network)
  bool_labels <- as.data.frame(read.csv(file = labels, header = TRUE))
  net_data <- list(network = bool_net, labels = bool_labels)
  return(net_data)
}

# We define a function called attractors which evaluates the network in asynchronous mode
# and labels the attractors obtained. The parameters to give are:
# bool_net = the network to analyze, labels_net = the .csv file with the 
# labeling rules, geneON = over expression of a gene if any. 
# labeling rules, geneON = over expression of a gene if applicable and genesOFF = loss of function of a gene if applicable. 
# loss of function of a gene, if any.

attractors <- function(bool_net, labels_net, type = "asynchronous", 
                         states = 1, geneON = NULL, geneOFF = NULL) {
  
  net_attr <- getAttractors(network = bool_net, 
                            type = type,
                            startStates = states, 
                            genesON = c(geneON),
                            genesOFF = c(geneOFF),
                            method = "random")
  labelsObj <- labelAttractors(attr = net_attr, label.rules = labels_net)
  labeled_attr <- attractorsToLaTeX(net_attr)
  labeled_attr <- as.data.frame(labeled_attr)
  #colnames(labeled_attr) <- labelsObj
  
  if (ncol(labeled_attr) == 0) {
    labeled_attr <- c("complex attractor")
    # Attractors with stable states and chaotic attractors
  } else if (ncol(labeled_attr)!= length(labelsObj)) {
    print(paste0("complex attractor ", geneON, geneOFF ) )
    # Change in the name of the columns omitting the chaotic attractors
    colnames(labeled_attr) <- labelsObj[1:length(labeled_attr)]
    # Attractors with stable states
  } else if (ncol(labeled_attr)== length(labelsObj)){
    colnames(labeled_attr) <- labelsObj
  }
  return(labeled_attr)
  #labeled_attr <- t(labeled_attr) # convertir columnas a filas de la tabla resultante
  #return(labeled_attr)
}

# We define a function to extract all the labels associated with the attractors
pheno_extract <- function(table) {
  phenotypes <- c()
  for (i in 1:length(tabla)) {
    names <- colnames(tabla[[i]])
    phenotypes[[length(phenotypes) + 1]] <- names
  }
  fenotipos <- unlist(fenotipos)
  return(fenotipos)
}

# We define a function for counting the phenotypes characterized by the perturbed network
counts <- function(vector) {
  Hypertension <- length(str_view(vector, "Hyper"))
  Hypotension <- length(str_view(vector, "Hypo")) 
  No_label <- length(str_view(vector, "NoLabel"))
  phenotypes <- tibble(Hypertension = Hypertension, Hypotension = Hypotension, No_label = No_label)
  return(phenotypes)
}

# We define a parallel workflow for calculating the phenotypic acquisition frequency
# We initialize an empty vector where we will add the results of simulating 1000 times 
# the perturbed network and repeat this experiment 50 times to have good statistical power.

replicates_WT <- c()

for (i in 1:50) {
  # We perform 50 realizations to obtain statistical power
  # We initialize the parallelization
  num.cores <- detectCores() - 1
  cl <- makeCluster(num.cores, type="FORK", outfile="Log.txt")
  registerDoParallel(cl)
  
  simulationsWT <- foreach(iterators::icount(1000)) %dopar% {
    
    perturbedNet <- perturbNetwork(network = hiper$network,
                                   perturb = "functions", 
                                   method = "bitflip", 
                                   maxNumBits = 1)
    perturbedNetAttr <- attractors_2(bool_net = perturbedNet, 
                                     labels_net = hiper$labels, geneON = c("Renin", "CPN", "ANP"))
                                     }
  
  # We stop the cluster
  stopCluster(cl)
  
  # We extract the phenotypes
  rep <- pheno_extract(simulationsWT)
  
  # We calculate the phenotypic frecuency
  rep <- conteo(rep)
  
  replicates_WT[[length(replicates_WT) + 1]] <- rep
}

# we extract the contabilization of phenotypes
data_WT <- tibble()
data_WT <- add_column(data_WT, replicates_WT[[1]])
for (i in 1:length(replicates_WT)) {
  data_WT <- data_WT %>% add_row(replicates_WT[[i]])
}

# We save the data
write_csv(x = data_WT, file = "name_data.csv", col_names = TRUE)

#### Statistical analysis with ARTOOL ####

#The file "all_simulations.csv" contains all the repeitions and counts for the network experiment.
#You can find it in the meta directory
all_data <- as_tibble(read_csv(file = "all_simulations.csv", col_names = TRUE))

# We define the factors for this test
genotypes <- c("WT", "DW_ACE", "DW_ANP", "DW_IFNG", "UP_ACE2", "UP_IL10", "UP_NEP")
phenotypes <- c("Hypertension", "Hypotension")
all_data$Genotype <- factor(x = all_data$Genotype, levels = genotypes)

# Change the data to longer format and add phenotypes as factors
all_data_long <- all_data %>% pivot_longer(cols = c(Hypertension, Hypotension), names_to = "Phenotypes") %>% select(Genotype, Phenotypes, value)
all_data_long$Phenotypes <- factor(x = all_data_long$Phenotypes, levels = phenotypes)

## Applying the aligned rank transform to a factorial model ##
library(ARTool)
m.art <- art(value ~ Phenotypes * Genotype, data = all_data_long)
anova(m.art)

# Computing pairwise (post hoc) contrast between factors
contrast.m.art <- art.con(m = m.art, formula = "Phenotypes:Genotype", adjust = "bonferroni")
contr_results <- as_tibble(contrast.m.art)
contr_results <- as.data.frame(contrast.m.art)

# We calculate compact letter display for paired comparisons
library(rcompanion)
cld <- cldList(p.value ~ contrast, data=contr_results)
cld

# We summarise the mean and sd for the error bars and add the compact letter display for paired comparisons
dt <- all_data_long %>%
  group_by(Genotype, Phenotypes) %>%
  summarise(w = mean(value), sd = sd(value))
  
dt$cld <- c("e", "l", "a", "f", "b", "g", "a", "h", "c", "i", "d", "j", "a", "k")

# Finally we plot the results
ggplot(dt, aes(fill = Phenotypes,x = Genotype, y = w)) + 
  geom_col(width = 0.5, position = position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin = w-sd, ymax=w+sd), width = 0.5) +
  labs(x = "Genotypes", y = "Frequency of phenotypic acquisition") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3)) +
  geom_text(aes(label = cld), vjust = -0.5) +
  scale_fill_manual('Phenotypes', values=c('grey','black')) +
  theme_bw()
  
