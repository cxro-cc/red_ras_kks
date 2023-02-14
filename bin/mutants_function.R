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

# LIBRARIES ----
library(BoolNet)
library(BoolNetPerturb)
library(tidyverse)
library(erer)

# FUNCTIONS ----
## Reading data ----
#' Upload data to get the Boolnet Network
#' @param network file with logic rules of the net in .txt format
#' @param labels file with logic rules of the attractors in .csv format

data_net <- function(network, labels) {
  bool_net <- loadNetwork(file = network)
  bool_labels <- as.data.frame(read.csv(file = labels, header = TRUE))
  net_data <- list(network = bool_net, labels = bool_labels)
  return(net_data)
}


## Getting attractors ----
#' This function evaluates the network in the asynchronous mode for default, 
#' it also labels the attractors.
#' @param bool_net network from loadNetwork Boolnet function. If data proceeds from data_net function 
#' then object$network
#' @param labels_net labels in .csv format. If proceeds from data_net function 
#' then object$labels
#' @param geneON overexpression of gene if necessary
#' @param geneOFF knockout of gene if necessary
#' @param type transitions in network, default "asynchronous"
#' @param states how many states needs to evaluate

attractors <- function(bool_net, labels_net, type = "asynchronous", 
                       states = 10000, geneON = NULL, geneOFF = NULL) {
  
  net_attr <- getAttractors(network = bool_net, 
                            type = type,
                            startStates = states, 
                            genesON = c(geneON),
                            genesOFF = c(geneOFF))
  
  #Labeling attractor
  labelsObj <- labelAttractors(attr = net_attr, label.rules = labels_net)
  #Converting attractors to dataframe
  labeled_attr <- attractorsToLaTeX(net_attr)
  labeled_attr <- as.data.frame(labeled_attr)
  
  # We can find different attractors and therefore different rules for labeling
    ## Complex/chaotic attractors
  if (ncol(labeled_attr) == 0) {
    labeled_attr <- c("complex attractor ")
    ## Steady states attractors and chaotic attractor
  } else if (ncol(labeled_attr)!= length(labelsObj)) {
    print(paste0("complex attractor ", geneON, geneOFF ) )
    # Changing column names ignoring chaotic attractors
    colnames(labeled_attr) <- labelsObj[1:length(labeled_attr)]
    ## Steady states attractors
  } else if (ncol(labeled_attr)== length(labelsObj)){
    colnames(labeled_attr) <- labelsObj
  }
  return(labeled_attr)
}

## Getting mutants ----
#' It applies the attractors function on each network element and 
#' saves the results in a list
#' @param UP_REG if TRUE it evaluates overexpression mutants
#' @param DOWN_REG if TRUE it evaluates knockout mutants
#' @param set_gene name of gene that is constitutive or knockout
#' @param gene_value value of gene from set_gene, if constitutive 1, if knockout 0
#' 
mutants <- function(bool_net, labels_net, UP_REG = FALSE, DOWN_REG = FALSE,
                    set_gene = NULL, gene_value = NULL, states = 100) {
  bool_net <- fixGenes(network = bool_net, 
                       fixIndices = set_gene, 
                       values = gene_value)
  if(UP_REG == TRUE) {
    mut <- lapply(X = bool_net$genes, 
                  function (x) attractors(bool_net = bool_net, 
                                          geneON = x, 
                                          labels_net = labels_net,
                                          states = states))
    names(mut) <- paste("UP", bool_net$genes, sep = "_")
  } else if (DOWN_REG == TRUE) {
    mut <- lapply(X = bool_net$genes, 
                  function (x) attractors(bool_net = bool_net, 
                                          geneOFF = x, 
                                          labels_net = labels_net,
                                          states = states))
    names(mut) <- paste("DW", bool_net$genes, sep = "_")
  }
  return(mut)
}

