### Functions for mutation and labels ###

##----Libraries----##
library (BoolNet)
library (BoolNetPerturb)
library (tidyverse)

##----Mutation Function with labels (on//off)----##
# Network with fixed genes

### Overexpression
mut_on_fix <- function(net, labels, type = "asynchronous", startStates =  100, file = "on.csv"){
  #Ganes in network
  genes <- net$genes
  
  #Reading labels data
  labelsNet <- read.csv(file = labels, header = TRUE)
  
  for (node in genes){
    #Obtain each attractor with overexpression of every gene
    atr_on <- getAttractors(network = net,
                            type = type,
                            startStates = startStates,
                            genesON = node)
    #labeling each attractor with overexpresion
    labelsObj <- labelAttractors(attr = atr_on, label.rules = labelsNet, sep = "/")
    #Making a table of the attractor
    table <- attractorsToLaTeX(atr_on)
    #turn into data frame
    df_table <- as.data.frame(table)
    #replacing the column names with the labels
    colnames(df_table) <- labelsObj
    #writing every data frame in a csv file
    write.table(df_table, file = file, append = T, sep = ',')
  }
}


####Knockout

mut_off_fix <- function(net, labels, type = "asynchronous", startStates =  100, file = "off.csv"){
  #Ganes in network
  genes <- net$genes
  
  #Reading labels data
  labelsNet <- read.csv(file = labels, header = TRUE)
  
  for (node in genes){
    #Obtain each attractor with knockout of every gene
    atr_off <- getAttractors(network = net,
                            type = type,
                            startStates = startStates,
                            genesOFF = node)
    #labeling each attractor with knockout
    labelsObj <- labelAttractors(attr = atr_off, label.rules = labelsNet, sep = "/")
    #Making a table of the attractor
    table <- attractorsToLaTeX(atr_off)
    #turn into data frame
    df_table <- as.data.frame(table)
    #replacing the column names with the labels
    colnames(df_table) <- labelsObj
    #writing every data frame in a csv file
    write.table(df_table, file = file, append = T, sep = ',')
  }
}

##----Labeling function----##

label_attr <- function(attractor, labels){
  #attractor from net
  #labels as csv  file
  labelsNet <- read.csv(file = labels, header = TRUE)
  labelsObj <- labelAttractors(attr =  attractor, label.rules = labelsNet, sep = "/")
  table <- as.data.frame(attractorsToLaTeX(attractor))
  colnames(table) <- labelsObj
}

##----ON/OFF in one function ----##
mut_fix <- function(net, labels, type = "asynchronous", startStates =  100, file = "on.csv"){
  #name of every gene
  genes <- net$genes
  
  #obtain label data
  labelsNet <- read.csv(file = labels, header = TRUE)
  
  ##overexpression loop
  for (node in genes){
    #Obtain each attractor with overexpression of every gene
    atr_on <- getAttractors(network = net,
                            type = type,
                            startStates = startStates,
                            genesON = node)
    #labeling each attractor with overexpresion
    labelsObj <- labelAttractors(attr = atr_on, label.rules = labelsNet, sep = "/")
    #Making a table of the attractor
    table <- attractorsToLaTeX(atr_on)
    #turn into data frame
    df_table <- as.data.frame(table)
    #replacing the column names with the labels
    colnames(df_table) <- labelsObj
    #writing every data frame in a csv file
    write.table(df_table, file = file, append = T, sep = ',')
  }
  
  ##knockout loop
  for (node in genes){
    #Obtain each attractor with knockout of every gene
    atr_off <- getAttractors(network = net,
                             type = type,
                             startStates = startStates,
                             genesOFF = node)
    #labeling each attractor with knockout
    labelsObj <- labelAttractors(attr = atr_off, label.rules = labelsNet, sep = "/")
    #Making a table of the attractor
    table <- attractorsToLaTeX(atr_off)
    #turn into data frame
    df_table <- as.data.frame(table)
    #replacing the column names with the labels
    colnames(df_table) <- labelsObj
    #writing every data frame in a csv file
    write.table(df_table, file = file, append = T, sep = ',')
  }
}


