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

#### LIBRERIAS ####
library(BoolNet)
library(BoolNetPerturb)
library(tidyverse)
library(erer)

#### FUNCIONES ####
# Definimos una función llamada data_net para cargar los datos que serán 
# necesarios para el análisis. 1.- Las funciones lógicas, 2.- Las estiquetas

data_net <- function(network, labels) {
  bool_net <- loadNetwork(file = network)
  bool_labels <- as.data.frame(read.csv(file = labels, header = TRUE))
  net_data <- list(network = bool_net, labels = bool_labels)
  return(net_data)
}

# Definimos una función llamada attractors que evalua la red en modo asíncrono
# y etiqueta los atractores obtenidos. Los parámetros que tienen que dar son:
# bool_net = la red a analizar, labels_net = el archivo .csv con las reglas de 
# etiquetado, geneON = sobre expresión de un gen si fuese el caso y genesOFF = 
# pérdida de función de algún gen si fuese el caso
attractors <- function(bool_net, labels_net, type = "asynchronous", 
                       states = 100, geneON = NULL, geneOFF = NULL) {
  
  net_attr <- getAttractors(network = bool_net, 
                            type = type,
                            startStates = states, 
                            genesON = c(geneON),
                            genesOFF = c(geneOFF))
  
  #Etiquetado de cada atractor obtenido
  labelsObj <- labelAttractors(attr = net_attr, label.rules = labels_net)
  #Conversión de atractores a dataframe
  labeled_attr <- attractorsToLaTeX(net_attr)
  labeled_attr <- as.data.frame(labeled_attr)
  
  #Para el etiquetado de los atractores obtenidos se pueden encontrar diferentes
  #atractores y por tanto diferentes reglas de etiquetado
    ## Atractores completamente caóticos/complejos
  if (ncol(labeled_attr) == 0) {
    labeled_attr <- c("complex attractor ")
    ## Atractores con estados estables y un atractor caótico
  } else if (ncol(labeled_attr)!= length(labelsObj)) {
    print(paste0("complex attractor ", geneON, geneOFF ) )
    # Cambio en el nombre de las columnas omitiendo a los atractores caóticos
    colnames(labeled_attr) <- labelsObj[1:length(labeled_attr)]
    ## Atractores con estados estables
  } else if (ncol(labeled_attr)== length(labelsObj)){
    colnames(labeled_attr) <- labelsObj
  }
  return(labeled_attr)
}

# Definimos una función llamada mutants que aplica la función attractors sobre c
# ada elemento de la red y el resultado lo guarda en una lista. 
# Si UP_REG = TRUE evalua las mutantes de sobre expresión
# Si DOWN_REG = TRUE evalua las mutantes de pérdida de función
# set_gene = nombre del gen que se quiere evaluar como constitutivo o pérdida de función entre comillas
# gene_value = valor del gen evaluado en set_value, este valor tiene que ser 0 o 1.
#   0 = pérdida de función
#   1 = ganancia de función
mutants <- function(bool_net, labels_net, UP_REG = FALSE, DOWN_REG = FALSE,
                    set_gene = NULL, gene_value = NULL) {
  bool_net <- fixGenes(network = bool_net, 
                       fixIndices = set_gene, 
                       values = gene_value)
  if(UP_REG == TRUE) {
    mut <- lapply(X = bool_net$genes, 
                  function (x) attractors(bool_net = bool_net, 
                                          geneON = x, 
                                          labels_net = labels_net))
    names(mut) <- paste("UP", bool_net$genes, sep = "_")
  } else if (DOWN_REG == TRUE) {
    mut <- lapply(X = bool_net$genes, 
                  function (x) attractors(bool_net = bool_net, 
                                          geneOFF = x, 
                                          labels_net = labels_net))
    names(mut) <- paste("DW", bool_net$genes, sep = "_")
  }
  return(mut)
}

# Escribir la lista que regresa la función mutants en un archivo .csv
write.list(z = objeto_a_escribir , file = "nombre_del_archivo.csv", row.names = TRUE)



