#Título: Funciones para el análsis de la red RAS

##################################################

#Cargar librerías
#library(BoolNet)
#library(tidyverse)

### FUNCIÓN PARA OBTENER MUTANTES DE SOBREEXPRESIÓN
mutantes_on <- function(rules, constitutives = NULL, type = "asynchronous", startStates = 100){
  #rules, archivo txt que contiene las reglas lógicas de la  red
  #constitutives puede ser un vector o un str 

  #Leer el archivo de reglas como una red de Boolnet
  net <- loadNetwork(rules)
  #Selecciona los nombres de los genes implicados en la red
  genes <- net$genes
  #print(genes)
  #Se abre un pdf en el cual se guardan los atractores creados por la mutación de cada nodo
  pdf("on.pdf")
  
  #Si la función contiene nodos constitutivos
  if(!is.null(constitutives)){
    #Se eliminan los nodos constitutivos de la lista de genes
    w_const <- genes[! genes %in% constitutives]
    #Cada uno de estos nodos (constitutivos y no constitutivos) es reemplazado en la función de mutación 'genesON'
    for(nodo in w_const){
      atr_on <- getAttractors(network = net, type = type, startStates = startStates, 
                              genesON = c(nodo,constitutives))
      plotAttractors(atr_on, title = paste0("overexpression ", nodo))
    }
  }
  #Si la función no contiene nodos constitutivos
  else{
    for(nodo in genes){
      #Cada uno de los nodos es reemplazado en la función de mutación 'genesON'
      atr_on <- getAttractors(network = net, type = type, startStates = startStates, 
                              genesON = nodo)
      plotAttractors(atr_on, title = paste0("overexpression ", nodo))
    }
  }
  #Se cierra el archivo PDF
  dev.off()
}


### FUNCIÓN PARA OBTENER MUTANTES DE INHIBICIÓN
mutantes_off <- function(rules, constitutives = NULL, type = "asynchronous", startStates = 100 ){
  #rules, archivo txt que contiene las reglas lógicas de la  red
  #constitutives puede ser un vector o un str

  #Leer el archivo de reglas como una red de Boolnet
  net <- loadNetwork(rules)
  #Selecciona los nombres de los genes implicados en la red
  genes <- net$genes
  #print(genes)
  #Se abre un pdf en el cual se guardan los atractores creados por la mutación de cada nodo
  pdf("off.pdf")
  
  #Si la función contiene nodos constitutivos
  if(!is.null(constitutives)){
    #Se eliminan los genes constitutivos del vector de nodos
    w_const <- genes[! genes %in% constitutives]
    
    #Cada uno de los nodos no constitutivos se mutan con 'genesOFF' y se mantienen los constitutivos con 'genesON'
    for(nodo in w_const){
      atr_on <- getAttractors(network = net, type = type, startStates = startStates, 
                              genesON = constitutives, genesOFF = nodo)
      plotAttractors(atr_on, title = paste0("knockout ", nodo))
    }
  }
  #En caso de que no hayan constitutivos solo se muta cada nodo  con 'genesOFF'
  else{
    for(nodo in genes){
      atr_on <- getAttractors(network = net, type = type, startStates = startStates, 
                              genesOFF = nodo)
      plotAttractors(atr_on, title = paste0("knock out ", nodo))
    }
  }
  dev.off()
}


#################################

### FUNCIONES PARA REALIZAR ANÁLISIS DE DERRIDA
derridaCurve <- function(network, numberOfStates)
{
  points <- sapply(1:numberOfStates,function(i)
  {
    # sample random initial states
    state1_init <- round(runif(n=length(network$genes)))
    state2_init <- round(runif(n=length(network$genes)))
    # calculate state transitions
    state1_next <- stateTransition(network,state1_init)
    state2_next <- stateTransition(network,state2_init)
    # return the Hamming distances (x and y value)
    return(c(sum(state1_init != state2_init),
             sum(state1_next != state2_next)))
  })
  # calculate the means of the y values for all possible x values,
  # and normalize by the number of genes
  curve <- sapply(1:length(network$genes), function(i)
  {
    indices <- which(points[1,] == i)
    return(c(i,mean(points[2,indices])))
  })/length(network$genes)
  # return a matrix with the x values in the first row
  # and the y values in the second row
  curve <- cbind(c(0,0),curve)
  return(curve)
}

compareToRandom <- function(network)
{
  # generate a random network with the same function in-degrees
  60
  randomNetwork <- generateRandomNKNetwork(n=length(network$genes),
                                           k=sapply(network$interactions,
                                                    function(int)length(int$input)))
  # calculate the two Derrida curves
  points1 <- derridaCurve(network, 10000)
  points2 <- derridaCurve(randomNetwork, 10000)
  # plot the diagonal
  plot(c(0,1), c(0,1), type="l", xlim=c(0,1), ylim=c(0,1), xlab=expression(h[t]),
       ylab=expression(h[t+1]), main = "Test de Derrida")
  # plot the curves
  lines(points1[1,], points1[2,], col="red")
  lines(points2[1,], points2[2,], col="green")
}



