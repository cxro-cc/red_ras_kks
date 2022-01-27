## Se modifica la red inicial de 32 nodos a una red de 29 nodos y se corre de manera sincróna y asincróna: 

# Cargar librerías
library(BoolNet)
library(BoolNetPerturb)

# Cargar la red 
net_ras <-  loadNetwork(file = "RAS292022.txt")
net_ras

# Obtener atractores en wild type.
# Colocar parametros. Con parámetro "sat.exhaustive" cuando es una red mayor a 29 genes:
# Colocar genes constitutivos
attr.ras.wt <- getAttractors(network = net_ras, 
                             type = "synchronous",
                             method = "sat.exhaustive",
                             genesON = c("Renin","CPN"))
                             
# Gráficar los atractores de la red en wild type
plotAttractors(attractorInfo = attr.ras.wt)

## Red de 29 genes se corre de manera asíncrona y se etiquetan los atractores:

## Se corre la red de 29 genes de manera asíncrona, primero se llama a la red que se guardo previamente, en un archivo de texto con terminación txt;
net <- loadNetwork("RAS292022.txt")
# Se imprime la red:
net
# Se fijan los genes constitutivos "CPN" &"Renin" y se imprime red:
net <- fixGenes(network = net, fixIndices = c("CPN","Renin"), values = c(1,1))
net
## Se corre la red de manera asíncrona, se inicia con un número de estados starStates de 100, 1000, 10,000, etc. :
atractor<- getAttractors(network = net, type = "asynchronous", startStates = 1000)
# Se imprimen los plots de los atractores obtenidos de manera asíncrona:
plotAttractors(atractor)

### Etiquetas (labels) de los atractores de manera asíncrona:

# Cargar la paquetería de BoolNetPerturb si aún no se cargo en un inicio y además, se carga la paquetería de devtools. 
# Nota: esta ultima solo se cargar una sola vez, no es necesario cargar cada vez que se abre el archivo:
library(BoolNetPerturb)
library(devtools)
# Se escriben las funciones lógicas, para cada fenotipo a estudiar en este caso tenemos los fenotipos "Hipertensión" e "Hipotensión". 
# Se guardan en un archivo de texto con la terminal csv y se llama al archivo con el siguiente código:
labelsNet <- read.csv(file = "labelKKS29.csv", header = TRUE)
# Imprime las funciones lógicas para verificar:
labelsNet
# Se coloca el siguiente código para etiquetar cada atractor, Nota: sep= / es por si existen los 2 fenotipos en un atractor y los pueda separar;
labelsObj <- labelAttractors(attr = atractor, label.rules = labelsNet, sep = "/")
# Nota: Puedes omitir el código table(labelsObj) dado que aparecen sin una numeración, y puedes pasar directamente a imprimir labelsObj:
table(labelsObj)
labelsObj
# Otra forma de visualizarlo los fenotipos por atractor, es uno por uno con el siguiente comando;
labelsObj[1]
labelsObj[2]
# Finalmente, cuando colocas un sapply visualizas los atractores en diferentes líneas, esto ayudara a resolver la situación,
# de que de manera asíncrona, coloca [1] hipo hiper, es decir, en una misma línea y sin separar con corchetes y números [1] y [2]
sapply(labelsObj,print)
#Para etiquetar la red de manera síncrona se siguen los mismos pasos anteriores, 
# pero se corre la red en un inicio de manera síncrona ejemplo:
net <- loadNetwork("RAS292022.txt")
net
net <- fixGenes(network = net, fixIndices = c("CPN","Renin"), values = c(1,1))
net
attractor<- getAttractors(network = net, type = "synchronous", method = "sat.exhaustive")
plotAttractors(attractor)
labels2 <- labelAttractors(attr = attractor, label.rules = labelsNet)
labels2
