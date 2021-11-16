# Cargar librerías
library(BoolNet)
library(BoolNetPerturb)

# Cargar la red 
net_ras <-  loadNetwork(file = "red_ras_kks.txt")
net_ras

# Obtener atractores en wild type.
# Colocar parametros. Con parametro "sat.exhaustive" por ser una red mayor a 29 genes
# Colocar genes constitutivos
attr.ras.wt <- getAttractors(network = net_ras, 
                             type = "synchronous",
                             method = "sat.exhaustive",
                             genesON = c("SERPING1","Renin","CPN"))
                             
# Gráficar los atractores de la red en wild type
plotAttractors(attractorInfo = attr.ras.wt)
