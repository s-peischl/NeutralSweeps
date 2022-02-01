library(coalesceR)
library(igraph)
library(ape) #this gives the plotting functions for phylogenies
library(phyclust) #this contains ms (Hudsons program)

setwd("C:/Users/Antoine/OneDrive - Universitaet Bern/SurfSweepCodeBlocks/trees/")

#tree <- sim.tree(method = "hudson",sample = 20,current = 100,ancestral = 1000,time = 50)

#draw.tree(tree,labels = TRUE)
#tree


#tree <- make_tree(20, 3)
#plot(as.undirected(tree),layout=layout_as_tree,vertex.size=0,vertex.label=NA,edge.width=2,edge.color="blue")

#tree



################################################################ build tree from list of coalescence times
tau <- 80
coaltimes <- read.table(paste("FocalCoalTreeN20tau",as.character(tau),".dat",sep=""))

node_i <- coaltimes$V1 + 1
node_j <- coaltimes$V2 + 1
coalij <- coaltimes$V3
node_is_in <- 1:40


#hist(coalij,col="blue",xlim=c(0,tau))


treet <- make_empty_graph()
treet <- add_vertices(treet,40,weight=0)



local_i <- node_i[coalij==min(coalij)]
local_j <- node_j[coalij==min(coalij)]


local_i
local_j

weight0 <- min(coalij)
weight0
treet <- add_vertices(treet,1,weight=weight0)

new_node <- 40+1



while (length(node_i) >= 1)
{
  local_i <- node_i[coalij==min(coalij)]
  local_j <- node_j[coalij==min(coalij)]
  new_node
  while (length(local_i) >= 1)
  {
    treet <- add.edges(treet,c(node_is_in[local_i[1]],new_node))
    is_in <- node_is_in[local_i[1]]
    node_is_in[ (1:40)[ node_is_in == is_in ] ] <- new_node
    treet <- add.edges(treet,c(node_is_in[local_j[1]],new_node))
    is_in <- node_is_in[local_j[1]]
    node_is_in[ (1:40)[ node_is_in == is_in ] ] <- new_node
    peche_a_ete_bonne <- 1
    
    while (peche_a_ete_bonne == 1)
    {
      remove_index <- vector()
      peche_a_ete_bonne <- 0
      for (k in 1:length(local_i))
      {
        if (is.element(local_i[k], (1:40)[ node_is_in == node_is_in[ local_i[1] ] ] ) )
        {
          if (node_is_in[local_j[k]] != new_node)
          {
            treet <- add.edges(treet,c(node_is_in[local_j[k]],new_node))
            is_in <- node_is_in[local_j[k]]
            node_is_in[ (1:40)[ node_is_in == is_in ] ] <- new_node
            peche_a_ete_bonne <- 1
          }
          remove_index <- c(remove_index,k)
        }
        
        if (is.element(local_j[k], (1:40)[ node_is_in == node_is_in[ local_i[1] ] ] ) )
        {
          if (node_is_in[local_i[k]] != new_node)
          {
            treet <- add.edges(treet,c(node_is_in[local_i[k]],new_node))
            is_in <- node_is_in[local_i[k]]
            node_is_in[ (1:40)[ node_is_in == is_in ] ] <- new_node
            peche_a_ete_bonne <- 1
          }
          remove_index <- c(remove_index,k)
        }
      }
    }
    
    treet <- add_vertices(treet,1,weight=weight0)
    new_node <- new_node + 1
    
    local_i <- local_i[-remove_index]
    local_j <- local_j[-remove_index]
  }
  
  node_i <- node_i[-(1:length(node_i))[coalij==min(coalij)]]
  node_j <- node_j[-(1:length(node_j))[coalij==min(coalij)]]
  coalij <- coalij[-(1:length(coalij))[coalij==min(coalij)]]
  
  if (length(coalij)!=0)
  {
    weight0 <- min(coalij)
  }
  
  treet <- set_vertex_attr(treet,"weight",new_node , weight0)
  
}

num_nodes <- new_node - 1

treet <- delete_vertices(treet,new_node)
###############################################################

#plot(as.undirected(treet),layout=layout.reingold.tilford,vertex.size=10,edge.width=2,edge.color="blue")






weights <- get.vertex.attribute(treet)$weight   
weights

################################################### check that it works by checking coalescence times
nodes_above <- list()

for(i in 1:40)
{
  above <- i
  nodes_above[[i]] <- i
  
  while (above < num_nodes)
  {
    neighbors <- ego(treet,order=1,above)[[1]]
    above <- neighbors[weights[neighbors]>weights[above]]
    nodes_above[[i]] <- c(nodes_above[[i]],above)
  }
}

coal <- vector()
k <- 1
for (i in 1:39)
{
  
  for (j in (i+1):40)
  {
    coal[k] <- min(weights[ nodes_above[[i]][ is.element(nodes_above[[i]],nodes_above[[j]]) ] ])
    k <- k+1
  }
}
## the min of the weight of the intersection between the set of vertices above i and the set of vertices above j 

coalij <- coaltimes$V3
length(coalij)
sum(abs(coal)-abs(coalij))
head(coal)
head(coalij)
#####################################################33




treet2 <- unroot(rtree(n = num_nodes))  ### object of class phylo



treet2$edge <- ends(treet,E(treet))[,c(2,1)]   ##### make sure the first node in the edge is above the second one
treet2$tip.label <- 1:40
treet2$Nnode <- num_nodes - 40
treet2$edge.length <- weights[treet2$edge[,1]]-weights[treet2$edge[,2]]    




################################## make sure that there is no branches crossing each other, (tips are plotted by plot.phylo in the order of appearance in the edge list)
num_edges<- length(treet2$edge[,1])
below1 <- num_nodes
below <- treet2$edge[(1:num_edges)[treet2$edge[,1] == below1] ,2 ]
below
while(sum(weights[below]) != 0)
{
  below1 <- below
  below <- vector()
  for (i in 1:length(below1))
  {
    if (weights[below1[i]] == 0)
    {
      below <- c(below,below1[i])
    }
    if (weights[below1[i]] != 0)
    {
      below <- c(below,treet2$edge[(1:num_edges)[treet2$edge[,1] == below1[i]] ,2 ])
    } 
  }
}
below
##################################





#####################################   internal nodes with lower labels at the top
reverseNodes <- c(1:40,num_nodes:41)  


weights <- weights[reverseNodes]   #### update weights accordingly

for (i in 1:num_edges)   ##### update edges accordingly
{
  treet2$edge[i,] <- c(reverseNodes[treet2$edge[i,1]],reverseNodes[treet2$edge[i,2]])
}
######################################



##########################
order <- below
treet2$edge[treet2$edge[,2] <= 40,]
disorderedEdges <- treet2$edge[treet2$edge[,2] <= 40,] 
orderedEdges <- disorderedEdges
for (i in 1:length(orderedEdges[,1]))
{
  orderedEdges[i,] <- c(disorderedEdges[  disorderedEdges[,2] == order[i]  ,1],disorderedEdges[ disorderedEdges[,2] == order[i]  ,2])
}
orderedEdges

treet2$edge[,1] <- c(orderedEdges[,1], treet2$edge[treet2$edge[,2] > 40,1] ) 
treet2$edge[,2] <- c(orderedEdges[,2], treet2$edge[treet2$edge[,2] > 40,2] ) 

treet2$edge.length <- weights[treet2$edge[,1]]-weights[treet2$edge[,2]]   #### update weights


treet2$root.edge <- tau - max(weights)

plot.phylo(treet2,type = "cladogram",edge.width = 3,direction="right",use.edge.length = TRUE,show.tip.label=FALSE,root.edge = TRUE,edge.color = "red")  ### clado or phylo
#nodelabels()

is.rooted(treet2)


