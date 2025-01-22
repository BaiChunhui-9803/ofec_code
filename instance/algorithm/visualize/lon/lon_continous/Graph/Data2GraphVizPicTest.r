library(igraph)
library(plyr)
library(rgl)     # 3D plots
library(foreach)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)




file= "E:/Diao_Yiya/code/OFEC/instance/algorithm/visualize/Result/Ackley4/data_filtered/data_Ackley4n2_p1.631000.txt"
output_folder="./"

#system(paste('mkdir', output_folder))

# PLOT3D
plotNet3D <- function(N, z, ewidth, vsize, asize, mylay) {
  mylayout3D <- cbind(mylay, z) # append z coordinate to the 2D layout
  rgl.open()      # open new window
  bg3d("white")   # set background to white
  rglplot(
    N,
    layout = mylayout3D,
    edge.width = ewidth,
    edge.arrow.size = asize,
    vertex.label = NA,
    vertex.size = vsize
  )
}

# Width of edges based on their weight
edgeWidth <- function(N, minw, maxw) {
  ewidth <- 1
  if (ecount(N) > 0) {
    ewidth <- (maxw * E(N)$Count) / max(E(N)$Count)
    ewidth = ifelse(ewidth > minw, ewidth, minw)
    ewidth = ifelse(ewidth < maxw, ewidth, maxw)
  }
  return (ewidth)
}

# Size of nodes based on their strength (weighted vertex degree)
nodeSizeStrength <- function(N, minsize, maxsize) {
  vsize <- maxsize
  if (ecount(N) > 0) {
    vsize =  2 * graph.strength(N, mode = "in")
    vsize =  ifelse(vsize > minsize, vsize, minsize)
    vsize =  ifelse(vsize < maxsize, vsize, maxsize)
  }
  return (vsize)
}


trace <-
  read.table(
    file,
    sep = "",
    header = TRUE,
    colClasses = c("integer", "integer", "character", "integer", "character")
  )


colnames(trace) <- c("run", "fit1", "node1", "fit2", "node2")


edgeFile = "E:/Diao_Yiya/paper/fitness_landscape/LON/TSPLIB/att532.tsp.s1_X10000_r1000_z1/att532.tsp.s1_X10000_r1000_z1.edges"

edgeTrace <-
  read.table(
    edgeFile,
    sep = "",
    header = TRUE,
    colClasses = c("integer", "integer",  "integer")
  )

colnames(edgeTrace) <- c("Start", "End", "Count")


testNodeRowNum  <- nrow(testNode)
for (r in 1:nrow(edgeTrace)){
  if ( edgeTrace[r,2]== -1) {
    edgeTrace[r,2]<- edgeTrace[r,1]
  }
}



tspEdge <- edgeTrace

nodeFile = "E:/Diao_Yiya/paper/fitness_landscape/LON/TSPLIB/att532.tsp.s1_X10000_r1000_z1/att532.tsp.s1_X10000_r1000_z1.nodes"

nodeTrace <-
  read.table(
    nodeFile,
    sep = "",
    header = TRUE,
    colClasses = c("integer", "integer")
  )
colnames(nodeTrace) <- c("node", "fitness")


tspNode <- matrix(0, nrow(nodeTrace), 3)
for (r in 1:nrow(tspNode)){
  tspNode[r,1] <- nodeTrace[r,1]
  tspNode[r,2] <- nodeTrace[r,2]
}

colnames(tspNode) <- c("Node", "Fitness", "Count")


for (r in 1:nrow(tspEdge)){
  nodeId <- tspEdge[r,1]
  tspNode[nodeId,3] <- tspNode[nodeId,3]+ tspEdge[r,3]
  nodeId <- tspEdge[r,2]
  tspNode[nodeId,3] <- tspNode[nodeId,3]+ tspEdge[r,3]
  
}


# rbind combine by rows
lnodes <-
  rbind(setNames(trace[c("node1", "fit1")], c("Node", "Fitness")),
        setNames(trace[c("node2", "fit2")], c("Node", "Fitness")))
ledges <- trace[c("node1", "node2")]



# ddply delete duplicand nodes or edges
nodes <- ddply(lnodes, .(Node, Fitness), nrow)
colnames(nodes) <- c("Node", "Fitness", "Count")
edges <- ddply(ledges, .(node1, node2), nrow)
colnames(edges) <- c("Start", "End", "Count")


nodes <- tspNode
edges <- tspEdge

# CREATE THE GRAPH
LON <- graph_from_data_frame(d = edges,
                             directed = T,
                             vertices = nodes)

# REMOVE SELF LOOPS
LON <- simplify(LON, remove.multiple = FALSE, remove.loops = TRUE)


best <- nodes[1,2]
for(r in 1:nrow(nodes)){
  if(best>nodes[r,2]){
    best <- nodes[r,2]
  }
}

# Compute the best fitness value (minimization problem)
#lineBest <- apply(nodes,2,min)
#best <- lineBest[2]
#best <- min(nodes$Fitness)



#lonlay <- layout_nicely(LON)
lonlay <- layout_with_fr(LON)

ew <- edgeWidth(LON, 1, 3)
vs <- nodeSizeStrength(LON, 2, 8)

# CMLON model
el <- as_edgelist(LON)
fits <- V(LON)$Fitness
names <- V(LON)$name
f1 <- fits[match(el[, 1], names)]
f2 <- fits[match(el[, 2], names)]
E(LON)[which(f2 < f1)]$type <-
  "improving"   # improving edges - Minimisation
E(LON)[which(f2 == f1)]$type <-
  "equal"      # equal fitness (neutral) edges
MLON <- LON
V(MLON)$Count <- 1
gnn <-
  subgraph.edges(LON, which(f2 == f1), delete.vertices = FALSE)     #TRUE Network with the same fitness value
nn_memb <- components(gnn, mode = "weak")$membership
CMLON <-
  contract.vertices(
    MLON,
    mapping = nn_memb,
    vertex.attr.comb = list(Fitness = "first", Count = "sum", "ignore")
  )
CMLON <- simplify(CMLON, edge.attr.comb = list(Count = "sum"))



#lonlay2 <- layout_nicely(CMLON)#layout_nicely(CMLON)#layout_with_kk

lonlay2 <- layout_with_fr(CMLON)#layout_nicely(CMLON)#layout_with_kk



ew2 <- edgeWidth(CMLON, 1, 3)
vs2 <- nodeSizeStrength(CMLON, 2, 8)


# Coloring Nodes
g_ncol <-  "Red"      # Global Optima Red "#d7191c"
l_ncol <-  "snow4"    # Local Optima Gray - LON and MLON model
ls_ncol <-
  "royalblue3"              # non-global sink (bottom of funnel)
myalpha = 255
palt <- c("pink", "lightskyblue")
palt = apply(col2rgb(palt), 2, function(x) {
  rgb(
    red = x[1],
    green = x[2],
    blue = x[3],
    alpha = myalpha,
    maxColorValue = 255
  )
})
pal <- c(g_ncol, ls_ncol,  palt)
pi_ecol <-  "gray25" # Opaque dark gray improvement edges

sinks_ids <- which(degree(CMLON, mode = "out") == 0)
sinks_fit <- vertex_attr(CMLON, "Fitness")[sinks_ids]
global_opt <- V(CMLON)[Fitness == best]
nsinks <- length(sinks_ids)
nglobals <- length(global_opt)
igs <- sinks_ids[sinks_fit == best]
ils <- sinks_ids[sinks_fit > best]
sing <-
  sum(strength(
    graph = CMLON,
    vids = igs,
    mode = "in",
    loops = FALSE
  ),
  na.rm = TRUE)
sinl <-
  sum(strength(
    graph = CMLON,
    vids = ils,
    mode = "in",
    loops = FALSE
  ),
  na.rm = TRUE)
gstrength <-  sing / (sing + sinl)
gbasins <-
  foreach(
    sink = igs,
    .packages = c("igraph", "foreach"),
    .inorder = TRUE
  ) %do% {
    subcomponent(graph = CMLON, v = sink, mode = "in")
  }
v <- unique(array(as.numeric(unlist(gbasins))))
gfunnel <-
  length(v) / vcount(CMLON)  # proportion of nodes that reach a global optimum
V(CMLON)$color <- pal[4]  # Basins of sub-optimal sinks
for (gb in gbasins) {
  V(CMLON)$color[V(CMLON) %in% gb] = pal[3]    # Color Optimal basin
}
V(CMLON)$color[V(CMLON) %in% sinks_ids] =   pal[2]  # Color of suboptimal sinks
V(CMLON)$color[V(CMLON)$Fitness == best] =  pal[1]  # Color of global sinks
E(CMLON)$color = pi_ecol # edges are all improving in CMLON


V(LON)$color <- pal[3]
V(LON)$color[V(LON)$Fitness == best] <- pal[1]
E(LON)$color = pi_ecol



################### SELECT THE MODEL AND THE TYPE OF VISUALIZATION

# 2D Visualization
png(file = paste0(output_folder, 'lon.png'))
plot(
  LON,
  vertex.label = NA,
  vertex.size = vs,
  edge.width = ew,
  edge.arrow.size = 0.2,
  layout = lonlay
)
dev.off()
png(file = paste0(output_folder, 'cmlon.png'))
plot(
  CMLON,
  vertex.label = NA,
  vertex.size = vs2,
  edge.width = ew2,
  edge.arrow.size = 0.2,
  layout = lonlay2
)
dev.off()
# 3D Visualization
zcoord <- V(LON)$Fitness
zcoord2 <- V(CMLON)$Fitness



######### LON #########

plotNet3D(
  N = LON,
  z = zcoord,
  ewidth = ew,
  vsize = vs,
  asize = 1,
  mylay = lonlay
)
# rotating to nice positioning
rgl.viewpoint(180, 90)
# size of image
par3d(windowRect = c(20, 30, 800, 800))
# saving
snapshot3d(paste0(output_folder, '3D_lon.png'), fmt = 'png')


######### CMLON #########
plotNet3D(
  N = CMLON,
  z = zcoord2,
  ewidth = ew2,
  vsize = vs2,
  asize = 1,
  mylay = lonlay2
)
# rotating to nice positioning
rgl.viewpoint(180, 90)
# size of image
par3d(windowRect = c(20, 30, 800, 800))
# saving
snapshot3d(paste0(output_folder, '3D_cmlon.png'), fmt = 'png')
# rotating to nice positioning
rgl.viewpoint(180, 90)

