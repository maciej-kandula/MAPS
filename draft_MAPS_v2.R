#install.packages(c("igraph", "NetIndices", "mcclust", "mixer", "blockmodels"))
#library(NetIndices)
#library(mcclust)
#library(mixer) # for SBM community detection

### load essential packages ###
library(igraph) # for graph preprocessing, vi, etc.
library(blockmodels) # for SBM methods (instead of Wmixnet; same author)

############################
### --- START CONFIG --- ###
############################

### General params ###
netDir <- "BPGs"
netSuffix <- ".bpg"
d1Suffix <- "Exp"
d2Suffix <- "CNV"
nodeNames <- read.table("/bi/aim/scratch/mkandula/Networks/DBs/KEGG_entrez_MSigDB_4_0.ids", stringsAsFactors=FALSE)[,1]

### Wmixnet params ###
threads <- 2

### Prepare data ###
# list all the input files (per patient, per data type)
d1Files <- list.files(paste0(netDir, "/"), pattern = paste0(d1Suffix, ".*", netSuffix))
d1Files <- sapply(strsplit(basename(d1Files),netSuffix), function(x) paste(x[1:(length(x)-1)]))
d2Files <- list.files(paste0(netDir, "/"), pattern = paste0(d2Suffix, ".*", netSuffix))
d2Files <- sapply(strsplit(basename(d2Files),netSuffix), function(x) paste(x[1:(length(x)-1)]))

sel <- 21 # length(d1Files)
allData <- list(d1=d1Files[1:sel], d2=d2Files[1:sel])
names(allData) <- c(d1Suffix, d2Suffix)
# initialize patient-patient networks
ppNet.G <- list(data.frame(p1=character(), p2=character(), vi=as.numeric()), data.frame(p1=character(), p2=character(), vi=as.numeric()))
ppNet.B <- list(data.frame(p1=character(), p2=character(), vi=as.numeric()), data.frame(p1=character(), p2=character(), vi=as.numeric()))
ppNet.P <- list(data.frame(p1=character(), p2=character(), vi=as.numeric()), data.frame(p1=character(), p2=character(), vi=as.numeric()))

names(ppNet.G) <- c(d1Suffix, d2Suffix)
names(ppNet.B) <- c(d1Suffix, d2Suffix)
names(ppNet.P) <- c(d1Suffix, d2Suffix)

### Set thresholds ###
#thr <- 0.01
thr <- 0

##########################
### --- END CONFIG --- ###
##########################

# testing

#d=1
#i <- sample(length(d1Files), 1)
#j=2

# force a specific no of groupings to make the patients comparable with the VI
minC=16
maxC=16

for (d in 1:length(allData)) {

  ## selecting the data type
  data <- allData[[d]]

  for (i in 1:(length(data)-1)) {

    ## first patient
    file <- data[i]
    p1 <- file

    # create an adjacency matrix from patient's bipartite graph
    dat1 <- as.matrix(read.table(paste0(netDir, "/", file, netSuffix), sep="\t", head=F, fill=TRUE, stringsAsFactors=FALSE))
    dat1 <- t(dat1) %*% dat1
    diag(dat1) <- 0 # get rid of self-loops
    # filter on a preselected threshold [optionally]
    # IDEA: thr could be mean(dat1)
    dat1[dat1 < thr] <- 0 #    dat1.A <- graph.adjacency(dat1, mode="undirected", weighted=TRUE, add.rownames=FALSE)

    # logit transform the data to make it analysable under Gaussian distro
    dat1 <- plogis(dat1)

    ## find communities in the first patient
    # Gaussian
    res1 <- BM_gaussian(membership_type='SBM', adj=dat1, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
    res1$estimate()

    selOpt <- which.max(res1$ICL)
    
    clustering1 <- as.data.frame(res1$memberships[[selOpt]]$Z, row.names=nodeNames)
    clustering1$cluster <- apply(clustering1, 1, function(x) as.numeric(which.max(x)))
    comm1.G <- clustering1$cluster
    
    js <- c((i+1):length(data))
    for (j in js) {

      ## next patient to be compared with the previous one
      file <- data[j]
      p2 <- file

      # create an adjacency matrix from patient's bipartite graph; filter on a preselected threshold [optionally]
      dat2 <- as.matrix(read.table(paste0(netDir, "/", file, netSuffix), sep="\t", head=F, fill=TRUE, stringsAsFactors=FALSE))
      dat2 <- t(dat2) %*% dat2
      diag(dat2) <- 0
      dat2[dat2 < thr] <- 0 #     dat2.A <- graph.adjacency(dat2, mode="undirected", weighted=TRUE, add.rownames=FALSE)

      # logit transform the data to make it analysable under Gaussian distro
      dat2 <- plogis(dat2)

      ## find communities in the next patient
      # Gaussian
      res2 <- BM_gaussian(membership_type='SBM', adj=dat2, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
      res2$estimate()
      clustering1 <- as.data.frame(res2$memberships[[selOpt]]$Z, row.names=nodeNames)
      clustering1$cluster <- apply(clustering1, 1, function(x) as.numeric(which.max(x)))
      comm2.G <- clustering1$cluster
      
      ## get a metric of similarity between the two community structures
      vi.G <- compare(comm1.G, comm2.G, method = "vi")

      ## save results into a patient-patient network
      ppRow.G <- data.frame(p1=p1, p2=p2, vi=vi.G)
      ppNet.G[[names(allData[d])]] <- rbind(ppRow.G, ppNet.G[[names(allData[d])]])
      
    }
    
  }
  
}
save.image()
save(ppNet, file="ppNets_v2.RData")









### test Gaussian distro - to figure out if we can use the clustering algorithm at all
adj <- dat1
rownames(adj) <- nodeNames
colnames(adj) <- nodeNames
spm <- write_spm(adj)

X11(); par(mfrow=c(4,4))
X11(); par(mfrow=c(2,2))
X11(); par(mfrow=c(5,4))
for (gr in 1:length(unique(clustering1$cluster))) {
  
  groupSel <- clustering1[clustering1$cluster==gr,]

  if (length(rownames(groupSel)) > 1) {
    
    groupPairs <- c()
    for (i in 1:(length(rownames(groupSel))-1)) {
      
      node1 <- rownames(groupSel)[i]
      
      js <- c((i+1):length(rownames(groupSel)))
      for (j in js) {
        
        node2 <- rownames(groupSel)[j]
        pair <- data.frame(node1=as.character(node1), node2=as.character(node2))
        groupPairs <- rbind(groupPairs, pair)
        
      }
      
    }
    
    groupPairs$weight <- 0

    for (i in 1:length(groupPairs$node1)) {
      
      node1 <- as.character(groupPairs$node1[i])
      node2 <- as.character(groupPairs$node2[i])
      
      groupPairs$weight[i] <- spm[spm$node1==node1 & spm$node2==node2, "weight"] # the order doesn't matter
      
    }
    
    hist(groupPairs$weight[groupPairs$weight!=0], main=paste0('force ', maxC,' gr; group ', gr, ' with ', length(rownames(groupSel)), ' unique nodes'), breaks=30)
#        hist(groupPairs$weight[groupPairs$weight!=0], main=paste0('opt ICL; group ', gr, ' with ', length(rownames(groupSel)), ' unique nodes'), breaks=30)
    
  }
}

norm <- 'logit.'
dev.copy2pdf(file=paste0("Plots/hist.",norm,maxC,"groups.patient359.pdf"))
dev.copy2pdf(file=paste0("Plots/hist.",norm,"optICL.patient359.pdf"))

### FUNCTIONS ###



write_spm <- function(adj_matrix)
{
    # This function take 1  argument
    # adj_matrix : adjacency matrix of the graph
    # and outputs a rewritten adjacency matrix

  spm <- c()
  n <- dim(adj_matrix)[1]
  v <- data.frame()

    for( i in 1:n )
    {
        for( j in 1:n )
        {
            if( i != j )
            {
                v<-data.frame(node1=as.character(rownames(adj_matrix)[i]),node2=as.character(colnames(adj_matrix)[j]),weight=as.numeric(adj_matrix[i,j]))
                spm <- rbind(spm, v)
              }
        }
    }

  return(spm)
  
}



## look for the optimal no of clusters for each distribution

patients <- sample(length(data), 10)

gaus <- c()
bern <- c()
pois <- c()

for (i in patients) {

  minC=1
  maC=Inf
  
  # Gaussian
  res1 <- BM_gaussian(membership_type='SBM', adj=dat1, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
  res1$estimate()
  gaus <- append(gaus, which.max(res1$ICL))
  
  # Bernoulli
  res1 <- BM_bernoulli(membership_type='SBM', adj=dat1, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
  res1$estimate()
  bern <- append(bern, which.max(res1$ICL))
  
  # Poisson
  res1 <- BM_poisson(membership_type='SBM', adj=dat1, verbosity=0, autosave='', plotting='', exploration_factor=1.5, explore_min=minC, explore_max=maxC, ncores=threads)
  res1$estimate()
  pois <- append(pois, which.max(res1$ICL))
  
}

mean(gaus)
mean(bern)
mean(pois)



### (B) Finding motifs (e.g., subgraphs) ###

D1 <- matrix(0, 5, 5)
D2 <- matrix(0, 5, 5)
D3 <- matrix(0, 5, 5)
D1[1:3, 1:3] <- 2
D2[3:5, 3:5] <- 3
D3[2:5, 2:5] <- 1
 
g <- simplify(graph.adjacency(D1 + D2 + D3, mode="undirected", weighted=TRUE))
V(g)$color <- "white"
E(g)$label <- E(g)$weight
E(g)$label.cex <- 2
E(g)$color <- "black"
layout(matrix(1:6, nrow=2, byrow=TRUE))
co <- layout.kamada.kawai(g)
par(mar=c(1,1,1,1))
plot(g, layout=co)
 
## Calculate graphlets
gl <- graphlets(g, niter=1000)
 
## Plot graphlets
for (i in 1:length(gl$cliques)) {
  sel <- gl$cliques[[i]]
  V(g)$color <- "white"
  V(g)[sel]$color <- "#E495A5"
  E(g)$width <- 1
  E(g)[ V(g)[sel] %--% V(g)[sel] ]$width <- 2
  E(g)$label <- ""
  E(g)[ width == 2 ]$label <- round(gl$Mu[i], 2)
  E(g)$color <- "black"
  E(g)[ width == 2 ]$color <- "#E495A5"
  plot(g, layout=co)
}
