require(igraph)

sna_orange <- rgb(255, 153, 0, maxColorValue=255)
sna_green <- rgb(153, 204, 0, maxColorValue=255)
grey <- rgb(200, 200, 200, maxColorValue=255)
black <- rgb(0, 0, 0, maxColorValue=255)
blue <- rgb(153,204,255, maxColorValue=255)
white <- rgb(255,255,255, maxColorValue=255)



actornetwork <- function(df){ # network data frame
  colnames(df) <- c("x","y","z")
  df$z[df$z<0] <- -1 # binarize
  df$z[df$z>0] <- 1
  dn <- df[df$z < 0,] # negative subnetwork
  dp <- df[df$z > 0,] # positive subnetwork
  n <- graph.data.frame(dn, directed=TRUE)
  V(n)$type <- bipartite.mapping(n)$type;
  V(n)$degree1mode <- degree(n);
  proj <- bipartite.projection(n);
  g1 <- proj[[1]];
  E(g1)$w2 <- -E(g1)$weight;
  g1e <- cbind(get.edgelist(g1,names=TRUE), E(g1)$w2);
  p <- graph.data.frame(dp, directed=TRUE)
  V(p)$type <- bipartite.mapping(p)$type;
  V(p)$degree1mode <- degree(p);
  proj <- bipartite.projection(p);
  g2 <- proj[[1]];
  E(g2)$w2 <- E(g2)$weight;
  g2e <- cbind(get.edgelist(g2,names=TRUE), E(g2)$w2);
  # Both graphs combined
  g <- graph.data.frame(data.frame(rbind(g1e,g2e)), directed=FALSE);
  E(g)$w2 <- as.numeric(E(g)$X3);
  E(g)$weight <- abs(as.numeric(E(g)$X3));
  return(g)
}

conceptnetwork <- function(df){ # network data frame
  colnames(df) <- c("x","y","z")
  df$z[df$z<0] <- -1 # binarize
  df$z[df$z>0] <- 1
  dn <- df[df$z < 0,] # negative subnetwork
  dp <- df[df$z > 0,] # positive subnetwork
  n <- graph.data.frame(dn, directed=TRUE)
  V(n)$type <- bipartite.mapping(n)$type;
  V(n)$degree1mode <- degree(n);
  proj <- bipartite.projection(n);
  g1 <- proj[[2]];
  E(g1)$w2 <- -E(g1)$weight;
  g1e <- cbind(get.edgelist(g1,names=TRUE), E(g1)$w2);
  p <- graph.data.frame(dp, directed=TRUE)
  V(p)$type <- bipartite.mapping(p)$type;
  V(p)$degree1mode <- degree(p);
  proj <- bipartite.projection(p);
  g2 <- proj[[2]];
  E(g2)$w2 <- E(g2)$weight;
  g2e <- cbind(get.edgelist(g2,names=TRUE), E(g2)$w2);
  # Both graphs combined
  g <- graph.data.frame(data.frame(rbind(g1e,g2e)), directed=FALSE);
  E(g)$w2 <- as.numeric(E(g)$X3);
  E(g)$weight <- abs(as.numeric(E(g)$X3));
  return(g)
}


actorconflictnetwork <- function(df){ # network data frame
  df[,3] <- ifelse(df[,3]<0, -1, 1)   # binarize data frame
  colnames(df) <- c("x","y","z")      # set standard column names
  m <- xtabs(z~x+y, data=df)          # create matrix from data frame
  # positive links
  m_p <- m                            # create matrix containing only the positive links
  m_p[m_p < 0] <- 0
  g1m <- m_p %*% t(m_p)               # bipartite projection for this matrix
  diag(g1m) <- 0                      # remove loops
  g1m[lower.tri(g1m)] <- 0            # remove redundant information
  g1 <- graph.adjacency(g1m, weighted = TRUE) # create network from matrix
  E(g1)$coalition <- "pos"            # set coalition value
  E(g1)$color <- "green"              # set edge color for this coalition
  g1e <- cbind(get.edgelist(g1,names=TRUE), E(g1)$weight, E(g1)$coalition, E(g1)$color); # create edgelist from network
  # negative links
  m_n <- m
  m_n[m_n > 0] <- 0
  g2m <- m_n %*% t(m_n)
  diag(g2m) <- 0
  g2m[lower.tri(g2m)] <- 0
  g2 <- graph.adjacency(g2m, weighted = TRUE)
  E(g2)$coalition <- "neg"
  E(g2)$color <- "red"
  g2e <- cbind(get.edgelist(g2,names=TRUE), -E(g2)$weight, E(g2)$coalition, E(g2)$color);
  # conflict links
  g3m <- (m_p %*% t(m_n) + m_n %*% t(m_p))*-1 # sum of the two disagreeement projections
  diag(g3m) <- 0
  g3m[lower.tri(g3m)] <- 0
  g3 <- graph.adjacency(g3m, weighted = TRUE)
  E(g3)$coalition <- "conf"
  E(g3)$color <- "grey"
  g3e <- cbind(get.edgelist(g3,names=TRUE), E(g3)$weight, E(g3)$coalition, E(g3)$color);
  # all graphs combined
  df_all <- data.frame(rbind(g1e,g2e,g3e)) # combine the three dataframes
  colnames(df_all) <- c("actor", "actor2", "weight", "coalition", "color") # set meaningful column names
  g <- graph.data.frame(df_all, directed=FALSE); # create graph from combined data frame
  E(g)$strength <- E(g)$weight          # re-create signed and unsigned weight values
  E(g)$weight <- abs(as.numeric(E(g)$strength))
  return(g)
}

conceptconflictnetwork <- function(df){ # network data frame
  # negative
  # create adjacency matrix
  colnames(df) <- c("x","y","z")
  m <- with(df, {
    out <- matrix(nrow=nlevels(x), ncol=nlevels(y),
                  dimnames=list(levels(x), levels(y)))
    out[cbind(x, y)] <- z
    out
  })
  m[is.na(m)] <- 0
  # positive
  m_p <- m
  m_p[m_p < 0] <- 0
  g1m <- t(m_p) %*% m_p
  diag(g1m) <- 0
  g1m[lower.tri(g1m)] <- 0
  g1 <- graph.adjacency(g1m, weighted = TRUE)
  E(g1)$coalition <- "pos"
  g1e <- cbind(get.edgelist(g1,names=TRUE), E(g1)$weight, E(g1)$coalition);
  # negative
  m_n <- m
  m_n[m_n > 0] <- 0
  g2m <- t(m_n) %*% m_n
  diag(g2m) <- 0
  g2m[lower.tri(g2m)] <- 0
  g2 <- graph.adjacency(g2m, weighted = TRUE)
  E(g2)$coalition <- "neg"
  g2e <- cbind(get.edgelist(g2,names=TRUE), -E(g2)$weight, E(g2)$coalition);
  # conflict
  g3m <- t(m) %*% m
  g3m[g3m > 0] <- 0
  g3m <- abs(g3m)
  diag(g3m) <- 0
  g3m[lower.tri(g3m)] <- 0
  g3 <- graph.adjacency(g3m, weighted = TRUE)
  E(g3)$coalition <- "conf"
  g3e <- cbind(get.edgelist(g3,names=TRUE), E(g3)$weight, E(g3)$coalition);
  # all graphs combined
  df_all <- data.frame(rbind(g1e,g2e,g3e))
  colnames(df_all) <- c("actor", "actor2", "weight", "coalition")
  g <- graph.data.frame(df_all, directed=FALSE);
  #  g <- simplify(g, remove.loops=TRUE, remove.multiple = FALSE)
  E(g)$strength <- E(g)$weight
  E(g)$weight <- abs(as.numeric(E(g)$strength))
  return(g)
}

nslice <- function(n,i){ # network (n) and slice (i)
  n <- delete.edges(n, E(n)[weight < i]) # 3-sclice;
  V(n)$degree <- degree(n);
  n <- delete.vertices(n, V(n)[degree == 0]);
  return(n)
}

nsplit <- function(n){ # network (n)
  neg <- delete.edges(n, E(n)[Agreement > 0]) # negative subnet
  V(neg)$degree <- degree(neg);
  neg <- delete.vertices(neg, V(neg)[degree == 0]);
  pos <- delete.edges(n, E(n)[Agreement < 0]) # positive subnet
  V(pos)$degree <- degree(pos);
  pos <- delete.vertices(pos, V(pos)[degree == 0]);
  n <- list("neg"=neg, "pos"=pos)
  return(n)
}

twoModeAttributes <- function(n){ # network (n)
  V(n)$type <- bipartite.mapping(n)$type;
  V(n)$degree <- degree(n, mode="all");
  V(n)$indegree <- degree(n, mode="in");
  V(n)$outdegree <- degree(n, mode="out");
  V(n)$id <- V(n)$name;
  V(n)[indegree>0]$shape <- "square";
  V(n)[outdegree>0]$shape <- "circle";
  V(n)[indegree>0]$color <- blue;
  V(n)[outdegree>0]$color <- white;
  E(n)[Agreement<0]$color <- sna_orange;
  E(n)[Agreement>0]$color <- sna_green;
  return(n)
}

#################
#################
