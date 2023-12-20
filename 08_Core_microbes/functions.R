### scoring functions
func.score <- function(x){
  v <- 0
  if(sum(x[1:2]) == 2) v <- x[4]
  else if(sum(x[1:2]) == 1) v <- 1
  return(v*x[3])
}

func.redundancy.positive <- function(x, y, alpha, beta){
  i <- which(alpha >= 0)
  v <- 0
  if(length(i) > 0){
  z <- rbind(x[i], y[i], alpha[i], beta[i])
  v <- sum(apply(z, 2, func.score))
  }
  return(v)
}

func.redundancy.negative <- function(x, y, alpha, beta){
  i <- which(alpha < 0)
  v <- 0
  if(length(i) > 0){
  z <- rbind(x[i], y[i], alpha[i], beta[i])
  v <- sum(apply(z, 2, func.score))
  }
  return(v)
}

func.balance <- function(x, y, alpha, gamma, delta){
  i <- which(alpha >= 0)
  v <- 0
  if(length(i) > 0){
    v <- 1 - gamma * abs(sum((delta + (1 - delta)*alpha[i]) * (x[i] - y[i])))/sum((delta + (1 - delta)*alpha[i]))
  }
  return(v)
}

score.cal <- function(x, y, alpha, beta, gamma, delta){
  v <- func.redundancy.positive(x, y, alpha, beta) * func.balance(x, y, alpha, gamma, delta) +
    func.redundancy.negative(x, y, alpha, beta)
  return(v)
}


func.btw.final <- function(g, f, alpha, beta, gamma = 1, delta = 1){
  ### improved betweeness centrality
  n <- vcount(g) # number of nodes
  m <- ncol(f) # number of function types
  b <- numeric(n)
  h <- matrix(0, n, m)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      r <- get.all.shortest.paths(g, i, j)
      if (length(r$res) > 0) {
        fscore <- score.cal(f[i, ], f[j, ], alpha, beta, gamma, delta)
        for (k in 1:length(r$res)) {
          l <- length(r$res[[k]])
          if (l > 1){
            node.id <- as.numeric(r$res[[k]][c(-1, -l)])
            b[node.id] <- b[node.id] + 1 / length(r$res) * fscore
            h[node.id, ] <- h[node.id, ] + matrix(rep(unlist(f[i,]), length(node.id)), ncol = m, byrow = T) + matrix(rep(unlist(f[j,]), length(node.id)), ncol = m, byrow = T)
          }
        }
      }
    }
  }
  result <- list(b, h)
  return(result)
}



### original betweenness centrality
func.btw.original <- function(g){
  ### improved betweeness centrality 
  n <- vcount(g) # number of nodes
  b <- numeric(n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      r <- get.all.shortest.paths(g, i, j)
      if (length(r$res) > 0) {
        for (k in 1:length(r$res)) {
          l <- length(r$res[[k]])
          if (l > 1) b[r$res[[k]][c(-1, -l)]] <- b[r$res[[k]][c(-1, -l)]] + 1 / length(r$res)
        }
      }
    }
  }
  return(b)
}

### community detection
community.imp <- function(g, n = NULL){
  data <- edge.betweenness.community(g)
  
  g.Q <- rep(0, vcount(g))
  for (i in 1:vcount(g)){
    memb <- cut_at(data, i)
    g.Q[i] <- modularity(g, memb)
  }
  opt.n <- which(g.Q == max(g.Q))
  
  if(is.null(n)){
    x <- cut_at(data, opt.n)
    plot(g.Q[1:(opt.n+5)], type = "l", xlab = "number of modules", ylab = "modularity")
    abline(v = opt.n, col = "blue", lty = 2)
  }else{
    x <- cut_at(data, n)
    plot(g.Q[1:(opt.n+5)], type = "l", xlab = "number of modules", ylab = "modularity")
    abline(v = opt.n, col = "blue", lty = 2)
    abline(v = n, col = "red", lty = 3)
  }
  return(x)
}

