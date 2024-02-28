rggd.fit<-function (xdat, r = dim(xdat)[2],n=20,ydat = NULL, mul = NULL, sigl = NULL, hl= NULL, 
                     mulink = identity, siglink = identity, shlink = identity, hlink = identity, show = TRUE, 
                     method = "Nelder-Mead", maxit = 10000, ...) 
{
  
  options(digits=8)
  z <- list()
  k <- list()
  
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  nph  <- length(hl) + 1
  
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
  }
  
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
  }
  
  if (is.null(hl)) {
    hmat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    hmat <- cbind(rep(1, dim(xdat)[1]), ydat[, hl])
  }
  
  xdatu <- xdat[, 1:r, drop = FALSE]
  init <- ginit2(xdatu,n=n)[,1:3]
  
  z$link <- deparse(substitute(c(mulink, siglink, hlink)))
  
  z1 <- apply(xdatu, 1, max, na.rm = TRUE) # z_1
  zr <- apply(xdatu, 1, min, na.rm = TRUE) # z_r
  
  ggd.lik <- function(a) {
    
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    h  <- hlink(hmat %*% (a[seq(npmu + npsc + 1, length = nph)]))
    
    y <- exp( - (xdat - mu)/sc )
   
    f <- 1 - h * exp(-(xdat - mu)/sc)
      
    if(max(f,na.rm=T)<0) return(10^6)
      
    if(any(y<= 0,na.rm=T) || any(sc<= 0,na.rm=T) || any(f^(1/h)<0,na.rm=T) || any(f^(1/h)>1,na.rm=T))  
        return(10^6)
    sum(log(sc)) - sum(log(y)) - sum(((1-h)/h)*log(f))

   }
  
  rggd.lik <- function(a) {
    
    mu <- mulink(drop(mumat %*% (a[1:npmu])))
    sc <- siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    h  <- hlink(drop(hmat %*% (a[seq(npmu + npsc + 1, length = nph)])))
    
    ## constraints ##
    
    if (r>=2){
      
      if(min(h) > (1/(r-1))) return(10^6)
      
    }
    
    ri  <- (r-seq(1:(r))) # r-i
    cr  <- (1-ri*h[1])    # c_r
    
    if (any(sc <= 0) || any(cr < 0) ) return(10^6) 
    
      y <- exp(-(xdatu - mu)/sc)
      
      f <- 1 - h * exp(-(zr - mu)/sc)
      
      if(any(f<0,na.rm=T)) return(10^6)
      
      if(any(y<= 0,na.rm=T) || any(sc<= 0,na.rm=T) || any(f^(1/h)<0,na.rm=T) || any(f^(1/h)>1,na.rm=T))  
        return(10^6)
      
      # constraints 2 #
      
      if (any( min(h) > 1/exp(-(zr - mu)/sc), na.rm= TRUE)) return(10^6)
      
      # constraints 2,3,4 #
      
      
      y <- log(sc) - log(y) - log(cr)  
      y <- rowSums(y, na.rm = TRUE)
      
      sum((r*h - 1)/h * log(f) + y)
      
  }
  
  if(r==1){
    
    tryCatch(
      
      for(i in 1:nrow(init)){       
        
        value <- try(optim(init[i,], ggd.lik,method = method, hessian = TRUE))
        
        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
        }
        
      }
    )
    
  }else{
    
    tryCatch(
      
      for(i in 1:nrow(init)){       
        
        value <- try(optim(init[i,], rggd.lik, method = method, hessian = TRUE))

        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
        }
        
      }
    )
    
  }
  
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, function(x) x$value[which.min(x$value)]))# %>% filter(value!=10^6)
  optim_value  <-optim_value[optim_value$value!=10^6,]
  
  
  if(r==1){optim_grad   <-sapply(optim_value$num, function(x) sum(abs(grad(ggd.lik,k[[x]]$par))))
  }else   {optim_grad   <-sapply(optim_value$num, function(x) sum(abs(grad(rggd.lik,k[[x]]$par))))}
  
  optim_value$grad <- optim_grad
  
  optim_table1 <-optim_value[order(optim_value$grad, optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <-k[[selc_num]]
  
  mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
  sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
  h  <- hlink(drop(hmat %*% (x$par[seq(npmu + npsc + 1, length = nph)])))
  
  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$data <- xdat
  z$mle  <- x$par
  
  if(r==1){z$grad <- grad(ggd.lik,x$par)
  }else   {z$grad <- grad(rggd.lik,x$par)}
  
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, h)
  z$r <- r
  
  z$rl20  <- ggd.rl_m(z$mle,z$cov,year=20)
  z$rl50  <- ggd.rl_m(z$mle,z$cov,year=50)
  z$rl100 <- ggd.rl_m(z$mle,z$cov,year=100)
  
  nr   <-nrow(xdat) 
  
  options(digits=8)
  
  z$rslt <- round(data.frame(r = r,
                             nllh = z$nllh,
                             tr_V = tr(z$cov),
                             logdet = log(det(z$cov)),
                             aic  = cal_stat(z$nllh,3,method="AIC"),
                             bic  = cal_stat(z$nllh,3,nr,method="BIC"),
                             mle  = t(z$mle),
                             se   = t(z$se),
                             rl20   = z$rl20$rl,
                             rl20se = z$rl20$rl_se,
                             rl50   = z$rl50$rl,
                             rl50se = z$rl50$rl_se,
                             rl100   = z$rl100$rl,
                             rl100se = z$rl100$rl_se),3)
  
  z$g   <-  data.frame(r=r,
                       nllh=z$nllh,
                       g1 = z$grad[1],
                       g2 = z$grad[2],
                       g3 = z$grad[3]) 
  
  return(z)
  invisible(z)
}