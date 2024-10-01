rm( list=ls() )

niter <- 1000
cor1 <- function(x,y) cor( x, y+rnorm(length(y))*1e-12)

pred1 <- function(y,g,e,trn){
  betas  <- coef( lm( y ~ g + e, subset=trn ) )
  betas[1] + g[-trn]*betas[2] + e[-trn]*betas[3]
}

pred2 <- function(y,g,e,trn){
  betas  <- coef( lm( y ~ g + e + g:e, subset=trn ) )
  betas[1] + g[-trn]*betas[2] + e[-trn]*betas[3] + (g*e)[-trn]*betas[4]
}

pred3 <- function(y,g,e,trn,f,finv){
  ystar <- finv(y)
  betas  <- coef( lm( ystar ~ g + e, subset=trn ) )
  f(betas[1] + g[-trn]*betas[2] + e[-trn]*betas[3])
}

pred4 <- function(y,g,e,trn,f,finv,b,y0){
  yhat   <- -y0 + g[-trn]*sqrt(b) + e[-trn]*sqrt(b)
  f(yhat)
}

bn <- 21
bs <- seq( 0, 1/2, len=bn )

for( N in c( 20, 100 ) ){
  pdf( paste0( '~/Desktop/scale-prediction_N=', N, '.pdf' ), w=12, h=6 )

par( mfcol=c(2,4) )
for( ftype in 1:4 )
  for( xtype in 1:2 )
    {
   
  if( ftype == 1 ){
    f <- function(x) x
    finv <- function(x) x
    main='No Transform'
  } else if( ftype==2 ){
    f <- function(x) x^2
    finv <- function(x) sqrt(x)
    main='Quadratic Transform'
  } else if( ftype==3 ){
    f <- function(x) exp(x)
    finv <- function(x) log(x)
    main='Exponential Transform'
  } else if( ftype==4 ){
    f <- function(x) exp(x)/(1+exp(x))
    finv <- function(x) log(x/(1-x))
    main='Inverse-Logit Transform'
  }
  
  if( xtype == 1 ){
    rg <- function(x) rnorm(N)
    re <- function(x) rnorm(N)
  } else {
    rg <- function(x) scale( rbinom(N,1,.5) )
    re <- function(x) scale( rbinom(N,1,.5) )
  }
  
  corrs <- array( NA, dim=c(4, bn, niter) )
  for( it in 1:niter )
    for( bi in 1:bn )
    {
      g <- rg()
      e <- rg()
      ystar <- g * sqrt(bs[bi]) + e * sqrt(bs[bi]) + rnorm(N)*sqrt(1-2*bs[bi])
      y0 <- min( ystar )
      ystar <- ystar - y0
      y <- f( ystar )
      
      ytst <- y[-(1:(N/2))]
      corrs[1,bi,it] <-  cor1( ytst, pred1(y,g,e,trn=1:(N/2)) )
      corrs[2,bi,it] <-  cor1( ytst, pred2(y,g,e,trn=1:(N/2)) )
      corrs[3,bi,it] <-  cor1( ytst, pred3(y,g,e,trn=1:(N/2),f,finv) )
      corrs[4,bi,it] <-  cor1( ytst, pred4(y,g,e,trn=1:(N/2),f,finv,bs[bi],y0) )
    }
  
  xs <- 2*bs
  ys <- apply( corrs, 1:2, mean, na.rm=T )^2
  plot( 0:1, 0:1, type='n', main=main, xlab='Var Explained by G+E on additive scale', ylab='Prediction R2' )
  for( j in 4:1 )
    lines( xs, ys[j,], col=j, lwd=2 )
  if(ftype==1 & xtype==1)
    legend( 'bottomright', bty='n', fill=1:4, title='Predictor', leg=c( 'Add', 'GxE', 'Descale', 'Oracle' ) )
  
}
dev.off()
}
