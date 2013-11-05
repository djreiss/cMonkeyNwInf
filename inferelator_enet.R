###########################################
# cv.glmnet
###########################################
# fxn derived from CRAN pkg 'elasticnet'
# to do CV on results from glmnet since
# glmnet itself doesn't do CV   (note as of version 1.4 it does have a CV function!)
# requires some internal fxns from pkg 'lars'
cv.glmnet <- function(x, y, lambda, K = 10, cv.reps=10, trace = FALSE, plot.it = TRUE, se = TRUE, weights=NA, ...) { ##lambda, s, 
  all.folds <- do.call( "c", lapply( 1:cv.reps, function( i ) cv.folds(length(y), K) ) )
  ##all.folds <<- cv.folds(length(y), K)
  ##print(length(all.folds))##;stop()
##  names(all.folds) <- 1:length( all.folds )
  ##residmat <- matrix(0, length(lambda), K * cv.reps)
  ##for(i in 1:length( all.folds ) ) { ##seq(K)) {
  apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
  if ( is.na( weights[ 1 ] ) ) weights <- rep( 1, length( y ) )
  residmat <- do.call( cbind, apply.func( 1:length( all.folds ), function( i ) { 
    if(trace) cat("CV Fold", i, "\n")
    omit <- all.folds[[i]]
    fit <- my.glmnet(x[ -omit,  ], y[ -omit], lambda=lambda, weights=weights[ -omit ], ... ) ##s=s,...)
    fit <- predict(fit, x[omit,  ,drop=FALSE], ...)
    if(length(omit)==1){fit<-matrix(fit,nrow=1)}
    ##residmat[, i] <-
    apply((y[omit] - fit)^2, 2, mean)
  } ) )
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(cv = cv, cv.error = cv.error, lambda=lambda, fraction=log(lambda)) ##s = s,  (fraction is for plotting via plotCVLars)
  if(plot.it) {
     plot(lambda, cv, type = "b", ylim = range(cv, cv + cv.error, cv - cv.error), ...)
     if (se) error.bars(lambda, cv + cv.error, cv - cv.error, width = 1/length(lambda)) ##length(s)
  }
  invisible(object)
}

## Stoopid - allow glmnet to not barf on ...'s
my.glmnet <- function(x, y, family=c("gaussian","binomial","poisson","multinomial","cox")[1],
                      weights, offset=NULL, alpha=1, nlambda=100,
                      lambda.min=1e-6, ##ifelse(nrow(x) < ncol(x), 0.05, 1e-04),
                      lambda=NULL, standardize=TRUE, intercept=TRUE, thresh=1e-04, 
                      dfmax=ncol(x) + 1, pmax=min(dfmax * 1.2, ncol(x)), exclude,
                      penalty.factor=rep(1, ncol(x)), lower.limits=-Inf, upper.limits=Inf,
                      maxit=100, type.gaussian=ifelse(ncol(x)<500,"covariance","naive"),
                      type.logistic=c("Newton","modified.Newton"), standardize.response=FALSE,
                      type.multinomial=c("ungrouped","grouped"),
                      ... ) ##c("covariance", "naive"), HessianExact=FALSE, 
  glmnet( x, y, family, weights, offset, alpha, nlambda, lambda.min, lambda, standardize, intercept, thresh, dfmax, pmax,
         exclude, penalty.factor, lower.limits, upper.limits, maxit, type.gaussian, type.logistic, standardize.response,
         type.multinomial ) ##HessianExact, 

inferelator.enet <- function( profile, predictor.mat, conds.use, col.map=NULL, tau=10, ##tau.best=NULL,
                             ratio.cutoff=3, coef.cutoff=0.02, cv.k=10, cv.choose="min+2se", ##logit.scale=0.25,
                             n.boot=1, boot.opt=c( "resample", "cv" ), rescale.coeffs=T, quiet=T, alpha=0.9, ##alpha="cv.choose",
                             ##alphas=seq( 0, 1, by=0.05 ), ## 0.5, 1, by=0.05 )
                             weights=NA, penalties=NA, max.coeffs=NA, min.coeffs=NA, ... ) { ##plot=T, 

  ##require( glmnet, quietly=T )

  ## this function is the MEAT - run lars and output coefficients	
  if ( cv.choose == "min" ) cv.choose <- "min+0se"

  ##profilez<<-profile;predictor.matz<<-predictor.mat;condz.use<<-conds.use;col.mapz<<-col.map
  ## Convert predicted variable to the difference equation, i.e. from y to tau/dt * (y_i - y_(i-1)) + y_(i-1)
  tmp <- get.input.matrix( profile, predictor.mat, conds.use, col.map=col.map, tau=tau, ratio.cutoff=ratio.cutoff,
                          quiet=quiet )

  ##tmpz<<-tmp
##   out.tmp <- profile
##   in.tmp <- predictor.mat
##   if ( ! is.null( col.map ) && ! is.na( tau ) && tau > 0 ) {
##     if ( ! quiet ) cat( "Time series data supplied, converting predictors to difference equation... \n" )
##     if ( ! quiet ) cat( "Tau =", tau, "\n\n" )
##     conds <- colnames( predictor.mat ) ##colnames( data )
##     cm <- col.map[ conds, ]
##     good.i <- ( ( cm$isTs == TRUE ) & ( cm$is1stLast %in% c( "m", "l" ) ) ) | ##| col.map[ conds, "" ] == "l" ) ) |
##                 ( cm$isTs == FALSE & cm$is1stLast == "e" )
##     curs <- as.character( col.map$condName[ good.i ] )
##     prevs <- as.character( col.map$prevCol[ good.i ] ) ##[ conds, "prevCol" ]
##     del.ts <- as.numeric( as.character( col.map$delta.t[ good.i ] ) ) ##[ conds, "delta.t" ]
##     del.ts[ del.ts < 1 ] <- 1
##     out.tmp <- ( ( tau / del.ts ) * ( out.tmp[ curs ] - out.tmp[ prevs ] ) ) + out.tmp[ prevs ]
##     out.tmp[ out.tmp > ratio.cutoff ] <- ratio.cutoff
##     out.tmp[ out.tmp < -ratio.cutoff ] <- -ratio.cutoff
##     in.tmp <- predictor.mat[ ,prevs ]
##     colnames( in.tmp ) <- names( out.tmp )
  
##   } else {
##     if ( ! quiet ) cat( "NO time series data supplied, assuming STEADY-STATE... \n\n" )
##   }

## ########### FIND PARAMS USING LARS

##   df.tmp <- t( in.tmp[ ,names( out.tmp ) %in% conds.use ] ) ##predictor.mat[ predictors, ] )
##   output <- as.numeric( out.tmp[ names( out.tmp ) %in% conds.use ] ) ##predictor.mat[ predicted[ 1 ], ] ) ### Linear output
  df.tmp <- tmp$inp[ ! is.na( tmp$outp ), ]
  output <- tmp$outp[ ! is.na( tmp$outp ) ] ## just in case (thanks to SD)
  rm( tmp )

  orig.penalties <- penalties
  in.penalties <- rep( 1, ncol( df.tmp ) ); names( in.penalties ) <- colnames( df.tmp )
  if ( exists( "penalties" ) && ! is.na( penalties ) && any( names( penalties ) %in% names( in.penalties ) ) ) {
    penalties <- penalties[ names( penalties ) %in% names( in.penalties ) ]
    in.penalties[ names( penalties ) ] <- penalties
  } else {
    orig.penalties <- in.penalties
  }
  
  ## TODO: if penalty[x] is not 1, set any penalty[x.y.min] to the average of penalty[x] and penalty[y] ??
  if ( any( orig.penalties != 1 ) ) {
    warning( "PENALTIES were set to non-1!!!" )
    tmp <- names( which( orig.penalties != 1 ) )
    for ( t in tmp ) {
      g <- grep( t, names( in.penalties ), val=T, fixed=T )
      g <- grep( "~~", g, val=T, fixed=T )
      if ( length( g ) <= 0 ) next
      gg <- strsplit( g, "~~", fixed=T )
      for ( i in 1:length( g ) ) {
        p1 <- orig.penalties[ gg[[ i ]][ 1 ] ]; if ( is.na( p1 ) ) p1 <- 1
        p2 <- orig.penalties[ gg[[ i ]][ 2 ] ]; if ( is.na( p2 ) ) p2 <- 1
        in.penalties[ g[ i ] ] <- mean( c( p1, p2 ), na.rm=T )
      }
    }
    rm( orig.penalties, g, gg, p1, p2, t, tmp )
  }

  if ( is.na( weights[ 1 ] ) ) weights <- rep( 1, length( output ) )
  else weights <- weights[ names( output ) ]
  names( weights ) <- names( output )

  ##df.tmp<<-df.tmp;output<<-output;weights<<-weights
  ##tmpz<<-cbind(apply(df.tmp,2,cor,output),apply(predictor.mat[colnames(df.tmp),rownames(df.tmp)],1,cor,profile[rownames(df.tmp)]))

  ##out.coe <- list() ##numeric()
  ##for ( boot in 1:n.boot ) {
  apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
  if ( ! quiet ) cat( "Alpha =", alpha, "\n" )
  out.coe <- apply.func( 1:n.boot, function( boot ) {
    cols <- 1:length( output )
    if ( boot > 1 && boot.opt == "resample" ) cols <- sample( cols, replace=T )

    ##tmpz<<-list(df.tmp=df.tmp,cols=cols,output=output,in.penalties=in.penalties,weights=weights,alpha=alpha)
### Run glmnet on predicted ~ predictors
    glmnet.obj <- ##try(
      my.glmnet( df.tmp[ cols, ], output[ cols ], penalty.factor=in.penalties, weights=weights[ cols ],
                                 alpha=if ( alpha == "cv.choose" ) 0 else alpha, ... )##, silent=quiet )
    if ( "try-error" %in% class( glmnet.obj ) ) {
      tries <- 1; while( tries <= 20 && "try-error" %in% class( glmnet.obj ) ) {
        glmnet.obj <- try( my.glmnet( df.tmp[ cols, ], output[ cols ], penalty.factor=in.penalties[ cols ],
                                     weights=weights[ cols ], alpha=if ( alpha == "cv.choose" ) 0 else
                                     alpha, ... ), silent=quiet )
        tries <- tries + 1
      } }
    if ( "try-error" %in% class( glmnet.obj ) ) return( numeric() )

### Do glmnet CV and ...
##     if ( alpha == "cv.choose" ) { ## Looks like alpha=0.9 is best in most cases but this is unstable for some...
##       ## Try multiple alphas and choose alpha that minimizes the cv error... then re-run glmnet/cv.glmnet with this
##       ##   alpha (below)
##       ##alphas <- seq( 0.25, 1, by=0.01 ) ## 0.5, 1, by=0.05 )
##       lambdas <- glmnet.obj$lambda; lambdas <- lambdas[ lambdas > 1e-8 ]
##       apply.func <- if ( multicore:::isChild() ) lapply else mclapply
##       se <- 0; if ( grepl( "+", cv.choose[ 1 ], fixed=T ) )
##         se <- as.numeric( gsub( "min+", "", gsub( "se", "", cv.choose[ 1 ] ) ) )
##       ##cv.mat <- cv.err.mat <- 0
##       ##n.iter <- 100
##       ##for ( iter in 1:n.iter ) {
##       cv.glmnet.obj <- apply.func( alphas, function( alph )
##                                   try( cv.glmnet( df.tmp[ cols, ], output[ cols ], lambda=lambdas, K=cv.k,
##                                                  cv.reps=100, alpha=alph, plot.it=F, trace=F ), silent=quiet ) )
##       cv.mat <- sapply( cv.glmnet.obj, "[[", "cv" )
##       cv.err.mat <- sapply( cv.glmnet.obj, "[[", "cv.error" )
##       if ( TRUE ) { ## plotting stuff
##         arr.min <- which( cv.mat == min( cv.mat ), arr=T )
##         alpha <- alphas[ arr.min[ 2 ] ]
##         min.lam <- lambdas[ arr.min[ 1 ] ]
##         par( mfrow=c( 2, 2 ) )
##         grey.image( cv.mat[ nrow( cv.mat ):1, ], x=log10( rev( lambdas ) ), y=alphas, xlab="log10 lambda",
##                    ylab="alpha", col=topo.colors( 256 ) )
##         best.s <- t( sapply( 1:ncol( cv.mat ), function( i ) {
##           min.i <- which.min( cv.mat[ ,i ] )
##           min.err <- cv.err.mat[ min.i, i ]
##           best.s <- min( which( cv.mat[ ,i ] <= min( cv.mat[ ,i ] ) + se * min.err ) )

##           glmnet.obj <- try( my.glmnet( df.tmp[ cols, ], output[ cols ], lambda=lambdas,
##                             alpha=alphas[ i ] ), silent=quiet )
##           coeffs.tmp <- as.matrix( coef( glmnet.obj, s=glmnet.obj$lambda[ best.s ] ) )
##           coeffs <- coeffs.tmp[ abs( coeffs.tmp ) >= coef.cutoff, ]
##           ##cat(i,best.s,glmnet.obj$lambda[ best.s ],"\n");print(coeffs)
          
##           c( best.s, cv.mat[ best.s, i ], cv.mat[ min.i, i ], length( coeffs ) )
##         } ) )
##         best.s<<-best.s;cv.glmnet.obj<<-cv.glmnet.obj;cv.mat<<-cv.mat

##         points( log10( min.lam ), alpha, col="yellow", pch=19 )
##         points( log10( lambdas[ best.s[ ,1 ] ] ), alphas, col="green", pch=19 )
##         points( log10( min.lam ), alpha, col="red", pch=19 )
##         plot( alphas, best.s[ ,2 ] )
##         plot( alphas, best.s[ ,4 ] )
##         plot( best.s[ ,2 ], best.s[ ,4 ] )
##         cat(alpha,min.lam,"\n")
##         ##}
##         stop()
##       }
      
##       cv.mat <- cv.mat / n.iter
##       cv.err.mat <- cv.err.mat / n.iter
##       ##cv.glmnet.obj<<-cv.glmnet.obj;cv.mat<<-cv.mat;cv.err.mat<<-cv.err.mat;alphas<<-alphas;lambdas<<-lambdas

##       arr.min <- which( tmp.arr == min( tmp.arr ), arr=T )
##       alpha <- alphas[ arr.min[ 2 ] ]
##       min.lam <- lambdas[ arr.min[ 1 ] ]
##       if ( TRUE ) { ## plotting stuff
##         grey.image( cv.mat[ nrow( cv.mat ):1, ], x=log10( rev( lambdas ) ), y=alphas, xlab="log10 lambda",
##                    ylab="alpha", col=topo.colors( 256 ) )
##         contour( log10( rev( lambdas ) ), alphas, cv.mat[ nrow( cv.mat ):1, ], nlev=10, col="red", add=T )
##         points( log10( min.lam ), alpha, col="red", pch=19 )
##         x11()
##       }
##       if ( ! quiet ) cat( "Using alpha =", alpha, "chosen via cross-validation with lambda =", min.lam,
##                          ", and", cv.choose[ 1 ], ".\n" )
##       glmnet.obj <- try( my.glmnet( df.tmp[ cols, ], output[ cols ], ... ), silent=quiet )
##       if ( "try-error" %in% class( glmnet.obj ) ) {
##         tries <- 1; while( tries <= 20 && "try-error" %in% class( glmnet.obj ) ) {
##           glmnet.obj <- try( my.glmnet( df.tmp[ cols, ], output[ cols ], alpha=alpha, ... ), silent=quiet )
##           tries <- tries + 1
##         } }
##       if ( "try-error" %in% class( glmnet.obj ) ) return( numeric() )
##     }

    cv.glmnet.obj <- try( cv.glmnet( df.tmp[ cols, ], output[ cols ], lambda=glmnet.obj$lambda, K=cv.k, trace=F,
                                    penalty.factor=in.penalties, weights=weights[ cols ],
                                    alpha=if ( alpha == "cv.choose" ) 1 else alpha, plot.it=F ), silent=quiet )
    if ( "try-error" %in% class( cv.glmnet.obj ) ) {
      tries <- 1; while( tries <= 20 && "try-error" %in% class( cv.glmnet.obj ) ) {
        cv.glmnet.obj <- try( cv.glmnet( df.tmp[ cols, ], output[ cols ], lambda=glmnet.obj$lambda, K=cv.k, trace=F,
                                        penalty.factor=in.penalties, weights=weights[ cols ],
                                        alpha=if ( alpha == "cv.choose" ) 1 else alpha, plot.it=F ), silent=quiet )
        tries <- tries + 1
      } }
    if ( "try-error" %in% class( cv.glmnet.obj ) ) return( numeric() )
    cv.glmnet.obj$alpha <- alpha

### ... choose min CV (either absolute min if cv.choose="min" or min+1se if cv.choose="min+1se")
    min.i <- which.min( cv.glmnet.obj$cv )
    min.err <- cv.glmnet.obj$cv.error[ min.i ]

    ##if ( cv.choose[ 1 ] == "min" ) best.s <- which.min( cv.glmnet.obj$cv )
    ##else
    se <- 1
    if ( grepl( "+", cv.choose[ 1 ], fixed=T ) ) {
      se <- as.numeric( gsub( "min+", "", gsub( "se", "", cv.choose[ 1 ] ) ) )
      best.s <- min( which( cv.glmnet.obj$cv <= min( cv.glmnet.obj$cv ) + se * min.err ) )
    } else best.s <- which.min( cv.glmnet.obj$cv )

    coeffs.tmp <- as.matrix( coef( glmnet.obj, s=glmnet.obj$lambda[ best.s ] ) )
    ##coeffs <- coeffs[ abs( coeffs ) >= coef.cutoff ]
    coeffs <- coeffs.tmp[ coeffs.tmp != 0, ]
    names( coeffs ) <- rownames( coeffs.tmp )[ coeffs.tmp != 0 ]
    orig.coeffs <- coeffs <- coeffs[ names( coeffs ) != "(Intercept)" ] ## Remove intercept - is this correct?
    sorted <- names( sort( abs( coeffs ), decreasing=T ) )
    coeffs <- coeffs[ abs( coeffs ) >= coef.cutoff ]
    if ( ! is.na( min.coeffs ) && length( coeffs ) < min.coeffs && length( orig.coeffs ) >= min.coeffs )
      coeffs <- orig.coeffs[ sorted[ 1:min.coeffs ] ]
    if ( ! is.na( max.coeffs ) && length( coeffs ) > max.coeffs ) coeffs <- orig.coeffs[ sorted[ 1:max.coeffs ] ]
    if ( ! quiet ) cat( boot, cv.choose[ 1 ], min.i, min.err, best.s, cv.glmnet.obj$cv[ best.s ],
                       glmnet.obj$lambda[ best.s ], length( coeffs ), "\n" )

    ## NOTE: rescaling coeffs needs to be done using glmnet, NOT lm b/c of correlated predictors
    if ( rescale.coeffs && length( coeffs ) > 0 ) { ## plug coeffs from LARS back into LM, coeffs will be larger.
      ins <- df.tmp[ ,names( coeffs ), drop=F ]
      ##inz<<-ins[cols,];outz<<-output[cols];alphaz<<-alpha;penaltyz<<-in.penalties;weightz<<-weights[cols]
      glmnet.obj2 <- my.glmnet( ins[ cols, ,drop=F ], output[ cols ], alpha=alpha, penalty.factor=in.penalties,
                               weights=weights[ cols ], ... ) ##lars( predictors, y, type="lasso", trace=F )
      ##tmp.out<<-glmnet.obj2
      coeffs.s <- t( as.matrix( coef( glmnet.obj2 ) ) )
      coeffs <- coeffs.s[ nrow( coeffs.s ), ]
      coeffs <- coeffs[ abs( coeffs ) >= coef.cutoff ]
      coeffs <- coeffs[ names( coeffs ) != "(Intercept)" ] ## Remove intercept - is this correct?
    }
    ##out.coe[[ boot ]] <- coeffs
    if ( boot == 1 ) { out <- list( coeffs=coeffs, lars.obj=glmnet.obj, cv.lars.obj=cv.glmnet.obj, best.s=best.s,
                           se=se, min.err=min.err ); return( out ) }
    coeffs
  } )

  lars.obj <- out.coe[[ 1 ]]$lars.obj
  cv.lars.obj <- out.coe[[ 1 ]]$cv.lars.obj
  best.s <- out.coe[[ 1 ]]$best.s
  se <- out.coe[[ 1 ]]$se
  min.err <- out.coe[[ 1 ]]$min.err
  out.coe[[ 1 ]] <- out.coe[[ 1 ]]$coeffs
  
  coeffs <- out.coe[[ 1 ]]
  coeffs <- coeffs[ order( abs( coeffs ), decreasing=T ) ]

  coef.quantiles <- NULL
  if ( n.boot > 1 ) {
    tmp <- unlist( out.coe )
    tmp2 <- sort( table( names( tmp ) ), decreasing=T )
    coef.quantiles <- t( sapply( names( tmp2 ), function( i ) {
      tmp3 <- tmp[ names( tmp ) == i ]
      tmp3 <- c( tmp3, rep( 0, n.boot - length( tmp3 ) ) )
      c( n=sum( names( tmp ) == i ) / n.boot, quantile( abs( tmp3 ), prob=c( 0.01, 0.05, 0.1, 0.5, 0.9, 0.95 ) ) *
        sign( mean( tmp3[ tmp3 != 0 ], na.rm=T ) ) ) } ) )
    coef.quantiles <- coef.quantiles[ ! apply( coef.quantiles, 1, function( i ) all( i[ -1 ] == 0 ) ), ]
    if ( ! quiet ) print( coef.quantiles, digits=3 )
  }
  
  ##if ( ! quiet ) { cat( "NONZERO COEFFS:\n" ); print( coeffs ) }
  ##if ( plot )
  return( list( coeffs=coeffs, coeffs.boot=out.coe, coef.quantiles=coef.quantiles,
               lars.obj=lars.obj, cv.lars.obj=cv.lars.obj, cv.choose=cv.choose, best.s=best.s,
               se=se, min.err=min.err, all.inputs=rownames( df.tmp ) ) )
  ##return( list( coeffs=coeffs, coeffs.boot=out.coe, coef.quantiles=coef.quantiles ) )
} ## end of inferelator function
