###########################################
## David J Reiss
###########################################

source( "inferelator.R" )
source( "predictelator.R" )
source( "write.inf.network.R" )
require( lars ) ## So it doesnt get loaded by each core separately
require( glmnet ) ## So it doesnt get loaded by each core separately
require( multicore )

## Good defaults: cv.choose="min+4se", aic.filter=15
## NOTE: aic.filter is just for speedup; seems to give same results   if use aic.filter=100 !

get.apply.func <- function( plot=F ) if ( multicore:::isChild() || plot || ( exists( "DEBUG" ) && DEBUG ) )
  lapply else mclapply

load.egrin.data <- function( path=".", ... ) {
  load( paste( path, "data_orig_EGRIN/egrin_newcode_workspace.RData", sep="/" ) ) ## env.map.egrin is f*cked up for some reason;
  load( paste( path, "data_orig_EGRIN/env_map_egrin.RData", sep="/" ) ) ## load the one from the original run and make the names match
  ##colnames( env.map ) <- gsub( "-", ".", colnames( env.map ), fixed=T )
  relevant.env <- c("oxygen", "illumination", "Fe", "Cu", "Co", "Mn", "Zn", "Ni", "gamma", "uv") ## from Rich code
  env.map <- env.map[ relevant.env, ]
  load( paste( path, "data_orig_EGRIN/col_map_egrin.RData", sep="/" ) ) ## col.map.egrin is also f*cked - so lets load the orig. one and transform it
  ##names( colMap ) <- gsub( "-", ".", names( colMap ), fixed=T ); col.map <- data.frame();
  col.map <- NULL
  for ( i in 1:( length( colMap ) - 1 ) ) {
    ##colMap[[ i ]]$condName <- gsub( "-", ".", colMap[[ i ]]$condName, fixed=T )
    ##colMap[[ i ]]$prevCol <- gsub( "-", ".", colMap[[ i ]]$prevCol, fixed=T )
    col.map <- rbind( col.map, as.data.frame( colMap[[ i ]] ) )
  }
  rownames( col.map ) <- names( colMap )[ 1:( length( colMap ) - 1 ) ]
  colnames( col.map )[ colnames( col.map ) == "del.t" ] <- "delta.t"
  pc <- as.character( col.map$prevCol )
  pc[ is.na( pc ) ] <- as.character( col.map$condName[ is.na( pc ) ] )
  col.map$prevCol <- as.factor( pc )
  col.map$delta.t[ is.na( col.map$delta.t ) ] <- 9999
  predictors <- c( readLines( paste( path, "data/halo/halo_tfs.txt", sep="/" ) ), rownames( env.map ) )
  data <- rbind( ratios.egrin, env.map )
  load( paste( path, "data_orig_EGRIN/egrin_coeffs.RData", sep="/" ) ) ## coeff.inf
  invisible( list( col.map=col.map, env.map=env.map, predictors=predictors, data=data,
                  clusterStack.egrin=clusterStack.egrin ) )
}


## runnit.newCM.egrin.data <- function( f, ks="all", tau=10, plot=T, coeffs=NULL, tf.groups=72, n.boot=1,
##                               boot.opt=c("resample.lars","resample.rows","resample","lars")[1], ... ) {
##   if ( is.character( f ) && file.exists( f ) ) {
##     load( f, envir=.GlobalEnv )
##     print( f )
##   } else if ( is.environment( f ) ) {
##     e <- f; rm( f )
##   }
##   ratios <- e$get.cluster.matrix()
##   ##attach( e ) ## e is environment output by cmonkey() as of version 4.3.1
##   if ( ks[ 1 ] == "all" ) ks <- 1:e$k.clust
##   if ( "egrin.data" %in% searchpaths() ) detach( egrin.data )
##   egrin.data <- load.egrin.data( ... )
##   egrin.data$data <- egrin.data$clusterStack.egrin <- NULL
##   attach( egrin.data )
##   data <- rbind( ratios, env.map )
##   out <- runnit( ks, data, col.map, predictors, e$clusterStack, tau=tau, plot=plot, coeffs=coeffs,
##          tf.groups=tf.groups, n.boot=n.boot, boot.opt=boot.opt, ... )
##   detach( egrin.data )
##   ##if ( ! is.null( f ) ) detach( e ) ##cm.detach()
##   invisible( out )
## }

### OK -- looks like good parameters to use are
## Updated (version 0.0.9, Aug. 2011)
## aic.filter <- Inf
## alph <- 0.8
## tau <- 10
## tf.groups <- Inf
## r.cutoff <- Inf
## r.filter <- Inf
## weighted <- TRUE
## cv.choose <- "min+2se"
## runnit.wrapper.halo("~/scratch/biclust/EGRIN2/EGRIN1_orig_clusters.RData",cv.choose="min+4se",tf.groups=999,alpha=0.8,tau=10,r.cutoff=2,weighted=T,aic.filter=15,plot=F)
runnit.wrapper.halo <- function( f, ks="all", ... ) {
  if ( is.character( f ) && file.exists( f ) && ( ! exists( "e" ) || e$tmp.file != f ) ) {
    load( f, envir=.GlobalEnv )
    print( f )
    assign( "tmp.file", f, env=e )
  } else if ( is.environment( f ) ) {
    e <- f; rm( f )
  }
  if ( ! exists( "ratios" ) ) ratios <- e$get.cluster.matrix() ## e is environment output by cmonkey() as of version 4.3.1
  ##if ( nrow( ratios ) == 0 ) ratios <- e$ratios ## HACK for "small-ified" env where get.cluster.matrix() doesnt work
  ##colnames( ratios ) <- gsub( "-", ".", colnames( ratios ), fixed=T )
  if ( ks[ 1 ] == "all" ) ks <- 1:e$k.clust
  if ( ! exists( "envMap" ) ) envMap <- NULL
  if ( ! exists( "colMap" ) ) colMap <- NULL
  if ( ! exists( "predictors" ) ) predictors <- readLines( "data/halo/halo_tfs.txt" )
  ## Remove variables from envMap that are not changing (or are nearly all zeroes or NAs)
  data <- ratios
  if ( ! is.null( envMap ) ) {
    envMap <- envMap[ , ! is.na( apply( envMap, 2, var, use="pair" ) ) & apply( envMap, 2, var, use="pair" ) > 0.01,
                     drop=F ]
    envMap <- envMap[ rownames( envMap ) %in% colnames( ratios ),, drop=F ]
    ratios <- ratios[ ,colnames( ratios ) %in% rownames( envMap ), drop=F ]
    data <- rbind( ratios, t( as.matrix( envMap ) ) )
    predictors <- c( predictors, colnames( envMap ) )
  }
  if ( ! is.null( colMap ) ) {
    ratios <- ratios[ ,colnames( ratios ) %in% rownames( colMap ), drop=F ]
  }
  if ( ! is.null( predictors ) ) predictors <- predictors[ predictors %in% rownames( data ) ]
  ##tmpz<<-list(data=data,colMap=colMap,predictors=predictors,envMap=envMap)
  ## Gene prefix is used to discriminate genetic predictors from env. predictors.
  out <- runnit( ks, data, colMap, predictors, e$clusterStack, gene.prefix=e$genome.info$gene.prefix, ... ) ##tau=tau, plot=plot, coeffs=coeffs,
  ##tf.groups=tf.groups, n.boot=n.boot, boot.opt=boot.opt, ... )
  invisible( out )
}

runnit <- function( ks, data, col.map, predictors, clusterStack, tau=10, plot=T, coeffs=NULL, tf.groups=Inf, n.boot=1,
                   boot.opt=c("resample.lars","resample.rows","resample","lars")[1], ... ) {
  ## Bootstrap options: "resample" -- resample cluster rows AND cols; "resample.rows" -- just resample cluster rows;
  ##   "resample.lars" -- run lars/cv.lars multiple times on resampled input matrices;
  ##   "lars" -- don't resample anything but just re-run lars/cv.lars multiple times on data (let the cv-ing in
  ##                 cv.lars be the thing that's sampled)
  in.args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments
               sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## nifty trick, eh?
  
  data <- mean.variance.normalize( data, filter=0.04 ) ## Removes most env. factors with little change
  ##rownames( data ) <- gsub( ".", "_", rownames( data ), fixed=T )
  ##predictors <- gsub( ".", "_", predictors, fixed=T )
  predictors <- predictors[ predictors %in% rownames( data ) ]
  ##colnames( data ) <- gsub( "-", ".", colnames( data ), fixed=T ) ## Halo-specific? Let's hope not!

  ## Note if tf.groups > length(tfs) OR tf.groups is 0 or NA, don't do preclustering
  if ( ! exists( "predictor.mats" ) ||
      ( ( is.na( tf.groups ) || tf.groups == 0 || tf.groups >= length( predictors ) ) &&
       length( predictor.mats$tf.groups ) != length( predictors ) ) ||
      ( ! is.na( tf.groups ) && tf.groups != 0 && tf.groups < length( predictors ) &&
       length( predictor.mats$tf.groups ) != tf.groups ) ) {
    predictor.mats <<- get.predictor.matrices( predictors, data, preclust.k=tf.groups, ... )
  }

  n.boot.lars <- 1; boot.opt.lars <- "resample"
  if ( n.boot > 1 && boot.opt %in% c( "resample.lars", "lars" ) ) {
    n.boot.lars <- n.boot; n.boot <- 1
    if ( boot.opt == "lars" ) boot.opt.lars <- "cv"
  }
  
  ##out <- list()
  ##for ( i in ks ) {
  apply.func <- get.apply.func( plot )
  if ( n.boot > 1 ) apply.func <- lapply
  out <- apply.func( ks, function( i ) {
    cluster <- clusterStack[[ i ]]
    k <- cluster$k
    ##for ( boot in 1:n.boot ) {

    apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
    if ( n.boot == 1 ) apply.func <- lapply
    out.k <- apply.func( 1:n.boot, function( boot ) {
      cat( "***      BICLUSTER:", k, boot, "\n" )
      clust <- cluster
      if ( boot > 1 ) {
        if ( boot.opt %in% c( "resample", "resample.rows" ) ) clust$rows <- sample( clust$rows, replace=T )
        if ( boot.opt == "resample" ) clust$cols <- sample( colnames( data ), length( clust$cols ), replace=F )
      }

      coeffs <- inferelate.one.cluster( clust, predictors, data, predictor.mats=predictor.mats, tau=tau,
                                       col.map=col.map, n.boot=n.boot.lars, boot.opt=boot.opt.lars, ##plot=plot, 
                                       quiet=n.boot>1, ... )

      clust.conds <- sort( coeffs$cluster.conds )

      observed <- apply( data[ clust$rows, ,drop=F ], 2, mean, na.rm=T )

      apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
      pred.ss <- do.call( rbind, apply.func( coeffs$coeffs.boot, function( b ) predictelate( clust$rows, b, data,
                                                           predictor.mats=predictor.mats, tau=tau, ... ) ) )
      pred.ts <- do.call( rbind, apply.func( coeffs$coeffs.boot, function( b ) predictelate( clust$rows, b, data,
                                               predictor.mats=predictor.mats, tau=tau, col.map=col.map, ... ) ) )

      if ( is.null( pred.ss ) ) pred.ss <- t( observed * 0 )
      if ( is.null( pred.ts ) ) pred.ts <- t( observed * 0 )

##       pred.ss <- pred.ts <- NULL
##       for ( b in 1:length( coeffs$coeffs.boot ) ) { ## If bootstrapping was done at inferelator level
##         pred.ss <- rbind( pred.ss, predictelate( clust$rows, coeffs$coeffs.boot[[ b ]], data,
##                                                 predictor.mats=predictor.mats, tau=tau, ... ) )
##         pred.ts <- rbind( pred.ts, predictelate( clust$rows, coeffs$coeffs.boot[[ b ]], data,
##                                                 predictor.mats=predictor.mats, tau=tau, col.map=col.map, ... ) )
##       }

      ##stop("NEED TO COMPUTE WEIGHTED RMSD IF WEIGHTED=T")
      if ( "weighted" %in% names( list( ... ) ) && list( ... )$weighted == TRUE ) {
        vars <- apply( data[ clust$rows, ,drop=F ], 2, var, na.rm=T )
        vars <- vars / ( abs( observed ) + 0.05 )
        vars[ is.na( vars ) | vars == 0 ] <- 1
        weights <- 1 / vars
        weights <- weights / sum( weights ) * length( weights )
      } else {
        weights <- rep( 1, ncol( data ) ); names( weights ) <- colnames( data )
      }

      rmsd.ss <- sqrt( weighted.mean( ( pred.ss[ nrow( pred.ss ), ] - observed )[ clust.conds ]^2,
                                     weights[ clust.conds ], na.rm=T ) )
      rmsd.ts <- sqrt( weighted.mean( ( pred.ts[ nrow( pred.ts ), ] - observed )[ clust.conds ]^2,
                                     weights[ clust.conds ], na.rm=T ) )
      not.clust.conds <- colnames( data )[ ! colnames( data ) %in% clust.conds ]
      rmsd.ts.out <- sqrt( weighted.mean( ( pred.ts[ nrow( pred.ts ), ] - observed )[ not.clust.conds ]^2,
                                         weights[ not.clust.conds ], na.rm=T ) )
      
      ##if ( plot ) {
      coeffs$plot.info$main <- paste( "Bicluster", cluster$k, cluster$nrows, "genes" )    
      coeffs$plot.info$clust.conds.plot <- c( clust.conds,
                                             sort( colnames( data )[ ! colnames( data ) %in% clust.conds ] ) )
      coeffs$plot.info$n.conds <- length( clust.conds )
      ##}
        
      if ( n.boot <= 1 ) cat( k, tau, rmsd.ss, rmsd.ts, rmsd.ts.out, "\n" )
      coeffs$pred.ss <- pred.ss ##; coeffs$rmsd.ss <- rmsd.ss
      coeffs$pred.ts <- pred.ts ##; coeffs$rmsd.ts <- rmsd.ts; coeffs$rmsd.ts.out <- rmsd.ts.out
      coeffs$rmsd <- c( ss=rmsd.ss, ts=rmsd.ts, ts.out=rmsd.ts.out )
      coeffs$observed <- observed
      coeffs$n.boot <- n.boot
      coeffs$boot.opt <- boot.opt
      
      attr( coeffs, 'class' ) <- 'coeff.obj'
      if ( boot > 1 ) coeffs$plot.info <- NULL
      coeffs
    } )
    names( out.k ) <- paste( k, 1:n.boot, sep="." )

    if ( n.boot > 1 ) {
      cc.tmp <- out.k ##coeffs[ grep( paste( "^", k, sep="" ), names( out ) ) ]
      nb <- max( n.boot, n.boot.lars )
      cc.tmp <- cc.tmp[ sapply( cc.tmp, length ) > 0 ]
      cc <- lapply( cc.tmp, "[[", "coeffs" )
      tmp <- cc; names( tmp ) <- NULL; tmp <- unlist( tmp )
      tmp2 <- sort( table( names( tmp ) ), decreasing=T )
      coef.quantiles <- t( sapply( names( tmp2 ), function( i ) {
        tmp3 <- tmp[ names( tmp ) == i ]
        tmp3 <- c( tmp3, rep( 0, nb - length( tmp3 ) ) )
        c( n=sum( names( tmp3 ) == i ) / nb, quantile( abs( tmp3 ), prob=c( 0.01, 0.05, 0.1, 0.5, 0.9, 0.95 ) ) *
          sign( mean( tmp3[ tmp3 != 0 ], na.rm=T ) ) ) } ) )
      coef.quantiles <- coef.quantiles[ ! apply( coef.quantiles, 1, function( i ) all( i[ -1 ] == 0 ) ), ]
      ##n.tot <- sapply( rownames( coef.quantiles ), function( n ) mean( sapply( sapply( cc.tmp, "[[", "coeffs" ),
      ##                                                     function( i ) sum( grepl( n ,names( i ) ) ) ) >= 1 ) )
      ##coef.quantiles <- cbind( n.tot, coef.quantiles )
      ##print( coef.quantiles, digits=3 )
      out.k[[ 1 ]]$coef.quantiles <- coef.quantiles
    } else if ( n.boot.lars > 1 ) {
      print( out.k[[ 1 ]]$coef.quantiles, digits=3 )
    }
    
    if ( plot ) {
      ##if ( n.boot == 1 ) try( plot.coeff.obj( out[[ k ]], ... ) )
      ## only plot the boot results for this k ...
      ##!else try( plot.coeff.obj( out[ grep( paste( "^", k, sep="" ), names( out ) ) ], ... ) )
      try( plot.coeff.obj( out.k, ... ) )
    }
    attr( out.k, 'class' ) <- 'coeff.obj'
    ##cat(k,class(out.k)," ");print(names(out.k[[1]]))
    ##if(class(out.k)=="character"){out.k<-list(out.k);cat("CHARACTER            ",k,"\n");save(k,out.k,file="qqqz")}
    out.k
    ##out <- c( out, out.k )
  } )
  out <- do.call( c, out )
  attr( out, "CALL" ) <- match.call( expand.dots=T )
  attr( out, 'class' ) <- 'coeff.obj'
  invisible( out )
}


