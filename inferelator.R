###################################################
## Inferelator version B, fewer false negatives
## Phu T Van and David J Reiss, ISB
###################################################

###########################################
## inferelator
###########################################

## 1. idea - instead of pre-clustering, use profiles of biclusters that have at least 1 TF in them
## 2. (DONE) only inferelate on the conditions IN THE BICLUSTER !!!
## 3. better pre-filtering of predictors for each bicluster - use correlation of TF prof. with dy/dt of bicluster
## 4. use exponentially weighted sum of dy/dt's from earlier t's than just t_{i-1}
## 5. (ONLY NEEDED FOR ORIG. EGRIN BICLUSTERS) add "loner" genes to clusterStack as single-gene clusters (esp. if any TFs are loners!)
## 6. (DONE - but has weird result - no combo predictors are included) use equal #s of neg and pos influences (determined by "cor") as potential predictors for each bicluster
## 7. (DONE) use glmnet instead of lars. (Note default args. for glmnet is SAME as lars.)
## 8. Now that we are using glmnet, we can (a) DONE add weights to each obs (based on variance in biclust profile),
##       (b) DONE upweight predictors using "penalty.factor" glmnet arg, and (c) do CV to choose alpha as well as
##       lambda

require( glmnet ) ## Just in case
require( lars ) ## Just in case
source( "inferelator_enet.R" )

combine.symbol <- "~~" ## "!!" ## "."

inferelate.one.cluster <- function( cluster, predictors, data, col.map=NULL, conds.use=c("clust","ALL")[1], 
                                   quiet=F, plot=T, shrink.opt=c("glmnet","lars")[1], predictor.mats=NULL,
                                   weighted=T, weights=NA, ... ) {

  if ( ! exists( "predictor.mats" ) || is.null( predictor.mats ) )
    predictor.mats <<- get.predictor.matrices( predictors, data, quiet=quiet, ... )
  
  cluster.rows <- cluster$rows
  cluster.conds <- colnames( data )
  ##cluster$cols <- gsub( "-", ".", cluster$cols )
  if ( conds.use == "clust" ) {
    if ( exists( "col.map" ) && ! is.null( col.map ) ) cluster.conds <- fill.in.time.series( cluster$cols, col.map,
                                         fill.all.ts.frac=1.25, remove.all.ts.frac=0.1, fill.ts.gap.size=1 )
    else cluster.conds <- cluster$cols
  }
  cluster.rows <- cluster.rows[ cluster.rows %in% rownames( data ) ]
  cluster.conds <- cluster.conds[ cluster.conds %in% colnames( data ) ]
  if ( is.null( cluster.rows ) || ( is.null( cluster.conds ) && conds.use == 'clust' ) ) {
    out <- list( k=cluster$k, coeffs=numeric(), possibly.regulates=NULL, cluster.conds=cluster.conds,
                 coeffs.boot=NULL, coef.quantiles=NULL, all.inputs=NULL, plot.info=NULL ) ##, params=in.args ) )
    return( out )
  }

  ##cluster.conds <- gsub( "-", ".", cluster.conds, fixed=T ) ## Halo-specific? Let's hope not!
  cluster.profile <- apply( data[ cluster.rows, ,drop=F ], 2, mean, na.rm=T )
  if ( any( is.na( cluster.profile ) ) )
    cluster.conds <- cluster.conds[ -which( is.na( cluster.profile ) ) ] ## If all are NA, ignore!
  cluster.weights <- NA
  if ( weighted ) {
    cluster.vars <- apply( data[ cluster.rows, ,drop=F ], 2, var, na.rm=T )
    cluster.vars <- cluster.vars / ( abs( cluster.profile ) + 0.05 )
    cluster.vars[ is.na( cluster.vars ) | cluster.vars == 0 ] <- 1
    cluster.weights <- 1 / cluster.vars
    cluster.weights[ is.na( cluster.weights ) | is.infinite( cluster.weights ) ] <- 1
    ##cluster.weights <- cluster.weights / sum( cluster.weights ) * length( cluster.weights )
  }

  if ( ! is.na( weights ) ) {
    if ( ! weighted ) cluster.weights <- weights
    else cluster.weights <- cluster.weights * weights
  }
  
  ## remove TF predictors (tfgroups) if any TF in the tfgroup is in the bicluster
  ## also remove TF predictors (tfgroups) with high ( >0.9) correlation with bicluster
  ## but make sure to add them back as "possiblyRegulates" predictors
  tmp <- get.cluster.predictors( cluster.rows, cluster.profile[ cluster.conds ],
                                predictor.mats$predictor.mat[ ,cluster.conds ],
                                predictor.mats$predictor.mat.ands[ ,cluster.conds ],
                                predictor.mats$tf.groups, predictor.mats$env.names, quiet=quiet, ... )
  possibly.regulates <- tmp$possibly.regulates
  predictor.mat <- rbind( predictor.mats$predictor.mat, predictor.mats$predictor.mat.ands )
  predictor.mat <- predictor.mat[ rownames( predictor.mat ) %in% tmp$predictors, ,drop=F ]
  predictor.mat <- mean.variance.normalize( predictor.mat, filter=NA )
  
  if ( ! quiet ) cat( "Inferelating on biclust #", cluster$k, "using", length( cluster.conds ), "conditions and",
                     nrow( predictor.mat ), "predictors:\n", paste( rownames( predictor.mat ), sep=", " ), "\n")

  ## this line actually runs Inferelator on the profile
  if ( shrink.opt == "glmnet" ) {
    coeffs <- inferelator.enet( cluster.profile, predictor.mat, cluster.conds, col.map=col.map, 
                               quiet=quiet, weights=cluster.weights, ... ) ##logit.scale=logit.scale, plot=T, 
  } else if ( shrink.opt == "lars" ) {
    coeffs <- inferelator( cluster.profile, predictor.mat, cluster.conds, col.map=col.map, 
                               quiet=quiet, ... ) ##logit.scale=logit.scale,plot=T, 
  }
  if ( ! quiet ) { cat( "BICLUSTER", cluster$k, ": NONZERO COEFFS:\n" ); print( coeffs$coeffs ) }
  singleresult <- coeffs$coeffs
  coeffs.boot <- coeffs$coeffs.boot
  coef.quantiles <- coeffs$coef.quantiles
  all.inputs <- coeffs$all.inputs
  coeffs$coeffs <- coeffs$coeffs.boot <- coeffs$coef.quantiles <- coeffs$all.inputs <- NULL
  
  ##if ( plot ) {
  ##conds <- c( cluster.conds, colnames( data )[ ! colnames( data ) %in% cluster.conds ] )
    ##coeffs$cluster.profile <- cluster.profile[ conds ]
  coeffs$main <- paste( "Bicluster", cluster$k, cluster$nrows, "genes" )
  ##}

  ##in.args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments
  ##             sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## nifty trick, eh?
  ##in.args$data <- sprintf( "Matrix of dimensions %d x %d", nrow( data ), ncol( data ) )
  ##in.args$predictor.mats <- NULL
  
  ## if there are non-zero coefficients, add result for single cluster to the bigger resultset
  if ( length( singleresult ) > 0 ) {
    ##if ( plot ) {
    conds <- c( cluster.conds, colnames( data )[ ! colnames( data ) %in% cluster.conds ] )
    coeffs$predictor.mat <- predictor.mat[ names( singleresult ), conds, drop=F ]
    coeffs$colors <- c( "red", ifelse( singleresult > 0, "#ffaaaa", "#aaffaa" ) )
    ##}
    return( list( k=cluster$k, coeffs=singleresult, possibly.regulates=possibly.regulates,
                 cluster.conds=cluster.conds, coeffs.boot=coeffs.boot, coef.quantiles=coef.quantiles,
                 all.inputs=all.inputs, plot.info=coeffs ) ) ##, params=in.args ) )
  } else {
    ## if there are null coefficients (null = didn't meet cut-off), still return a cluster number so we know it was processed ... also there can be non-null bootstrapped coeffs
    return( list( k=cluster$k, coeffs=numeric(), possibly.regulates=possibly.regulates, cluster.conds=cluster.conds,
                 coeffs.boot=coeffs.boot, coef.quantiles=coef.quantiles, all.inputs=all.inputs, plot.info=coeffs ) ) ##, params=in.args ) )
  }
} ## end of inferelate function

get.predictor.matrices <- function( predictors, data, gene.prefix="DVU", ##predictor.mat=NULL, 
                                   preclust.k=NA, funcs=NULL, quiet=F, ... ) { ##"min"
  if ( ! quiet ) cat( "Computing predictor matrices...\n" )

  ## make sure all the TFs given exist in the ratios given
  predictors <- predictors[ predictors %in% rownames( data ) ]
  predictors.genetic <- grep( paste( "^", gene.prefix, sep="" ), predictors, ignore.case=T, value=T )
  predictors.genetic <- predictors.genetic[ ! grepl( 'knockout', predictors.genetic, ignore.case=T ) ]
  env.names <- setdiff( predictors, predictors.genetic )

  ## Note if preclust.k > length(predictors) OR tf.groups is 0 or NA, don't do preclustering
  if ( ! is.na( preclust.k ) && preclust.k != 0 && preclust.k < length( predictors ) ) {
    ## "precluster" ONLY THE GENETIC predictors into TFGROUPs (don't want to mix environmental factors)
    ## but append the environmental ones so we can have env+TFGROUP combos later
    tmp <- preclust.tfs.kmeans( data=data, tfs=predictors.genetic, clust.count=preclust.k, ... )
    
    predictor.mat <- tmp$result
    tf.groups <- tmp$tf.groups ## Need to make sure no TFgroups are included as predictors of biclusters where there is a TFgroup member in the bicluster
    rm( tmp )
  } else {
    predictor.mat <- data[ predictors, ] ##tf.matrix
    rownames( predictor.mat ) <- predictors ##paste( "TFGROUP", 1:nrow( predictor.mat ), sep="" )
    tf.groups <- as.list( predictors ); names( tf.groups ) <- rownames( predictor.mat )
  }

  ## add the real -environmental- predictors to predictor matrix...
  if ( length( env.names ) > 0 && any( ! env.names %in% rownames( predictor.mat ) ) )
    predictor.mat <- rbind( predictor.mat, data[ env.names[ ! env.names %in% rownames( predictor.mat ) ], ] )

  tmp <- unique( predictor.mat )
  if ( any( ! rownames( tmp ) %in% rownames( predictor.mat ) ) )
    cat( "Predictor", rownames( tmp )[ ! rownames( tmp ) %in% rownames( predictor.mat ) ], "is not unique. Removing.\n" )
  predictor.mat <- tmp; rm( tmp )
  
  predictor.mat.ands <- NULL
  if ( ! is.na( funcs ) && length( funcs ) > 0 ) {
    if ( ! quiet ) cat( "Computing combined predictor matrix...\n" )
    ## make combinatory predictors from TFGROUPs 
    predictor.mat.ands <- make.combined.predictors( predictor.mat, ##predictors=rownames( predictor.mat ),
                                                   funcs=funcs, ... )
  }

  list( predictor.mat=predictor.mat, predictor.mat.ands=unique( predictor.mat.ands ),
       genetic.names=unique( predictors.genetic ), tf.groups=tf.groups, env.names=unique( env.names ) )
}

## old good: aic.filter.cutoff=5; now using 15. NOTE: aic.filter is just for speedup; seems to give same results
##   if use aic.filter=100 !
get.cluster.predictors <- function( cluster.rows, cluster.profile, predictor.mat, predictor.mat.ands,
                                   tf.groups, env.names, r.cutoff=NA, aic.filter.cutoff=NA, quiet=F,
                                   force.include=env.names, ... ) {
  ## remove TF predictors (tfgroups) if any TF in the tfgroup is in the bicluster
  ## also remove TF predictors (tfgroups) with high ( >r.cutoff) correlation with bicluster
  ## but make sure to add them back as "possiblyRegulates" predictors
  ## Can force certain regulators to be included; default is to include all env factors (useful if we are up-weighting certain TFs via penalty in glmnet -- TODO need to make this work for combination factors.

  tf.groups <- tf.groups[ names( tf.groups ) %in% rownames( predictor.mat ) ] ## May have gotten filtered
  force.include <- force.include[ force.include %in% names( tf.groups ) ]
  predictor.mat.env <- predictor.mat[ force.include, ,drop=F ]

  possibly.regulates <- numeric() ## make it a named vector with the value of the correlation or # in. ##character()
  is.in <- sapply( tf.groups, function( j ) mean( j %in% cluster.rows, na.rm=T ) )
  high.cor <- apply( predictor.mat[ names( tf.groups ), ], 1, cor, cluster.profile, use="pairwise" )
  high.cor[ is.na( high.cor ) ] <- 0
  ##tmpz<<-list(predictor.mat=predictor.mat,cluster.profile=cluster.profile,is.in=is.in,high.cor=high.cor);stop()

  ##possibly.regulates <- NULL
  ## Hack to get rid of predictors (e.g. "gamma") that dont change over bicluster's conditions and thus high.cor is NA
  if ( any( is.na( high.cor ) ) ) {
    tmp <- which( is.na( high.cor ) )
    ##possibly.regulates <- high.cor[ tmp ]
    ##predictor.mat <- predictor.mat[ -tmp, ]
    is.in[ tmp ] <- 1 ## <- is.in[ -tmp ]
    high.cor[ tmp ] <- 1 ## <- high.cor[ -tmp ]
  }

  ## Get rid of predictors that are either in the cluster or highly (positively) correlated with cluster -
  ##    but add them as "possibly.regulates" factors
  if ( any( is.in > 0 | ( ! is.na( r.cutoff ) && high.cor > r.cutoff ) ) ) {
    possibly.regulates <- c( is.in[ is.in > 0 ], high.cor[ is.in <= 0 & high.cor > r.cutoff ] )
    is.in <- which( is.in > 0 | high.cor > r.cutoff )
    ##possibly.regulates <- names( is.in )
    predictor.mat <- predictor.mat[ -is.in, ]
    is.in <- unique( sort( unlist( lapply( names( is.in ),
                                          function( j ) grep( j, rownames( predictor.mat.ands ) ) ) ) ) )
    predictor.mat.ands <- predictor.mat.ands[ -is.in, ]
  }
  if ( ! quiet && length( possibly.regulates ) > 0 ) cat( "POSSIBLY REGULATES (but removed from regression):",
                                               paste( unique( names( possibly.regulates ) ), sep=", " ), "\n" )
  ##possibly.regulates <- possibly.regulates[ ! is.na( possibly.regulates ) ]

  if ( ! is.na( aic.filter.cutoff ) && ! is.infinite( aic.filter.cutoff ) && aic.filter.cutoff != 0 ) {
    ## filter single TFGROUPs and combinatory TFGROUPs -separately- by AIC, keeping the best in both cases
    ## Make sure to keep all singleton env. predictors though!
    best.singletons <- filter.by.aic( cluster.profile, predictor.matrix=predictor.mat,
                                     top.aic.to.keep=aic.filter.cutoff, ... )
    ##best.singletons <- predictor.mat
    best.combos <- NULL
    if ( ! is.null( predictor.mat.ands ) ) {
      best.combos <- filter.by.aic( cluster.profile, predictor.matrix=predictor.mat.ands,
                                   top.aic.to.keep=nrow( best.singletons ) * 20, ... ) ##filter.by.aic * 20)
      ## Filter out combos that have a high correlation with any of the singletons (over these conditions)
      if ( ! is.na( r.cutoff ) && r.cutoff < 1 ) {
        tmp.cor <- t( cor( t( best.singletons ), t( best.combos ) ) )
        to.elim <- which( tmp.cor > r.cutoff - 0.2, arr=T )
        best.combos <- best.combos[ ! rownames( best.combos ) %in% rownames( to.elim ), ,drop=F ]
      }
      if ( nrow( best.combos ) > aic.filter.cutoff ) best.combos <- best.combos[ 1:aic.filter.cutoff, ]
      if ( nrow( best.combos ) > nrow( best.singletons ) + length( env.names ) )
        best.combos <- best.combos[ 1:( nrow( best.singletons ) + length( env.names ) ), ]
    }
    predictor.mat <- unique( rbind( predictor.mat[ rownames( best.singletons ), ,drop=F ],
                           predictor.mat.env,
                           predictor.mat.ands[ rownames( best.combos ), ,drop=F ] ) )
    ##predictor.mat <- rbind( best.singletons, best.combos )
  } else {
    predictor.mat <- unique( rbind( predictor.mat, predictor.mat.env, predictor.mat.ands ) )
  }

  list( possibly.regulates=possibly.regulates, predictors=rownames( predictor.mat ) ) ##.mat=predictor.mat )
}  

###########################################
## preclust.tfs
###########################################
preclust.tfs.kmeans <- function( data, tfs, clust.count, n.iter=200, n.start=25, seed=31337, r.cutoff=0.85, ... ) {
### cluster tf's w/ kmeans before feeding into LARS
### otherwise too many predictors and LARS will fail 
### manually set seed so will get same clusters each time 
  
  if ( ! is.na( seed ) ) set.seed( seed )

  tf.matrix <- data[ tfs[ tfs %in% rownames( data ) ], ]
  data.c <- kmeans( tf.matrix, clust.count, iter.max=n.iter, nstart=n.start )
  result <- data.c$centers
  
### name the kmean-clustered centers
  rownames( result ) <- paste( "TFGROUP", 1:clust.count, sep="" )

  tf.groups <- lapply( 1:length( data.c$size ), function( i ) names( which( data.c$cluster == i ) ) )
  names( tf.groups ) <- rownames( result )
  
  cat( "Preclustered with k-means, predictor matrix is", nrow( result ), "x", ncol( result ), "\n" )

  ## Get correlations between each tfgroup profile and all tfs -- see if there are any tfs outside of the
  ##    tfgroup cluster that are highly correlated with it; if so, add it.
  tmp <- apply( result, 1, function( i ) apply( data[ tfs, ], 1, cor, i ) )
  if ( any( tmp > r.cutoff ) ) {
    high.cors <- apply( tmp, 2, function( i ) which( i > r.cutoff ) )
    to.be.added <- lapply( names( tf.groups ), function( i )
                          names( high.cors[[ i ]] )[ ! names( high.cors[[ i ]] ) %in% tf.groups[[ i ]] ] )
    for ( i in 1:length( tf.groups ) ) tf.groups[[ i ]] <- unique( c( tf.groups[[ i ]], to.be.added[[ i ]] ) )
    result <- t( sapply( tf.groups, function( i ) apply( data[ i, ,drop=F ], 2, mean, na.rm=T ) ) )
  }

  ## Change names of singleton tfgroups to the name of the gene they contain
  for ( i in 1:length( tf.groups ) ) if ( length( tf.groups[[ i ]] ) == 1 )
    names( tf.groups )[ i ] <- rownames( result )[ i ] <- tf.groups[[ i ]][ 1 ]
  
  return( list( result=result, tf.groups=tf.groups ) )
} #end of preclust function

## ## New function using Rich's old method (cor.cutoff)
## preclust.tfs.orig <- function( ratios, tfs, cor.cut=0.80 ) {
##   ## join tfs
##   num.tfs <- length( tfs )
##   cat( "Starting with ", num.tfs, " tfs\n preparing to reduce to cor sets\n" )

##   still.in <- rep( TRUE, num.tfs )
##   names( still.in ) <- tfs

##   tf.merge.cors <- list()
##   tf.cor <- cor( t( ratios[ tfs, ] ) )

##   for ( i in 1:(num.tfs-1) ) {
##     if ( ! still.in[ i ] ) next
##     tf.merge.cors[[ tfs[ i ] ]] <- tfs[ i ]
##     for ( j in (i+1):num.tfs ) {
##       if ( tf.cor[ i, j ] > cor.cut ) {
##         still.in[ j ] <- FALSE
##         tf.merge.cors[[ tfs[ i ] ]] <- c( tf.merge.cors[[ tfs[ i ] ]], tfs[ j ] ) 
##       }
##     }
##   }

##   tfs.new <- tfs[ still.in ]

##   invisible( list( tfs=tfs.new, tf.cors=tf.merge.cors) )
## ##   loners <- as.list( tfs[ still.in ] )
## ##   out <- c( loners, tf.merge.cors )
## ##   names( out ) <- paste( "TFGROUP", 1:length( out ), sep="" )

## ##   centers <- t( sapply( out, function( i ) apply( ratios[ i, ,drop=F ], 2, mean, na.rm=T ) ) )
## ##   cat( "Preclustered with cor <", cor.cut, ", predictor matrix is", nrow( centers ), "x", ncol( centers ), "\n" )
## ##   rownames( centers ) <- names( out )

## ##   invisible( list( tf.groups=out, result=centers ) )
## }

###########################################
## filter.by.aic
## TODO: replace AIC as measure with mutual information using minet package!!!
## ALSO: let's make force.pos.neg=0.5 by default to include more repressors
##  NOTE with glmnet this is not the way to add repressors - should just down-penalize them instead!
###########################################
filter.by.aic <- function( mean.profile, predictor.matrix, top.aic.to.keep, ##r.min.cutoff=0.2,
                          force.pos.neg=NA, ... ) {
  ## this function takes a matrix of ratios and using AIC, keeps only the top <specified number>
  ## i.e.  ones with lowest AIC's , throwing out the rest
  ## Don't bother computing AIC for predictors with cor < r.min.cutoff (for speed)

  cors <- apply( predictor.matrix, 1, cor, mean.profile, use="pairwise" )
  
  ##paic <- rep( Inf, nrow( predictor.matrix ) )
  
  ##for ( j in which( abs( cors ) >= r.cutoff ) )
  ##  paic[ j ] <- AIC( lm( mean.profile ~ as.numeric( predictor.matrix[ j, ] ) ) )
  ## require( multicore )
  ## apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
  ## tmp <- unlist( apply.func( which( abs( cors ) >= r.min.cutoff ), function( j )
  ##                           AIC( lm( mean.profile ~ as.numeric( predictor.matrix[ j, ] ) ) ) ) )
  mi <- log( 1 - cors^2 + 1e-99 ) ## continuous mi measure (see minet package for more details)
  ##paic[ which( abs( cors ) >= r.min.cutoff ) ] <- tmp
  paic <- mi ##[ which( abs( cors ) >= r.min.cutoff ) ] <- tmp
  names( paic ) <- rownames( predictor.matrix )
  paic <- paic[ ! is.na( paic ) ]
  
  ## keep only the best <specified number> of predictors
  if ( is.na( force.pos.neg ) ) {
    best.preclusts <- names( sort( paic ) )[ 1:min( length( paic ), top.aic.to.keep ) ]
  } else { ## force.pos.neg=0.2 -> Force there to be 20% positive predictors and 80% negative (otherwise pos. predictors will dominate)
    best.preclusts <- unique(
      c( names( sort( cors, decreasing=T ) )[ 1:min( length( cors ), round( top.aic.to.keep * force.pos.neg ) ) ],
        names( sort( cors ) )[ 1:min( length( cors ), round( top.aic.to.keep * ( 1 - force.pos.neg ) ) ) ] ) )
  }
  result <- predictor.matrix[ best.preclusts, ,drop=F ]

  ##cat( "Filtered by AIC, predictor matrix is now ", nrow( result ), "x", ncol( result ), "\n" )	
  return( result )	
} #end of filter.by.aic function

###########################################
## make.combined.predictors
## r.filter -- threshold by which we remove pair-wise predictors that are highly correlated with either of the
##    two single predictors (note this is DIFFERENT from r.cutoff and r.min.cutoff. Aint that confusing?)
###########################################
make.combined.predictors <- function( predictor.mat, predictors=rownames( predictor.mat ), funcs="min",
                                     r.filter=0.7, ... ) {
  if ( is.null( funcs ) || is.na( funcs ) ) return( predictor.mat )
### Make array of min()s and max()s:
  result <- NULL
  tmp <- t( combn( predictors, 2 ) )
  tmp <- tmp[ tmp[ ,1 ] != tmp[ ,2 ], ]
  apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
  for ( func in funcs ) {
    ##tmp2 <- t( apply( tmp, 1, function( i ) apply( predictor.mat[ i, ], 2, func, na.rm=T ) ) )
    tmp2 <- do.call( rbind, apply.func( 1:nrow( tmp ), function( i ) apply( predictor.mat[ tmp[ i, ], ], 2,
                                                                           func, na.rm=T ) ) )
    tmp2[ is.infinite( tmp2 ) ] <- NA
    rownames( tmp2 ) <- paste( tmp[ ,1 ], tmp[ ,2 ], rep( func, nrow( tmp ) ), sep=combine.symbol )
    result <- rbind( result, tmp2 ); rm( tmp2 )
  }

  cat( "Combined", funcs, "predictor matrix is", nrow( result ), "x", ncol( result ), "\n" )
  out <- result

### Also Filter out all combined profiles that have a correlation > r.filt with either of the 2 indiv. profiles
  if ( ! is.na( r.filter ) && r.filter > 0 && r.filter < 1 ) {
    apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
    all.cors <- cor( t( predictor.mat ), t( result ), use="pairwise" )
    tmp <- apply.func( 1:nrow( result ), function( i ) {
      nm <- strsplit( rownames( result )[ i ], combine.symbol, fixed=T )[[ 1 ]]
      tmp.out <- NULL
      ttmp <- all.cors[ nm[ 1:2 ], i ]; ttmp[ is.na( ttmp ) ] <- 0
      if ( ! any( ttmp > r.filter ) ) tmp.out <- rownames( result )[ i ]
      tmp.out
    } )
    tmp <- do.call( "c", tmp )
    out <- result[ tmp, ] ##do.call( rbind, tmp )
    cat( "Filtered for cor <=", r.filter, ", combined predictor matrix is now ", nrow( out ), "x", ncol( out ), "\n" )
  }

  attr( out, "r.filter" ) <- r.filter
  return( out )
} #end of make.combined.predictors function

fill.in.time.series <- function( cols, col.map,
                                fill.all.ts.frac=1.25, ## If >1.25 of a ts in a clust (relative to fraction of all cols that are in the cluster), put the whole thang in there
                                remove.all.ts.frac=0.1, ## If <0.75 of a ts in a cluster (again, relative), get rid of it all.
                                fill.ts.gap.size=1, ... ) { ## Otherwise if there's a gap in a ts in the cols in a cluster, fill in the gap
  ##cols <- gsub( "-", ".", cols, fixed=T )
  out.cols <- cols
  all.cols <- rownames( col.map )
  firsts <- which( col.map$is1stLast == "f" )
  lasts <- which( col.map$is1stLast == "l" )
  ts <- lapply( 1:length( firsts ), function( i ) all.cols[ firsts[ i ]:lasts[ i ] ] )
  cols.frac <- length( cols ) / length( all.cols )
  for ( j in 1:length( ts ) ) {
    frac <- mean( ts[[ j ]] %in% cols )
    if ( frac >= cols.frac * fill.all.ts.frac ) { ## && frac < 1 ) {
      new.cols <- ts[[ j ]][ ! ( ts[[ j ]] %in% cols ) ]
      out.cols <- unique( c( out.cols, new.cols ) )
    } else if ( frac < cols.frac * remove.all.ts.frac ) {
      remove.cols <- ts[[ j ]][ ts[[ j ]] %in% cols ]
      out.cols <- out.cols[ ! ( out.cols %in% remove.cols ) ]
    } else {
      for ( q in 1:fill.ts.gap.size ) {
        is.in <- which( ts[[ j ]] %in% out.cols )
        expand <- unique( sort( c( is.in, is.in-1, is.in+1 ) ) )
        expand <- expand[ expand > 0 & expand <= length( ts[[ j ]] ) ]
        out.cols <- unique( c( out.cols, ts[[ j ]][ expand ] ) )
      }
    }
  }
  out.cols
}  

## Convert predicted variable to the difference equation, i.e. from y to tau/dt * (y_i - y_(i-1)) + y_(i-1)
get.input.matrix <- function( profile, predictor.mat, conds.use="ALL", col.map=NULL, tau=10, ratio.cutoff=3,
                             quiet=F ) {
  out.tmp <- profile
  in.tmp <- predictor.mat

  ##  tmpz<<-list(out.tmp=out.tmp,in.tmp=in.tmp,profile=profile,predictor.mat=predictor.mat,col.map=col.map);stop()
  if ( ! is.null( col.map ) && ! is.na( tau ) && tau > 0 ) {
    if ( ! quiet ) cat( "Time series data supplied, converting predictors to difference equation... \n" )
    if ( ! quiet ) cat( "Tau =", tau, "\n" )
    conds <- colnames( predictor.mat ) ##colnames( data )
    cm <- col.map[ conds, ]
    good.i <- ( ( cm$isTs == TRUE ) & ( cm$is1stLast %in% c( "m", "l" ) ) ) | ##| col.map[ conds, "" ] == "l" ) ) |
                ( cm$isTs == FALSE & cm$is1stLast == "e" )
    curs <- as.character( cm$condName[ good.i ] )
    prevs <- as.character( cm$prevCol[ good.i ] ) ##[ conds, "prevCol" ]
    prevs[ is.na( prevs ) ] <- curs[ is.na( prevs ) ]
    del.ts <- as.numeric( as.character( cm$delta.t[ good.i ] ) ) ##[ conds, "delta.t" ]
    del.ts[ del.ts < 1 ] <- 1
    tmp <- curs %in% names( out.tmp ) & prevs %in% names( out.tmp ) & prevs %in% colnames( predictor.mat )
    prevs <- prevs[ tmp ]; curs <- curs[ tmp ]; del.ts <- del.ts[ tmp ]
    out.tmp <- ( ( tau / del.ts ) * ( out.tmp[ curs ] - out.tmp[ prevs ] ) ) + out.tmp[ prevs ]
    in.tmp <- predictor.mat[ ,prevs ]
    ##tmpz<<-list(conds=conds,cm=cm,good.i=good.i,curs=curs,prevs=prevs,del.ts=del.ts,out.tmp=out.tmp,in.tmp=in.tmp)
    colnames( in.tmp ) <- names( out.tmp )
  } else {
    if ( ! quiet ) cat( "Time series data NOT supplied, using Tau = 0.\n" )
  }

  out.tmp[ out.tmp > ratio.cutoff ] <- ratio.cutoff
  out.tmp[ out.tmp < -ratio.cutoff ] <- -ratio.cutoff
  in.tmp[ in.tmp > ratio.cutoff ] <- ratio.cutoff
  in.tmp[ in.tmp < -ratio.cutoff ] <- -ratio.cutoff

  if ( conds.use[ 1 ] == "ALL" ) conds.use <- names( out.tmp )
  df.tmp <- t( in.tmp[ ,names( out.tmp ) %in% conds.use ] ) ##predictor.mat[ predictors, ] )
  df.tmp <- df.tmp[ , ! is.na( apply( df.tmp, 2, var, use="pair" ) ) &
                   apply( df.tmp, 2, var, use="pair" ) > 0.01 ]
  output <- as.numeric( out.tmp[ names( out.tmp ) %in% conds.use ] ) ##predictor.mat[ predicted[ 1 ], ] ) ### Linear output
  names( output ) <- names( out.tmp )[ names( out.tmp ) %in% conds.use ]
  df.tmp[ is.na( df.tmp ) ] <- 0 ## why! -- this is OK - just means NA's contribute NOTHING to the sum.
  list( inp=df.tmp, outp=output ) ##, profile=profile, predictor.mat=predictor.mat, col.map=col.map, conds.use=conds.use )
}
                                      
###########################################
## inferelator
###########################################
## old good: cv.choose="min+2.5se"
inferelator <- function( profile, predictor.mat, conds.use, col.map=NULL, tau=10, ##tau.best=NULL,
                        ratio.cutoff=3, coef.cutoff=0.02, cv.k=10, cv.choose="min+2se", ##logit.scale=0.25,
                        n.boot=1, boot.opt=c( "resample", "cv" ), rescale.coeffs=T, quiet=T,
                        max.coeffs=NA, min.coeffs=NA, ... ) { ##plot=T, 

  ##require( lars, quietly=T )

  ## this function is the MEAT - run lars and output coefficients	
  if ( cv.choose == "min" ) cv.choose <- "min+0se"

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
  df.tmp <- tmp$inp
  output <- tmp$outp
  rm( tmp )

  ##tmpz<<-cbind(apply(df.tmp,2,cor,output),apply(predictor.mat[colnames(df.tmp),rownames(df.tmp)],1,cor,profile[rownames(df.tmp)]))

  ##out.coe <- list() ##numeric()
  ##for ( boot in 1:n.boot ) {
  apply.func <- get.apply.func() ##if ( multicore:::isChild() ) lapply else mclapply
  out.coe <- apply.func( 1:n.boot, function( boot ) {
    cols <- 1:length( output )
    if ( boot > 1 && boot.opt == "resample" ) cols <- sample( cols, replace=T )

### Run lars on predicted ~ predictors
    lars.obj <- try( lars( df.tmp[ cols, ], output[ cols ], type="lasso", trace=F ), silent=quiet )
    if ( class( lars.obj ) == "try-error" ) {
      tries <- 1; while( tries <= 20 && class( lars.obj ) == "try-error" ) {
        lars.obj <- try( lars( df.tmp[ cols, ], output[ cols ], type="lasso", trace=F ), silent=quiet )
        tries <- tries + 1
      } }
    if ( class( lars.obj ) == "try-error" ) return( numeric() )

### Do lars CV and ...
    cv.lars.obj <- try( cv.lars( df.tmp[ cols, ], output[ cols ], K=cv.k, type="lasso", plot.it=F, trace=F ),
                       silent=quiet )
    if ( class( cv.lars.obj ) == "try-error" ) {
      tries <- 1; while( tries <= 20 && class( cv.lars.obj ) == "try-error" ) {
        cv.lars.obj <- try( cv.lars( df.tmp[ cols, ], output[ cols ], K=cv.k, type="lasso", plot.it=F, trace=F ),
                           silent=quiet )
        tries <- tries + 1
      } }
    if ( class( cv.lars.obj ) == "try-error" ) return( numeric() )

### ... choose min CV (either absolute min if cv.choose="min" or min+1se if cv.choose="min+1se")
    min.i <- which.min( cv.lars.obj$cv )
    min.err <- cv.lars.obj$cv.error[ min.i ]

    ##if ( cv.choose[ 1 ] == "min" ) best.s <- which.min( cv.lars.obj$cv )
    ##else
    if ( grepl( "+", cv.choose[ 1 ], fixed=T ) ) {
      se <- as.numeric( gsub( "min+", "", gsub( "se", "", cv.choose[ 1 ] ) ) )
      best.s <- min( which( cv.lars.obj$cv <= min( cv.lars.obj$cv ) + se * min.err ) )
    } else best.s <- which.min( cv.lars.obj$cv )

    orig.coeffs <- coeffs <- coef.lars( lars.obj, s=cv.lars.obj$fraction[ best.s ], mode="fraction" )
    sorted <- names( sort( abs( coeffs ), decreasing=T ) )
    coeffs <- coeffs[ abs( coeffs ) >= coef.cutoff ]
    if ( ! is.na( min.coeffs ) && length( coeffs ) < min.coeffs ) coeffs <- orig.coeffs[ sorted[ 1:min.coeffs ] ]
    if ( ! is.na( max.coeffs ) && length( coeffs ) > max.coeffs ) coeffs <- orig.coeffs[ sorted[ 1:max.coeffs ] ]
    if ( ! quiet ) cat( boot, cv.choose[ 1 ], min.i, min.err, best.s, cv.lars.obj$cv[ best.s ],
                       cv.lars.obj$fraction[ best.s ], length( coeffs ), "\n" )

    if ( rescale.coeffs && length( coeffs ) > 0 ) { ## plug coeffs from LARS back into LM, coeffs will be larger.
      ins <- df.tmp[ ,names( coeffs ), drop=F ]
      coeffs.s <- coef( lm( output[ cols ] ~ ins[ cols, ] - 1 ) )
      names( coeffs.s ) <- names( coeffs )
      coeffs <- coeffs.s[ abs( coeffs.s ) >= coef.cutoff ]
    }
    ##out.coe[[ boot ]] <- coeffs
    if ( boot == 1 ) { out <- list( coeffs=coeffs, lars.obj=lars.obj, cv.lars.obj=cv.lars.obj, best.s=best.s,
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
    tmp2 <- table( names( tmp ) )
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
               se=se, min.err=min.err, all.inputs=rownames( df.tmp ) ) ) ##, params=in.args ) )
  ##return( list( coeffs=coeffs, coeffs.boot=out.coe, coef.quantiles=coef.quantiles ) )
} ## end of inferelator function

get.boot.coef.quantiles <- function( coef.obj, use.grep=F ) {
  coef.quantiles <- NULL
  n.boot <- length( coef.obj$coeffs.boot )
  if ( n.boot > 1 ) {
    tmp <- unlist( coef.obj$coeffs.boot )
    tmp2 <- table( names( tmp ) )
    coef.quantiles <- t( sapply( names( tmp2 ), function( i ) {
      if ( ! use.grep ) tmp3 <- tmp[ names( tmp ) == i ]
      else tmp3 <- tmp[ grepl( i, names( tmp ) ) ]
      tmp3 <- c( tmp3, rep( 0, max( 0, n.boot - length( tmp3 ) ) ) )
      c( n=sum( names( tmp ) == i ) / n.boot, quantile( abs( tmp3 ), prob=c( 0.01, 0.05, 0.1, 0.5, 0.9, 0.95 ) ) ) } ) )
    coef.quantiles <- coef.quantiles[ ! apply( coef.quantiles, 1, function( i ) all( i[ -1 ] == 0 ) ), ]
    print( coef.quantiles, digits=3 )
  }
  invisible( coef.quantiles )
}

###########################################
## mean.variance.normalize
###########################################
mean.variance.normalize <- function( input.matrix, filter=0.04 ) {
  if ( ! is.na( filter ) ) { ## Remove rows that have < 4% nonzero (e.g. env.factors)
    ##which.good <- which( apply( input.matrix, 1, function( i ) mean( i != 0, na.rm=T ) ) >= filter )
    which.good <- which( apply( input.matrix, 1, function( i ) mean( ! is.na( i ), na.rm=T ) ) >= filter )
    input.matrix <- input.matrix[ which.good, ]
  }
  means <- apply( input.matrix, 1, mean, na.rm=T )
  sds <- apply( input.matrix, 1, sd, na.rm=T )
  sds[ sds == 0 | is.na( sds ) ] <- 1
  input.matrix <- apply( input.matrix, 2, "-", means )
  input.matrix <- apply( input.matrix, 2, "/", sds )
  return( input.matrix )
} ## end of mean.variance.normalize function

