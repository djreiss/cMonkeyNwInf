###########################################
## predictelator
## Phu T Van and David J Reiss, ISB
###########################################

###########################################
## build.profiles.from.weights
###########################################
## !!!! VERY IMPORTANT: IF FORMAT OF COEFFS OUTPUT CHANGES, MUST CHANGE THIS FUNCTION ACCORDINGLY !!!!
## TODO: make a predictelate.coeff.obj function that can do predictelating on a coeff
##    obj... allowing for knockout or knock-in of any predictor(s).
predictelate <- function( cluster.rows, coeffs, ratios, predictor.mats=NULL, tf.groups=NULL, col.map=NULL,
                         tau=10, max.coeffs=length( coeffs ), ... ) {
  ## this function takes a nested list of coeffs from LARS and construct new gene expression profiles
  ## if a list of tf.groups are provided, builds profiles for them too
  ## Note can use predictor.mats as pre-calculated by inferelate() but if that doesn't exist, it will
  ##   generage the predictor profiles (including combo's) on the fly.
  if ( length( coeffs ) <= 0 ) { out <- ratios[ 1, ] * 0; return( out ) }
  if ( max.coeffs < length( coeffs ) ) coeffs <- sort( coeffs, decreasing=T )[ 1:max.coeffs ]

  coeff.names <- unique( unlist( strsplit( names( coeffs ), combine.symbol, fixed=T ) ) )

  if ( is.null( predictor.mats ) && ! is.null( tf.groups ) && any( coeff.names %in% names( tf.groups ) ) ) {
    tfgroup.ratios <- t( sapply( tf.groups[ which( names( tf.groups ) %in% coeff.names ) ],
                                function( i ) apply( ratios[ i, ,drop=F ], 2, mean ) ) )
    rownames( tfgroup.ratios ) <- names( tf.groups )[ names( tf.groups ) %in% coeff.names ]
    ## tfgroup.ratios is a matrix of the ratios for the Tf.Groups
    ratios <- rbind( ratios, tfgroup.ratios )
  } else {
    ratios <- rbind( ratios, predictor.mats$predictor.mat[ ,colnames( ratios ) ],
                    predictor.mats$predictor.mat.ands[ ,colnames( ratios ) ] )
  }
  out.ss <- 0

  for ( j in 1:length( coeffs ) ) { ##[[i]]$coeffs)){
    if ( coeffs[ j ] == 0 ) next
    nodes <- unlist( strsplit( names( coeffs )[ j ], combine.symbol, fixed=T ) ) ##[[i]]$coeffs[j]), "\\."))
    ## there is a single predictor, multiply the predictor profile by the weight
    if ( length( nodes ) == 1 ) { ##2) {
      ##cat(paste("Single predictor", nodes[2], "regulates", nodes[1], "with a weight of", weight, "\n") )
      if ( ! nodes[ 1 ] %in% rownames( ratios ) ) next
      tmp <- ratios[ nodes, ] ##* weight	   	
    } else if ( length( nodes ) == 3 ) { ##){			    	
      ## there are 2 predictors, get their weighted ratios as well as the combining function (min or max)
      ##cat(paste(nodes[2], nodes[4], "with", nodes[3], "regulates", nodes[1], "with a weight of", weight, "\n"))
      if ( names( coeffs )[ j ] %in% rownames( ratios ) ) {
        tmp <- ratios[ names( coeffs )[ j ], ]
      } else {
        if ( ! all( nodes[ 1:2 ] %in% rownames( ratios ) ) ) next
        tmp <- apply( ratios[ c( nodes[ 1 ], nodes[ 2 ] ), ], 2, FUN=nodes[ 3 ], na.rm=T ) ##func) ##* weight
      }
    }
    tmp[ is.na( tmp ) ] <- 0 ## Ignore NAs in the sum
    out.ss <- out.ss + tmp * coeffs[ j ] ##[[i]]$coeffs[j] ##weight
  }
  out <- out.ss
  
  if ( ! is.null( col.map ) && ! is.na( tau ) && tau > 0 ) {
    ## if have col.map (i.e. we have time-series data),  put in the time component
    conds <- colnames( ratios )
    prevs <- as.character( col.map[ conds, "prevCol" ] )
    
    del.ts <- as.numeric( as.character( col.map[ conds, "delta.t" ] ) )
    del.ts[ del.ts < 1 ] <- 1
    
    ## trim the new profiles to include only the biclusts that have non-zero solutions
    
    ## OLD: out.tmp <- tmp3[,conds ] - tau * ( tmp3[ ,nexts ] - tmp3[ ,conds ]) / del.ts
    ##y_m = ( tau * y_{m-1} + del.t * sum(beta*z) ) / ( tau + del.t )
    ##  where out.ss is the sum(beta*z)
    cluster.prof <- apply( ratios[ cluster.rows, ,drop=F ], 2, mean, na.rm=T )
    tmp1 <- cluster.prof[ prevs ]; tmp1[ is.na( tmp1 ) ] <- 0 ## Ignore NAs in the sum
    tmp2 <- out.ss[ conds ]; tmp2[ is.na( tmp2 ) ] <- 0 ## Ignore NAs in the sum
    out.ts <- ( tau * tmp1 + del.ts * tmp2 ) / ( tau + del.ts )

    names( out.ts ) <- conds
    out <- out.ts
  }
  out[ out > 3.0 ] <- 3.0
  out[ out < -3.0 ] <- -3.0
  out
} ## end of predictelate function
