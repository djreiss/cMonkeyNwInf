###########################################
## Inf network visual output
## Phu T Van and David J Reiss, ISB
###########################################

## coeffs <- runnit.egrin.data( 66, plot=T, n.boot=1 ); plot.coef.obj( coeffs[[ 66 ]] )
## Plots the lars/cv.lars output plus the network, and pred. vs observed profiles
plot.coeff.obj <- function( coeffs, do.scattersmooth=T, ... ) {
  layout( matrix( c( 1,1,1,2,2,2,4,4,
                     1,1,1,2,2,2,4,4,
                     3,3,3,3,5,5,5,5,
                     3,3,3,3,5,5,5,5,
                     3,3,3,3,5,5,5,5 ),nrow=5,ncol=8,byrow=T ) )
  my.plotCVLars <- function (cv.lars.object, se = TRUE, ...) {
    attach(cv.lars.object)
    plot(fraction, cv, type = "b", ylim = range(cv, cv + cv.error, cv - cv.error), ...)
    if (se) error.bars(fraction, cv + cv.error, cv - cv.error, width = 1/length(fraction))
    detach(cv.lars.object)
    invisible()
  }  
  if ( ! is.null( coeffs$n.boot ) ) n.boot <- 1 ##coeffs$n.boot
  else if ( ! is.null( coeffs[[ 1 ]]$n.boot ) ) {
    n.boot <- coeffs[[ 1 ]]$n.boot  ## bootstrapped at the runnit level
    coeb <- coeffs; coeffs <- coeffs[[ 1 ]]
  }
  pi <- coeffs$plot.info
  require( lars ); require( glmnet )
  if ( "glmnet" %in% class( pi$lars.obj ) ) plot( pi$lars.obj, "lambda" )
  else plot( pi$lars.obj ) ##, main=paste( predicted, collapse=" ", sep=" " ) )
  lines( rep( pi$cv.lars.obj$fraction[ pi$best.s ], 2 ), c( -999, 999 ), col=2, lty=2, lwd=3 )
  my.plotCVLars( pi$cv.lars.obj, se=TRUE, main=class( pi$lars.obj )[ 1 ] )
  if ( "glmnet" %in% class( pi$lars.obj ) ) legend( "bottomleft", pi$cv.choose )
  else legend( "topright", pi$cv.choose )
  lines( rep( pi$cv.lars.obj$fraction[ pi$best.s ], 2 ), c( -999, 999 ), col=2, lty=2, lwd=3 )
  if ( grepl( "+", pi$cv.choose, fixed=T ) ) {
    lines( rep( pi$cv.lars.obj$fraction[ which.min( pi$cv.lars.obj$cv ) ], 2 ),
          rep( min( pi$cv.lars.obj$cv ), 2 ) + c( 0, pi$se * pi$min.err ), col=2, lty=2, lwd=1 )
    lines( c( pi$cv.lars.obj$fraction[ pi$best.s ],
             pi$cv.lars.obj$fraction[ which.min( pi$cv.lars.obj$cv ) ] ),
          rep( min( pi$cv.lars.obj$cv ), 2 ) + pi$se * pi$min.err, col=2, lty=2, lwd=1 )
  }

  if ( length( coeffs$coeffs ) > 0 ) {
    matplot( t( rbind( coeffs$observed[ pi$clust.conds.plot ], pi$predictor.mat[ ,pi$clust.conds.plot ] ) ), 
            col=pi$colors, ylab="Normalized expression", xlab="Conditions", type="l", main=pi$main )
    legend( "bottomright", c( "biclust", names( coeffs$coeffs ) ), col=pi$colors, lty=1, cex=0.5 ) 
    lines( pi$cluster.profile, col="red" )
  } else {
    plot( coeffs$observed[ pi$clust.conds.plot ], col="red", ylab="Normalized expression", xlab="Conditions",
         type="l", main=pi$main ) 
    legend( "bottomright", "biclust", col="red", lty=1, cex=0.5 )
  }

  if ( ! is.null( coeffs$pred.ts ) && nrow( coeffs$pred.ts ) > 1 ) { ## bootstrapped at the inferelator level
    matlines( t( apply( coeffs$pred.ss[ ,pi$clust.conds.plot ], 2, quantile, prob=c( 0.1, 0.9 ) ) ),
             col=rep( "lightblue", 2 ), lty=1, lwd=3 )
    matlines( t( apply( coeffs$pred.ts[ ,pi$clust.conds.plot ], 2, quantile, prob=c( 0.1, 0.9 ) ) ),
             col=rep( "gray", 2 ), lty=1, lwd=3 )
  }

  if ( n.boot > 1 ) { ##exists( 'coeb' ) && length( coeb ) > 1 && coeb[[ 1 ]]$n.boot > 1 ) {
    coeb <- coeb[ sapply( coeb, length ) > 0 ]
    pred.ts <- t( sapply( coeb, "[[", "pred.ts" ) )
    tmp <- t( apply( pred.ts, 2, quantile, prob=c( 0.05, 0.5, 0.95 ) ) )
    rownames( tmp ) <- colnames( coeb[[ 1 ]]$pred.ts )
    matlines( tmp[ pi$clust.conds.plot, ], typ="l", lty=1, col=c( "gray", "red", "gray" ), lwd=3 )
  }

  lines( coeffs$pred.ss[ 1, pi$clust.conds.plot ], col="blue" )
  lines( coeffs$pred.ts[ 1, pi$clust.conds.plot ], col="black" )
  lines( rep( pi$n.conds, 2 ), c( -999, 999 ), col="gray", lty=2, lwd=3 )
  legend( "bottomleft", c( "pred.ss", "pred.ts" ), col=c("blue", "black"), lty=1, cex=0.5 )
  lines( coeffs$observed[ pi$clust.conds.plot ], col="red" ) ## overplot it to be on top
  
  out.net <- plot.cluster.coeffs( list( coeffs ) )

  if ( do.scattersmooth && ! is.null( coeffs$pred.ts ) ) {
    if ( ! exists( "scattersmooth" ) ) source( "~/scratch/halo/generic_scripts/scattersmooth.R" )
    scattersmooth( coeffs$observed[ coeffs$cluster.conds ][ ! is.na( coeffs$pred.ts[ 1, coeffs$cluster.conds ] ) ],
                  coeffs$pred.ts[ 1, coeffs$cluster.conds ][ ! is.na( coeffs$pred.ts[ 1, coeffs$cluster.conds ] ) ] )
  }
  invisible( out.net )
}

## coeffs <- runnit.egrin.data( 66 )
## plot.cluster.coeffs( list( coeffs ) )
## Allows multiple coef sets for multiple biclusters to be included in same network plot

plot.cluster.coeffs <- function( coefs, scale=1, cex=0.5, ... ) {
  require( igraph )
  network <- data.frame()
  comb.cnt <- 1
  node.types <- character()
  for ( coe in coefs ) {
    if ( length( coe$coeffs ) <= 0 ) {
      network <- rbind( network, data.frame( n1=sprintf( "bic%s", coe$k ), n2=sprintf( "bic%s", coe$k ),
                                            weight=NA, mode="-" ) )
    } else {
      for ( i in 1:length( coe$coeffs ) ) {
        n <- strsplit( names( coe$coeffs )[ i ], combine.symbol, fixed=T )[[ 1 ]]
        if ( length( n ) == 1 ) {
          network <- rbind( network, data.frame( n1=n, n2=sprintf( "bic%s", coe$k ),
                                                weight=coe$coeffs[ i ], mode=">" ) )
        } else {
          n2 <- paste( "AND", comb.cnt, sep="" )
          network <- rbind( network, data.frame( n1=n2, n2=sprintf( "bic%s", coe$k ),
                                                weight=coe$coeffs[ i ], mode=">" ) )
          network <- rbind( network, data.frame( n1=n[ 1 ], n2=n2, weight=0, mode="-" ) )
          network <- rbind( network, data.frame( n1=n[ 2 ], n2=n2, weight=0, mode="-" ) )
          comb.cnt <- comb.cnt + 1
        }
      }
    }
    if ( ! is.null( coe$possibly.regulates ) && length( coe$possibly.regulates ) > 0 ) {
      for ( i in 1:length( coe$possibly.regulates ) ) {
        network <- rbind( network, data.frame( n1=names( coe$possibly.regulates )[ i ], n2=sprintf( "bic%s", coe$k ),
                                              weight=0, mode="*" ) )
      }
    }
  }
  gr <- graph.edgelist( as.matrix( network[ ,1:2 ] ), directed=T )
  gr.layout <- layout.fruchterman.reingold.grid( gr, niter=3000 * length( coefs )^2, coolexp=0.5, ... )
  gr.layout <- layout.norm( gr.layout, -1, 1, -1, 1 )
  node.names <- get.vertex.attribute( gr, "name" )
  node.sizes <- rep( 15, length( node.names ) ); names( node.sizes ) <- node.names
  node.sizes[ grepl( "^bic", node.names ) ] <- 25
  node.sizes[ grepl( "^AND", node.names ) ] <- 10
  node.sizes <- node.sizes * scale / length( coefs )
  node.colors <- rep( "lightgreen", length( node.names ) ); names( node.colors ) <- node.names
  node.colors[ grepl( "^bic", node.names ) ] <- "lightblue"
  node.colors[ grepl( "^AND", node.names ) ] <- "gray"
  node.frame.colors <- rep( "black", length( node.names ) ); names( node.frame.colors ) <- node.names
  if ( exists( "predictor.mats" ) ) node.frame.colors[ ! node.names %in% names( predictor.mats$tf.groups ) ] <- "red" ##grepl( "^TFGROUP", node.names ) ] <- "red"
  node.frame.colors[ grepl( "^bic", node.names ) ] <- "blue"
  node.frame.colors[ grepl( "^AND", node.names ) ] <- "gray"
  node.shapes <- rep( "circle", length( node.names ) ); names( node.shapes ) <- node.names
  node.shapes[ grepl( "^bic", node.names ) ] <- "square"
  ##node.shapes[ grepl( "^AND", node.names ) ] <- "triangle"
  node.names[ grepl( "^AND", node.names ) ] <- ""
  node.names <- gsub( "TFGROUP", "tf", node.names )
  edge.colors <- ifelse( is.na( network$weight ), "white", ifelse( network$weight > 0, "red",
                                                                  ifelse( network$weight < 0, "green", "blue" ) ) )
  edge.colors[ as.character( network$mode ) == "*" ] <- "black"
  edge.widths <- abs( network$weight ) * 6 + 0.25; edge.widths[ is.na( edge.widths ) ] <- 0.25
  edge.widths[ as.character( network$mode ) == "*" ] <- 0.25 ##network$n1 %in% names( coe$possibly.regulates ) ] <-
    ##coe$possibly.regulates[ as.character( network$n1[ network$n1 %in% names( coe$possibly.regulates ) ] ) ] / 10
  tmp <- as.character( network$mode ); tmp[ tmp == "*" ] <- "-"; network.mode <- as.factor( tmp )
  plot( gr, layout=gr.layout, axes=F, margin=0, rescale=F, 
       vertex.label=node.names, vertex.size=node.sizes, vertex.color=node.colors,
       vertex.shape=node.shapes, edge.arrow.size=0.5, vertex.frame.color=node.frame.colors,
       vertex.label.cex=cex, ##vertex.label.family="Arial", 
       edge.color=edge.colors, edge.width=edge.widths, edge.arrow.mode=as.character( network$mode ) )
  ##invisible( gr )
  invisible( cbind( network, edge.colors, edge.widths ) )
}

###########################################
## write.cyoscape.files
###########################################
write.cytoscape.files <- function(inf.result, clusterStack, sif.filename){
  ## this function takes an Inferelator result and a cM clusterStack
  ## and outputs a network file and the associated edge and node
  ## attribute files for visualizing the network in Cytoscape

  ## IMPORTANT: this function depends on the output returned by inferelator()
  ## change the below accordingly if the output there changes !

  out<- unlist(inf.result)
  ## write the weights of the TFGROUPs to an edge-attribute file 
  ## at the same time create the actual network .sif
  write("weight (java.lang.Double)", "weights.eda")
  
  gatecount = 1
  for (j in 1:length(out)){
    nodes<- strsplit(names(out)[j], combine.symbol, fixed=T)[[1]]
                                        #print(out[j])
    ## make a .sif file for the actual network
    ## .sif format : [node]<tab>[relationship]<tab>[node]

                                        #print(gatecount)
    ## TODO : load this into a data frame and write.table() or something similar 
    ## if there is only 1 TFGROUP, create two nodes
    if (length(nodes) == 2) {
      write(paste(nodes[2], "activates", nodes[1]), sif.filename, append=T)
      
      write(paste(nodes[2], "(activates)", nodes[1],  "=", out[j]), "weights.eda", append=T)
      
    } else if (length(nodes) == 4){
      ## there are 2 TFGROUP's, create an AND gate (a Y-shaped segment with 4 nodes and 3 edges)
      
      write(paste(nodes[2], "combines", paste("AND-", gatecount, sep="")), sif.filename, append=T)
      write(paste(nodes[3], "combines", paste("AND-", gatecount, sep="")), sif.filename, append=T)
      write(paste(paste("AND-", gatecount, sep=""), "activates", nodes[1]), sif.filename, append=T)
      
      write(paste(paste("AND-", gatecount, sep=""),"(activates)", nodes[1],  "=", out[j]), "weights.eda", append=T)
      
      write(paste(paste("AND-",gatecount,sep=""), "=", "(logicGate)"), "types.noa", append=T)
      gatecount = gatecount + 1
      
      ## there is only a cluster (no significant coeffs from inferelator), just write the cluster as a node to file
    } else if (length(nodes) == 1){
      write(nodes[1], sif.filename, append=T)
    }
    
  }
  ## get attributes of the biclusts from clusterStack (e.g. no. of genes, conds, p-vals)
  ## and write to appropriate node-attribute files
  write("clusterGenes", "clusterGenes.noa")
  write("clusterConditions", "clusterConditions.noa")	
  write("clusterGeneCount", "clusterGeneCount.noa")
  write("clusterConditionCount", "clusterConditionCount.noa")
  write("clusterMotifPValues", "clusterMotifPValues.noa")
  write("clusterMotifs", "clusterMotifs.noa")
  write("clusterResiduals", "clusterResiduals.noa")
  
  for (i in 1:length(clusterStack)){
    write( paste(paste(i, " = ", "(", sep="") ,paste(clusterStack[[i]]$rows, collapse="::"), ")", sep="") , "clusterGenes.noa", append=T)	
    write( paste(paste(i, " = ", "(", sep="") ,paste(clusterStack[[i]]$cols, collapse="::"), ")", sep="") , "clusterConditions.noa", append=T)
    write( paste(i, " = ", clusterStack[[i]]$nrows , sep="") , "clusterGeneCount.noa", append=T)
    write( paste(i, " = ", clusterStack[[i]]$ncols, sep="") , "clusterConditionCount.noa", append=T)
    write( paste(i, " = ", clusterStack[[i]]$e.val, sep="") , "clusterMotifPValues.noa", append=T)
    write( paste(i, " = ", clusterStack[[i]]$resid, sep="") , "clusterResiduals.noa", append=T)
    write( paste(i, " = (cluster)",sep=""), "types.noa", append=T )
    
    ## this bit of code requires pssm.to.string() from cMonkey (in R_scripts/motif_utils.R)
    ## if motif(s) exist, write it in the .noa file
    if(length(clusterStack.redox[[i]]$motif.out$pssms) >0){
      write(paste(i, " = (", gsub(" ", "::", paste(lapply(clusterStack.redox[[i]]$motif.out$pssm, pssm.to.string), collapse=" ")), ")", sep="")
            , "clusterMotifs.noa",append=T)
    } else {
      ## no motifs exist, leave the list blank
      write(paste(i, "=", "()"),	 "clusterMotifs.noa", append=T)
    }	
    
  } ## end of clusterStack traversal loop
  
} ## end of make.network.files function

## count predictors that are in the inferred network
## count.predictors = function(inf_result){
##   out = NULL;
##   for (i in 1:length(inf_result)){
##     tmp = names(inf_result[[i]])
##                                         #strip bcnumber
##     tmp = gsub("^[0-9]+\\.", "", tmp)
##                                         #strip "min" from combo predictors
##     tmp = gsub(".min", "",tmp)
##     out = c(out,tmp)
##   }
##   return(unique(out))
## } ## end of count.predictors function

## plot mean profile of bicluster K (black), predicted steady-state (red), predicted time-series (blue)
## !!!! WARNING : requires mean.old.profiles, tmp3.steady and tmp3.timeseries structures from Inf
## plot.prediction = function(k){
##   plot(mean.old.profiles[k,],type="l",main=paste("bicluster", k, "obs: black, pred.std: red, pred.ts: blue"), xlab="conditions", ylab="normalized expression")
##   lines(tmp3.steady[k,], col="red")
##   lines(tmp3.timeseries[k,], col="blue")
## } ## end of plot.prediction function
