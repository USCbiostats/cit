## Unchanged from cit 2.3.2.

##### Function to compute F test given continuous outcome and full vs reduced sets of covariates
linreg = function( nms.full, nms.redu=NULL, nm.y, mydat ){
   
   mydat = na.exclude( mydat )
   
   vrs.2 = paste( nms.full, collapse="+" )
   formula2 = paste( nm.y, " ~ ", vrs.2, sep="")
   fit.full = glm( formula2 , data=mydat )
   
   if( is.null( nms.redu ) ){
      formula2 = paste( nm.y, " ~ 1 ", sep="")
      fit.redu  = glm( formula2 , data=mydat )
   } else {
   
      vrs.1 = paste( nms.redu, collapse="+" )
      formula1 = paste( nm.y, " ~ ", vrs.1, sep="")
      fit.redu  = glm( formula1 , data=mydat )
      
   } # End if null redu
   
   tmp = anova( fit.full, fit.redu, test="F" )
   pval.f = tmp$"Pr(>F)"[2]
   return( pval.f )

} # End function linreg
