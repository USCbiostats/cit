## FDR helpers. fdr.od, iuq, fdr.q.perm and fdr.q.para are unchanged from cit 2.3.2;
## fdr.cit is generalised below.

# fdr function w/ overdispersion parameter
# In order to make CIs estimable, use the conservative approximation that at least 1 positve test was observed among permuted
fdr.od = function ( obsp, permp, pnm, ntests, thres, cl = 0.95, od = NA ) {
    z_ = qnorm(1 - (1 - cl)/2)
    pcount = rep(NA, length(permp))
    for (p_ in 1:length(permp)) {
        permp[[p_]][, pnm] = ifelse(permp[[p_]][, pnm] <= thres, 
            1, 0)
        pcount[p_] = sum(permp[[p_]][, pnm], na.rm = TRUE)
    }
    p = mean(pcount, na.rm = TRUE)/ntests
    e_vr = ntests * p * (1 - p)
    o_vr = var(pcount, na.rm = TRUE)
    if (is.na(od)) {
        od = o_vr/e_vr
        if (!is.na(od)) 
            if (od < 1) 
                od = 1
    }
    if (is.na(od)) od = 1
    nperm = length(permp)
    mo = ntests
    ro = sum(obsp <= thres)
    vp = sum(pcount)
    vp1 = vp
    rslt = rep(NA, 4)
    if (ro > 0) {
        if (vp == 0) 
            vp = 1
        mean.vp = vp / nperm
        fdr0 = mean.vp / ro
        pi0 = (mo - ro)/(mo - (vp/nperm))
        if( is.na(pi0) ) pi0 = 1
        if( pi0 < 0.5 ) pi0 = 0.5    # updated pi0 to limit its influence
        if( pi0 > 1 ) pi0 = 1
        fdr = fdr0 * pi0    # updated calculation of fdr to be robust to ro = mtests
        
        # variance of FDR
        mp = nperm * mo
        t1 = 1 / vp
        denom = mp - vp
        t2 = 1 / denom
        t3 = 1 / ro
        denom = ntests - ro
        if( denom < 1 ) denom = 1
        t4 = 1 /  denom
        s2fdr = (t1 + t2 + t3 + t4) * od
        ul = exp(log(fdr) + z_ * sqrt(s2fdr))
        ll = exp(log(fdr) - z_ * sqrt(s2fdr))
        
        rslt = c(fdr, ll, ul, pi0)
        rslt = ifelse(rslt > 1, 1, rslt)
        rslt = c(rslt, od, ro, vp1)
        names(rslt) = c( "fdr", "fdr.ll", "fdr.ul", "pi.0", "od", "s.obs", "s.perm" )
    }
    return(rslt)
} # End fdr.od


# function to combine q-values into an omnibus q-value that represents the intersection of alternative hypotheses and the union of null hypotheses
iuq = function( qvec ){
	qvec1 = 1 - qvec
	tmp = 1
	for( i in 1:length(qvec1) ) tmp = tmp * qvec1[ i ]
	qval = 1 - tmp
	return( qval )
} # End iuq

# wrapper function for fdr.od, gives q-values for input observed and permuted data
fdr.q.perm = function(obs.p, perml, pname, ntests, cl=.95, od=NA){
	# set q.value to minimum FDR for that p-value or larger p-values
	m = length(obs.p)
	new.order = order(obs.p)
	po = obs.p[new.order]
	qvals = rep(NA, m)
	for( tst in 1:m ){
		thresh = po[ tst ]
		thresh = ifelse( is.na(thresh), 1, thresh )
		thresh = ifelse( is.null(thresh), 1, thresh )
		if( thresh < 1){
			qvals[ tst ] = fdr.od(obs.p, perml, pname, ntests, thresh, cl=cl, od=od)[1]
			qvals[ 1:tst ] = ifelse( qvals[ 1:tst ] > qvals[ tst ], qvals[ tst ], qvals[ 1:tst ] )
		} else qvals[ tst ] = 1
	} # End tst loop
	qvals1 = qvals[order(new.order)]
	return( qvals1 )
} # End fdr.q.perm


# Millstein FDR (2013) parametric estimator, gives q-values
fdr.q.para = function( pvals ){
	# set q.value to minimum FDR for that p-value or larger p-values
	m = length( pvals )
	new.order = order(pvals)
	po = pvals[new.order]
	qvals = rep(NA, m)
	for( tst in 1:m ){
		thresh = po[ tst ]
		if( thresh > .99 ) qvals[ tst ] = 1
		if( thresh < 1 ){
			S = sum( pvals <= thresh )
			Sp = m * thresh
			prod1 = Sp / S
			prod2 = (1 - S/m) / (1 - Sp/m) 
			prod2 = ifelse(is.na(prod2), .5, prod2)
			prod2 = ifelse(prod2 < .5, .5, prod2)
			qvals[ tst ] = prod1 * prod2
		} # End if thresh
		qvals[ 1:tst ] = ifelse( qvals[ 1:tst ] > qvals[ tst ], qvals[ tst ], qvals[ 1:tst ] )
	} # End for tst
	qvals1 = qvals[order(new.order)]
	qvals1 = ifelse(qvals1 > 1, 1, qvals1)
	return( qvals1 )
} # End fdrpara



# compute FDR qvalues from output of cit.bp, cit.cp or cit.cox, organized in a list with each element the output for a specific test.
# Generalised from 2.3.2: uses whichever component p-values are present (cit.cox has no test 3)
# and accepts the name p_TindLgvnG used by development versions for test 4.
fdr.cit = function( cit.perm.list, cl=.95, c1=NA ){

	all.pnms = c( "p_TassocL","p_TassocGgvnL","p_GassocLgvnT","p_LindTgvnG" )
	all.qnms = list( c("q.TaL", "q.ll.TaL", "q.ul.TaL"),
		c("q.TaGgvL", "q.ll.TaGgvL", "q.ul.TaGgvL"),
		c("q.GaLgvT", "q.ll.GaLgvT", "q.ul.GaLgvT"),
		c("q.LiTgvG", "q.ll.LiTgvG", "q.ul.LiTgvG") )
	cit.perm.list = lapply( cit.perm.list, function(x){
		x = as.data.frame( x )
		names(x)[ names(x) == "p_TindLgvnG" ] = "p_LindTgvnG"
		x
	} )
	present = all.pnms %in% names( cit.perm.list[[1]] )
	pnms = all.pnms[ present ]
	pnm.lst = all.qnms[ present ]

	nperm = nrow(cit.perm.list[[1]]) - 1
	if( nperm < 1 ) stop( "fdr.cit needs results from cit.cp, cit.bp or cit.cox run with n.perm > 0." )
	ntest = length( cit.perm.list )
	perml = vector('list', nperm )
	obs = as.data.frame(matrix( NA, nrow=0, ncol=ncol(cit.perm.list[[1]] ) ) )
	names(obs) = names( cit.perm.list[[1]] )
	for( i in 1:ntest ){
		obs[ i, ] = cit.perm.list[[ i ]][ 1, ]
		for( j in 1:nperm ){
			if( i == 1 ) perml[[ j ]] = obs[ 0, ]
			perml[[ j ]][ i, ] = cit.perm.list[[ i ]][ j+1, ]
		}
	}
	## set 0 p-values to 1e-16
	for( pnm in pnms ) obs[, pnm ] = ifelse( obs[, pnm ] < 1e-16, 1e-16, obs[, pnm ] )
	for( perm in 1:nperm ){
		for( pnm in pnms ) perml[[ perm ]][, pnm ] = ifelse( perml[[ perm ]][, pnm ] < 1e-16, 1e-16, perml[[ perm ]][, pnm ] )
	}

	fdrmat = as.data.frame( matrix( NA, nrow=0, ncol=4 + 3 * length(pnms) ) )
	names( fdrmat ) = c( "p.cit", "q.cit", "q.cit.ll", "q.cit.ul", unlist( pnm.lst ) )

	for( tst in 1:nrow(obs) ){
		for( pind in 1:length(pnms) ) {
			pname = pnms[ pind ]
			cutoff = obs[ tst, pname ]
			cutoff = ifelse( is.na(cutoff), 1, cutoff )
			cutoff = ifelse( is.null(cutoff), 1, cutoff )
			if( cutoff < 1){
				fdrmat[ tst, pnm.lst[[ pind ]]  ] = fdr.od(obs[, pname], perml, pname, nrow(obs), cutoff, cl=cl, od=c1)[ 1:3 ]
			} else fdrmat[ tst, pnm.lst[[ pind ]]  ] = c(1,1,1)
		}
	}

	fdrmat[ , pnms  ] = obs[, pnms ]

	# make each component q-value monotone in its p-value
	for( pind in 1:length(pnms) ){
		qn = pnm.lst[[ pind ]]
		op = order(fdrmat[ , pnms[ pind ] ])
		for(tst in 1:nrow(fdrmat)){
			aa = fdrmat[ op[1:tst], qn[1] ] > fdrmat[ op[tst], qn[1] ]
			for( q_ in qn ) fdrmat[ op[1:tst], q_ ] = ifelse( aa, fdrmat[ op[tst], q_ ], fdrmat[ op[1:tst], q_ ] )
		}
	}

	# p.cit
	qcols = sapply( pnm.lst, `[`, 1 ); llcols = sapply( pnm.lst, `[`, 2 ); ulcols = sapply( pnm.lst, `[`, 3 )
	for( tst in 1:nrow(obs) ){
		fdrmat[ tst, "p.cit"  ] = obs[ tst, "p_cit" ]
		fdrmat[ tst, "q.cit"  ] = iuq( fdrmat[ tst, qcols ] )
		fdrmat[ tst, "q.cit.ll"  ] = iuq( fdrmat[ tst, llcols ] )
		fdrmat[ tst, "q.cit.ul"  ] = iuq( fdrmat[ tst, ulcols ] )
	}

	op = order(fdrmat[ , "p.cit" ])
	for(tst in 1:nrow(fdrmat)){
		aa = fdrmat[ op[1:tst], "q.cit" ] > fdrmat[ op[tst], "q.cit" ]
		fdrmat[ op[1:tst], "q.cit" ] = ifelse( aa, fdrmat[ op[tst], "q.cit" ], fdrmat[ op[1:tst], "q.cit" ] )
		fdrmat[ op[1:tst], "q.cit.ll" ] = ifelse( aa, fdrmat[ op[tst], "q.cit.ll" ], fdrmat[ op[1:tst], "q.cit.ll" ] )
		fdrmat[ op[1:tst], "q.cit.ul" ] = ifelse( aa, fdrmat[ op[tst], "q.cit.ul" ], fdrmat[ op[1:tst], "q.cit.ul" ] )
	}

	return( fdrmat )
} # End fdr.cit
