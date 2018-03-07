#' @export
postHocTukey=
function(data.in, id.col='entry', rep.col='n_tot', value.col='est_blue', sigma2, sigma2.df){
	t1 = outer(data.in[, value.col], data.in[, value.col], '-'); colnames(t1)=rownames(t1)=data.in[, id.col]
	t2 = sqrt(sigma2) * sqrt(.5 * outer(1/data.in[, rep.col], 1/data.in[, rep.col], '+')); colnames(t2)=rownames(t2)=data.in[, id.col]
	t3 = t1/t2
	pval = ptukey(abs(t3), nrow(t3), sigma2.df, lower.tail=F)
	ordering = as.character(taRifx::sort.data.frame(data.in[, c(id.col, value.col)], f=as.formula(paste('~', '-', value.col)))[, id.col])
	letterReport = multcompLetters(pval[ordering, ordering])
	out = if(length(letterReport) == 3){
				data.frame(id.var=names(letterReport[[2]]), TukeyGrouping=letterReport[[2]])
			}else{
				data.frame(id.var=names(letterReport[[1]]), TukeyGrouping=letterReport[[1]])
			}
	colnames(out)=c(id.col, 'TukeyGrouping')
	return(out)
}

#' @export
postHocHsu=
function(data.in, id.col='entry', rep.col='n_tot', value.col='est_blue', sigma2, sigma2.df){
	tmp.rep = droplevels(data.in[, c(id.col, rep.col)]); colnames(tmp.rep)=c('id', 'rep')
	tmp.ety = droplevels(data.in[, c(id.col, value.col)]); colnames(tmp.ety)=c('id', 'val')
	
	tmp.1 = do.call(rbind,
					lapply(1:nrow(tmp.ety), function(x){
						tmp1 = cbind(tmp.ety[x, ], M=max(tmp.ety[-x, 'val']), m=min(tmp.ety[-x, 'val']))
						tmp2 = cbind(tmp1
									, dM=qt(0.05/(nrow(tmp.ety)), sigma2.df, lower.tail=F)*
											sqrt(sigma2)*sqrt(
														1/subset(tmp.rep, id==tmp1$id, select=rep, drop=T) + 
														1/mean(subset(tmp.rep, id==subset(tmp.ety, val==tmp1$M, select=id, drop=T), select=rep, drop=T))
														)
									, dm=qt(0.05/(nrow(tmp.ety)), sigma2.df, lower.tail=F)*
											sqrt(sigma2)*sqrt(
														1/subset(tmp.rep, id==tmp1$id, select=rep, drop=T) + 
														1/mean(subset(tmp.rep, id==subset(tmp.ety, val==tmp1$m, select=id, drop=T), select=rep, drop=T))
														)
									)
						tmp3 = cbind(tmp2, with(tmp2, cbind(largest=val >= M-dM, smallest=val <= m+dm)))
					}))
	out = tmp.1[, c('id', 'largest', 'smallest')]
	colnames(out) = c(id.col, 'HsuMCB.Largest', 'HsuMCB.Smallest')
	return(out)
}