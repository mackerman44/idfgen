#
# Calculate population genetic F statistics
#
# Intraclass estimators following Pons and Chaouche (1995)
# and ANOVA estimators following Cockerham and Weir (1984)
# as described in Excoffier (2001)
#
# Data is of form: pop loc1 loc2...locN
#
# Genotype at locus is string eg "241243" a la FSTAT or "99/103"
#
# 20081103 Fixed bug in estimation of H_t
#  < ht[n]<- 2*(sum(1-p2[lower.tri(p2)]))/npop/(npop-1)
#  > ht[n]<- 2*(1-sum(p2[lower.tri(p2)]))/npop/(npop-1)
#


fstats <- function (dat) {
  nloc<-ncol(dat)-1
  pops<-table(dat[,1])
  npop<-length(pops)
  Fis<-Fit<-Fst<-ht<-hs<-h0<-double(nloc)
  n<-0
  for (i in 2:ncol(dat)) {
    n<-n+1
    g<-as.character(dat[!is.na(dat[,i]),i])
    p<-dat[!is.na(dat[,i]),1]
    div<-regexpr("/",g)
    wid<-nchar(g)
    div[div<0]<-wid[div<0]/2+0.5
    g1<-substr(g,1,ceiling(div)-1)
    g2<-substr(g,div+1,wid)
    hom<-as.integer(g1==g2)
    p.tab<-table(c(p,p),c(g1,g2))
    p.tab<-p.tab/apply(p.tab,1,sum)
    p2 <- p.tab %*% t(p.tab)
    h.tab<-table(p,hom)
    h.tab<-h.tab[,2]/apply(h.tab,1,sum)
    h0[n]<-1-mean(h.tab)
    hs[n]<- (pops/(pops-1)) %*% (1-apply(p.tab^2,1,sum)-h0[i-1]/(2*pops)) / npop
    ht[n]<- 2*(sum(1-p2[lower.tri(p2)]))/npop/(npop-1)
  }
  Fis<-(hs-h0)/hs
  Fit<-(ht-h0)/ht
  Fst<-(ht-hs)/ht
  allstats<-as.data.frame(round(cbind(h0,hs,ht,Fis,Fit,Fst),3))
  rownames(allstats)<-colnames(dat)[2:ncol(dat)]
  allstats
}
boot.fstats <- function (dat,R=300) {
  require(boot)
  fsi <- function(dat,idx) {
           fstats(dat[idx,])
  }
  boot(dat,statistic=fsi,R=R)
}
boot.fis <- function (dat,R=300) {
  require(boot)
  fis <- function(dat,idx) {
           fstats(dat[idx,])[,4]
  }
  boot(dat,statistic=fis,R=R)
}
boot.fit <- function (dat,R=300) {
  require(boot)
  fit <- function(dat,idx) {
           fstats(dat[idx,])[,5]
  }
  boot(dat,statistic=fit,R=R)
}
boot.fst <- function (dat,R=300) {
  require(boot)
  fst <- function(dat,idx) {
           fstats(dat[idx,])[,6]
  }
  boot(dat,statistic=fst,R=R)
}
cockerham.loc <- function (dat,loc) {
  g<-as.character(dat[!is.na(dat[,loc]),loc])
  deme<-dat[!is.na(dat[,loc]),1]
  pops<-table(deme)
  npop<-length(pops)
  n<-length(deme)
  div<-regexpr("/",g)
  wid<-nchar(g)
  div[div<0]<-wid[div<0]/2+0.5
  g1<-substr(g,1,ceiling(div)-1)
  g2<-substr(g,div+1,wid)
  g<-c(g1,g2)
  deme<-factor(c(deme,deme))
  indiv<-factor(c(1:n,1:n))
  alleles<-table(g,deme)
  allele.names<-as.integer(rownames(alleles))
  n.all<-length(allele.names)
  Capf<-Theta<-Smallf<-va<-vb<-vw<-double(n.all)
  nprime<- (n-sum(pops^2)/n)/max(1,npop-1)
  for(s in 1:n.all) {
    y<-(g==allele.names[s])
    m<-aov(y ~ deme+indiv)
    ms<-anova(m)[,3]
    vw[s]<-ms[3]
    vb[s]<-0.5*(ms[2]-ms[3])
    va[s]<-0.5*(ms[1]-ms[2])/nprime
  }
  m.va<-mean(va)
  m.vb<-mean(vb)
  m.vw<-mean(vw)
  return(alleles=as.data.frame(cbind(alleles, var.stats(va,vb,vw),va,vb,vw)),
         total=as.data.frame(cbind(t(pops),n.all,var.stats(m.va,m.vb,m.vw),
                                   va=m.va,vb=m.vb,vw=m.vw)))
}
var.stats <- function(va,vb,vw) {
  Capf<-(va+vb)/(va+vb+vw)
  Theta<-va/(va+vb+vw)
  Smallf<-vb/(vb+vw)
  cbind(Capf,Theta,Smallf)
}
cockerham.tot <- function (dat)  {
  var.ests<-cockerham.loc(dat,2)$alleles
  for (i in 3:ncol(dat)) {
    var.ests<-rbind(var.ests,cockerham.loc(dat,i)$alleles)
  }
  m.va<-mean(var.ests$va)
  m.vb<-mean(var.ests$vb)
  m.vw<-mean(var.ests$vw)
  return(c(var.stats(m.va,m.vb,m.vw),va=m.va,vb=m.vb,vw=m.vw))
}
cockerham <- function (dat)  {
  require(boot)
  var.ests<-cockerham.loc(dat,2)$total
  for (i in 3:ncol(dat)) {
    var.ests<-rbind(var.ests,cockerham.loc(dat,i)$total)
  }
  boot.c <- boot(var.ests, statistic=global.stats, R=200)
  return(markers=var.ests, total=boot.c)
}
global.stats<-function (ests,idx) {
  data.idx<-ests[idx,]
  m.va<-sum(data.idx$n.all*data.idx$va)/sum(data.idx$n.all)
  m.vb<-sum(data.idx$n.all*data.idx$vb)/sum(data.idx$n.all)
  m.vw<-sum(data.idx$n.all*data.idx$vw)/sum(data.idx$n.all)
  var.stats(m.va,m.vb,m.vw)
}
anova.strp <- function (dat) {
  nloc<-ncol(dat)-1
  pops<-table(dat[,1])
  n<-sum(pops)
  npop<-length(pops)
  Ris<-Rit<-Rst<-va<-vb<-vw<-double(nloc)
  s<-0
  for (i in 2:ncol(dat)) {
    s<-s+1
    g<-as.character(dat[!is.na(dat[,i]),i])
    deme<-dat[!is.na(dat[,i]),1]
    nobs<-length(deme)
    nprime<- (n-sum(pops^2)/n)/max(1,npop-1)
    div<-regexpr("/",g)
    wid<-nchar(g)
    div[div<0]<-wid[div<0]/2+0.5
    g1<-substr(g,1,ceiling(div)-1)
    g2<-substr(g,div+1,wid)
    alength<-as.integer(c(g1,g2))
    deme<-factor(c(deme,deme))
    indiv<-factor(c(1:nobs,1:nobs))
    print(length(alength))
    print(length(deme))
    print(length(indiv))
    m<-aov(alength ~ deme+indiv)
    ms<-anova(m)[,3]
    vw[s]<-ms[3]
    vb[s]<-0.5*(ms[2]-ms[3])
    va[s]<-0.5*(ms[1]-ms[2])/nprime
  }
  cbind(va,vb,vw)
}

