################### 23/03/26 ######################
############### Useful functions ##################
############## Achille RIVOLTELLA #################



# 1) FST HUDSON 
# 2) calcola_normalized_foldedSFS
# 3) calcola_TD_folded
# 4) sliding_fst_spectre


####################################################### FST Hudson from Pierre Lesturgie. Pairwise and per site
##########################################################################################################################
fst.hudson <- function(vcf, pop_list){
  require(rlist)
  require(vcfR)
  
  gt2alleles = function(gt){
    alleles = unlist(strsplit(gt,"[/|]"))
    alleles[which(alleles==".")] = NA
    while(length(alleles) != length(gt)*2){
      alleles <- c(alleles,NA)
    }
    return(as.numeric(alleles))
  }
  
  
  mpd = function(an,ac){
    n_pairs = an * (an - 1) / 2
    n_same = rowSums(ac * (ac - 1) / 2)
    
    n_diff = n_pairs - n_same
    
    return(n_diff / n_pairs)
  }
  
  
  mpd.between = function(an1,an2,ac1,ac2){
    n_pairs = an1 * an2 
    n_same = rowSums(ac1 * ac2)
    
    n_diff = n_pairs - n_same
    
    return(n_diff / n_pairs)
  }
  
  GT = extract.gt(vcf,as.numeric = F,convertNA = F)
  
  mx = max(as.numeric(unlist(strsplit(max(GT),"[/|]"))))
  
  a.counts = function(site){
    c <- c()
    for(i in 1:(mx+1)){
      c <- c(c,length(which(site == (i-1))))
    }
    return(c)
  }
  
  alleles_counts = allele_num = list()
  for(i in 1:length(pop_list)){
    alleles = t(apply(GT[,pop_list[[i]]],1,gt2alleles))
    
    alleles_counts[[i]] = t(apply(alleles,1,a.counts))
    
    allele_num[[i]] = rowSums(alleles_counts[[i]])
  }
  
  
  fst_nuc = list()
  fst = nms = c()
  for(i in 1:(length(pop_list)-1)){
    for(j in (i+1):length(pop_list)){
      within = (mpd(allele_num[[i]],alleles_counts[[i]]) +
                  mpd(allele_num[[j]],alleles_counts[[j]]))/2
      between = mpd.between(allele_num[[i]],allele_num[[j]],alleles_counts[[i]],alleles_counts[[j]])
      fst_nuc = list.append(fst_nuc,(between - within)/between)
      fst = c(fst,(sum(between,na.rm = T) - sum(within,na.rm = T))/sum(between,na.rm = T))
      nms = c(nms,paste0(names(pop_list)[i],"/",names(pop_list)[j]))
    }
  }
  names(fst_nuc) = nms
  fst = data.frame(FST=fst)
  rownames(fst) = nms
  
  result = list(fst,fst_nuc)
  names(result) = c("pairwise_fst","pairwise_fst_nuc")
  return(result)
}
########################################################################################



################# FUNZIONE 15  #############################
##############################################################################################################################

calcola_normalized_foldedSFS<-function(vettore_sfs){
  somma<-sum(vettore_sfs)
  eta_2<-vettore_sfs/somma # normalizzo per numero di SNPs
  ind<-length(eta_2)
  
  ### trasformo per ottenere una curva flat per una pop costante (formula Lapierre et al. 2017, Genetics)
  for (i in 1:(ind-1)){
    eta_2[i]<-eta_2[i]*i*((2*ind)-i)/(2*ind)
  }
  eta_2[ind]<-eta_2[ind]*ind
  eta_2<-as.vector(eta_2)
  
  asse_x<-as.vector(seq(1,ind, by=1)/(2*ind)) #### calcolo gli i/2N per plottare il folded normalizzato
  ssss=as.matrix(t(rbind(asse_x,eta_2)))
  
  return(ssss)
}



################# FUNZIONE 20  ############################# calcola la D di tajima a partire d'un SFS folded
##############################################################################################################################
calcola_TD_folded<-function(folded_sfs){  #### ha bisogno della funzione calcola_MPD
  sample_size<-2*length(folded_sfs)
  theta_P<-mpd_from_sfs(folded_sfs)
  S<-sum(folded_sfs)
  
  ###### 
  calcola_a1<-function(sample_size){
    a<-0
    for (i in 1:(sample_size-1)){
      a<-a+(1/i)
    }
    return(a)
  }
  ##################
  
  ###### 
  calcola_a2<-function(sample_size){
    a<-0
    for (i in 1:(sample_size-1)){
      a<-a+(1/(i*i))
    }
    return(a)
  }
  ##################
  
  
  ###### 
  calcola_b1<-function(sample_size){
    
    a<-(sample_size+1)/(3*(sample_size-1))
    
    return(a)
  }
  ##################
  
  ###### 
  calcola_b2<-function(sample_size){
    
    a<-(2*((sample_size*sample_size)+sample_size+3))/(9*sample_size*(sample_size-1))
    
    return(a)
  }
  ##################
  
  ###### 
  calcola_c1<-function(sample_size){
    a<-calcola_b1(sample_size)-(1/calcola_a1(sample_size))
    return(a)
  }
  ##################
  
  
  ###### 
  calcola_c2<-function(sample_size){
    a<-calcola_b2(sample_size)-((sample_size+2)/(calcola_a1(sample_size)*sample_size))+(calcola_a2(sample_size)/(calcola_a1(sample_size)^2))
    
    return(a)
  }
  ##################
  
  calcola_e1<-function(sample_size){
    a<-calcola_c1(sample_size)/calcola_a1(sample_size)
    return(a)
  }
  
  calcola_e2<-function(sample_size){
    a<-calcola_c2(sample_size)/((calcola_a1(sample_size)^2)+calcola_a2(sample_size))
    return(a)
  }
  
  TD<-(theta_P-(S/calcola_a1(sample_size)))/(((calcola_e1(sample_size)*S)+(calcola_e2(sample_size)*S*(S-1)))^0.5)
  risultati<-list(TD,theta_P,S/calcola_a1(sample_size),S)
  return(risultati)
}




################# FUNZIONE 27  ############################# SLIDING WINDOWS Fst pairwise e theta x tutte pop de lista_pop
############################################################################################################################## 
sliding_fst_spectre<-function(dati, lista_pop, window, slide, maf){
  dim_locus<-window
  interval<-slide
  vettore_snps_pos<-getPOS(dati)
  n_pop<-length(lista_pop)
  low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
  upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
  per_plot<-(upper_bound+low_bound)/2
  
  ris_fst_sliding<-c()
  ris_spectre_sliding<-c()
  for (i in 1:length(low_bound)) {
    sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
    if (length(sss)>0) {
      dati_temp<-dati[sss,]
      fst_temp<-calcola_fst_pairwise_nobootstrap(dati_temp,lista_pop, maf) ### METTERE MAF IN FUNZIONE!
      ris_fst_sliding<-rbind(ris_fst_sliding, fst_temp)
      
      temp_spectre<-c()
      for (z in 1:length(lista_pop)){
        dati_bin<-vcfR2DNAbin(dati_temp[, c("FORMAT",lista_pop[[z]])], extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
        sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
        temp_spectre<-c(temp_spectre,sfs_tot)
      }
      ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
    } else {
      ris_fst_sliding<-rbind(ris_fst_sliding, rep(NaN,n_pop*(n_pop-1)/2))
      ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,length(unlist(lista_pop))))
    }
  }
  
  finale_fst<-cbind(per_plot,ris_fst_sliding)
  finale_fst<-na.omit(finale_fst)
  
  finale_spectre<-cbind(per_plot,ris_spectre_sliding)
  finale_spectre<-na.omit(finale_spectre)
  
  lista_risultati<-list(finale_fst, finale_spectre)
  
  return(lista_risultati)
}


########################## FOLD DAF SPECTRUM from COALA ##################################
fold_dafSFS_coala<-function(dati){
  dati2<- c(rep(0, (length(dati)+1)/2))
  for (i in 1:(length(dati2)-1)) {
    dati2[i] <- dati[i] + dati[length(dati)+1-i]
  }
  dati2[length(dati2)]<-dati[(length(dati)+1)/2]
  return(dati2)
}