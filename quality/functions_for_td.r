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
##########################################################################################################################

################# FUNZIONE 17  #############################
##############################################################################################################################

mpd_from_sfs<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 1:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

################# FUNZIONE 17a  #############################
##############################################################################################################################

mpd_from_sfs_singleton<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 2:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

################# FUNZIONE 17b  #############################
##############################################################################################################################

mpd_from_sfs_doubleton<-function(sfs, folded=TRUE){
### check if it is folded or not
if (isTRUE(folded)){n_ind<-2*length(sfs)}
else {n_ind<-length(sfs)}
###
mpd=0
for (i in 3:length(sfs)){
mpd=mpd+(((n_ind-i)*i)*sfs[i])
}
mpd=mpd/((n_ind*(n_ind-1))/2)

return(mpd)
}
#########################################################################################################################

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
#########################################################################################################################




################# FUNZIONE 21  ############################# calcola la Fst di Reynolds et al. 1983
##############################################################################################################################
fst_reynolds_sfs2D<-function(sfs_2D,maff){ ### input uno spettro 2D e un valore di maf. Restituisce un vettore con Fst di Reynolds et al. 1983(unweigheted, weighted)
dati<-sfs_2D
maf<-maff ### definisco maf
n_pop1<-(dim(dati)[2]-1)/2 ## diploidi
n_pop2<-(dim(dati)[1]-1)/2 ## diploidi
somma_as<-0
somma_bs<-0
fst_unw<-c()
for (i in 1:nrow(dati)){
for (j in 1:ncol(dati)){

if (((i+j-2)/(2*(n_pop1+n_pop2)))>=maf){
freq_min_pop1<-(j-1)/(2*n_pop1)
freq_min_pop2<-(i-1)/(2*n_pop2)
freq_min_all<-(i+j-2)/(2*(n_pop1+n_pop2))

bs<-((n_pop1*2*freq_min_pop1*(1-freq_min_pop1))+(n_pop2*2*freq_min_pop2*(1-freq_min_pop2)))/(n_pop1+n_pop2-1)
as<-((4*n_pop1*((freq_min_pop1-freq_min_all)^2))+(4*n_pop2*((freq_min_pop2-freq_min_all)^2))-bs)/(2*((2*n_pop1*n_pop2)/(n_pop1+n_pop2)))

somma_as<-somma_as+(as*dati[i,j])
somma_bs<-somma_bs+(bs*dati[i,j])
if (as+bs==0) {
temp_fst=NA
} else {
temp_fst<-as/(as+bs)
}
fst_unw<-c(fst_unw, rep(temp_fst,dati[i,j]))
}

}
}
fst<-somma_as/(somma_as+somma_bs) # weighted
fst_unw_finale<-mean(fst_unw, na.rm=T) # unweighted
risultato<-c(fst_unw_finale,fst)
return(risultato)
}
####################################################################

################# FUNZIONE 22  ############################# scrivi file GENO
############################################################################################################################## 
fun_geno<-function(data){
for (i in 1:length(data)) {
if (data[i]=="./.") {data[i]<-9}
if (data[i]=="0/0") {data[i]<-0}
if (data[i]=="0/1") {data[i]<-1}
if (data[i]=="1/1") {data[i]<-2}
}
return(data)
}
##############################################################################################################################

################# FUNZIONE 22a ############################# scrivi file GENO mod per NaN
############################################################################################################################## 
fun_geno_mod<-function(data){
for (i in 1:length(data)) {
if (data[i]=="./.") {data[i]<-NaN}
if (data[i]=="0/0") {data[i]<-0}
if (data[i]=="0/1") {data[i]<-1}
if (data[i]=="1/1") {data[i]<-2}
}
return(data)
}
##############################################################################################################################

################# FUNZIONE 23  ############################# conta allele ALT in ogni popolazione
############################################################################################################################## 
################## conta allele ALT in ogni pop. Argomento: lista pop. Oggetto: matrice geno (0,1,2)
fun_conta_ALT_2<-function(data, lista_pop){
##################  prendo gli estremi di ogni popolazione per poter poi calcolare i dati mancanti per pop per filtrare dopo
born_inf_fin<-c()
born_sup_fin<-c()
n_pop<-length(lista_pop)
for (z in 1:n_pop){
####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'é in almeno un individuo
if (z==1) {
born_inf=1
born_sup=length(lista_pop[[z]])
}
else {
vett_lungh_sup<-c()
for (s in 1:z) {
vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
vett_lungh_inf<-c()
for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}

born_inf=1+sum(vett_lungh_inf)
born_sup=sum(vett_lungh_sup)
}
born_inf_fin<-c(born_inf_fin,born_inf) ## i due vettori con gli indici degli estremi delle popolazioni
born_sup_fin<-c(born_sup_fin,born_sup)
}
####################################
risultato<-c()
#for (i in 1:length(data)) {
for (j in 1:length(born_inf_fin)) {
risultato<-c(risultato, sum(data[born_inf_fin[j]:born_sup_fin[j]], na.rm=T))
}
risultato<-c(risultato,sum(data, na.rm=T))
#}
return(risultato)
}
#################################### FINE FUNZ CONTA ALT ####
#########################################################################################################################################

################# FUNZIONE 24  ############################# calcola 2S-SFS pairwise per fastsimcoal. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_sfs2D_pairwise<-function(dati,lista_pop){
n_ind<-2*length(unlist(lista_pop)) ### n° of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF

########## creo matrici di output
new_list<-list()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
righe=2*length(lista_pop[[i]])+1
colonne=2*length(lista_pop[[j]])+1
matrice<-matrix(rep(0,righe*colonne),nrow=righe, ncol=colonne)
print(i)
print(j)
new_list<-list.append(new_list,matrice)
}
}
}
########

for (z in 1:nrow(mat_con_alt)){
if (mat_con_alt[z,ncol(mat_con_alt)]<(n_ind/2)) { ### il MAF é l'allele ALT
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else if (mat_con_alt[z,ncol(mat_con_alt)]>(n_ind/2)) { ### il MAF é l"allele REF
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else {  ### i due alleli hanno freq 0.5 nella pop totale, quindi seguo Laurent che aumenta di 0.5 la entries per entrambi
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+0.5
new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+0.5
indice_list<-indice_list+1
}
}
}
}
}


return(new_list)
}
####################################################################################################################################################################################


################# FUNZIONE 25  ############################# calcola fst di reynolds pairwise. Ha busogno di matrice prodotta da FUNZIONE 23 come input.
##############################################################################################
calcola_fst_da_geno<-function(ogg_geno, lista_pop, maf){
########## creo matrici di output
new_list<-list()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
righe=2*length(lista_pop[[i]])+1
colonne=2*length(lista_pop[[j]])+1
matrice<-matrix(rep(0,righe*colonne),nrow=righe, ncol=colonne)
print(i)
print(j)
new_list<-list.append(new_list,matrice)
}
}
}
########
n_ind<-ncol(ogg_geno)-1
for (z in 1:nrow(ogg_geno)){
if (ogg_geno[z,ncol(ogg_geno)]<(n_ind/2)) { ### il MAF é l'allele ALT
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]<-new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else if (ogg_geno[z,ncol(ogg_geno)]>(n_ind/2)) { ### il MAF é l"allele REF
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]+1
indice_list<-indice_list+1
}
}
}
} else {  ### i due alleli hanno freq 0.5 nella pop totale, quindi seguo Laurent che aumenta di 0.5 la entries per entrambi
indice_list<-1
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
n_ind_pop_row<-2*length(lista_pop[[i]])
n_ind_pop_col<-2*length(lista_pop[[j]])
new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-ogg_geno[z,i]+1,n_ind_pop_col-ogg_geno[z,j]+1]+0.5
new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]<-new_list[[indice_list]][ogg_geno[z,i]+1,ogg_geno[z,j]+1]+0.5
indice_list<-indice_list+1
}
}
}
}
}
fst_parwise_obs<-c()
for (i in 1: length(new_list)){
temp=fst_reynolds_sfs2D(new_list[[i]],maf)[2]
fst_parwise_obs<-c(fst_parwise_obs,temp)
}
return(fst_parwise_obs)
}
##############################################################################################""


################# FUNZIONE 26  ############################# calcola fst pairwise e bootstrap. Ha bisogno di funzione 23 e 24
############################################################################################################################## 
calcola_fst_pairwise_bootstrap<-function(dati,lista_pop, maf, boot){
n_ind<-2*length(unlist(lista_pop)) ### n° of chromosomes, to know who is the MAF
dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
qqq=apply(genotype, 1, fun_geno_mod)
genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
colnames(genotype_num)<-unlist(lista_pop)
mat_con_alt<-t(apply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
pairwise_fst_obs<-calcola_fst_da_geno(mat_con_alt, lista_pop, maf)

fst_bootstrapped<-c()
for (z in 1: boot) {
pairwise_fst_temp<-c()
for (i in 1:(length(lista_pop)-1)){
for (j in 2:length(lista_pop)){
if (i!=j && i<j) {
new_lista_pop_temp<-sample(unlist(list(lista_pop[[i]],lista_pop[[j]])))
new_lista_pop<-list(new_lista_pop_temp[1:length(lista_pop[[i]])],new_lista_pop_temp[(length(lista_pop[[i]])+1):length(new_lista_pop_temp)])
print(new_lista_pop)
mat_con_alt_temp<-t(apply(genotype_num[,new_lista_pop_temp], 1, fun_conta_ALT_2, new_lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
temp<-calcola_fst_da_geno(mat_con_alt_temp, new_lista_pop, maf)###metti maf
pairwise_fst_temp<-c(pairwise_fst_temp, temp)
}
}
}
fst_bootstrapped<-rbind(fst_bootstrapped,pairwise_fst_temp)
}

lista_finale<-list(pairwise_fst_obs, fst_bootstrapped)
return(lista_finale)
}
####################################################################################################################################################################################

###################### CREATE SQUARE DISTANCE MATRIX FROM VECTOR OF DISTANCE COMPUTED FROM MY REYNOLDS FUNCTION
#### for (i in 1:4) {
#### a=1
#### b=N-1
#### mat[,1] <- 1:N-1
####           (N-1)+1: (N-1)+(N-2)
####           (N-1)+(N-2)+1: N-1+(N-2)+(N-3)
####}
create_matrix_Fst_from_vector <- function(fst_pairwise,lista_pop_all){
#### Extract distances from the vector and put in a list
N=length(lista_pop_all)
risss <- list()
a=1
b=N-1
for (i in 1:(N-1)) {
risss<-list.append(risss,fst_pairwise[a:b])
a=b+1
b=b+(N-i-1)
}
#### add 0 to each element of the list to fill the column of the future matrix
for (i in 1:length(risss)){
risss[[i]] <- c(rep(0,i),risss[[i]])
}
#### bind the colomn + the last one with only zeros
almost_final <- c()
for (i in 1:length(risss)){
almost_final <- cbind(almost_final,risss[[i]])
}
almost_final <- cbind(almost_final,rep(0,N))
########### translate col
final <- almost_final+t(almost_final)
return(final)
}
####################################################################################################################################


######## FOLD DAF SPECTRUM from COALA
fold_dafSFS_coala<-function(dati){
dati2<- c(rep(0, (length(dati)+1)/2))
for (i in 1:(length(dati2)-1)) {
dati2[i] <- dati[i] + dati[length(dati)+1-i]
}
dati2[length(dati2)]<-dati[(length(dati)+1)/2]
return(dati2)
}
##########################