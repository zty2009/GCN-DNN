library(stringr)
setwd('C:/Users/zty20/Documents/R/drug target/demo')
t=read.table('demodata.txt',stringsAsFactor=F)
target=data.frame(unique(t[,2]))

feature=read.table('protein_feature_26.txt',stringsAsFactor=F)
inter=read.table('protein_interact.txt',stringsAsFactor=F)
################sd_proteinnet#########################
index1=apply(target,1,function(x){which(inter[,1]==x)})
index2=apply(target,1,function(x){which(inter[,2]==x)})
sd_proteinnet=matrix(0,nrow(target),nrow(target))
for (i in 1:length(index1)){
	for (j in 1:length(index2)){
		index=intersect(index1[[i]],index2[[j]])
		if (length(index)!=0){
			sd_proteinnet[i,j]=inter[index[1],3]
			sd_proteinnet[j,i]=inter[index[1],3]
			}
		}
}
for (i in 1:nrow(sd_proteinnet)){
	sd_proteinnet[i,i]=1
	}
write.table(file='sd_proteinnet.txt',sd_proteinnet,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(file='sd_proteinname.txt',target,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
################st_proteinnet#########################
index1=apply(target,1,function(x){which(inter[,1]==x)})
index2=apply(target,1,function(x){which(inter[,2]==x)})
index=c(inter[unlist(index1),2],inter[unlist(index2),1])
a=as.data.frame(table(index))
protein=intersect(as.character(a[which(a[,2]>2),1]),feature[,1])
protein=as.matrix(unique(c(protein,as.matrix(target))))

index1=apply(protein,1,function(x){which(inter[,1]==x)})
index2=apply(protein,1,function(x){which(inter[,2]==x)})
st_proteinnet=matrix(0,nrow(protein),nrow(protein))
for (i in 1:length(index1)){
	for (j in 1:length(index2)){
		index=intersect(index1[[i]],index2[[j]])
		if (length(index)!=0){
			st_proteinnet[i,j]=inter[index[1],3]
			st_proteinnet[j,i]=inter[index[1],3]
			}
		}
}
for (i in 1:nrow(st_proteinnet)){
	st_proteinnet[i,i]=1
	}
write.table(file='st_proteinnet.txt',st_proteinnet,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(file='st_proteinname.txt',protein,row.names = FALSE,col.names = FALSE,quote = FALSE)
############sd_proteinfeature###############
index=apply(target,1,function(x){which(feature[,1]==x)[1]})
sd_proteinfeature=feature[index,2:27]
write.table(file='sd_proteinfeature.txt',sd_proteinfeature,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
#############st_proteinfeature###############

index=apply(protein,1,function(x){which(feature[,1]==x)[1]})
st_proteinfeature=feature[index,2:27]
write.table(file='st_proteinfeature.txt',st_proteinfeature,row.names = FALSE,col.names = FALSE,quote = FALSE)

###############st_label####################
drug=data.frame(unique(t[,1]))
st_label=matrix(0,length(protein),nrow(drug))

for (i in 1:length(protein)){
		p=which(t[,2]==protein[i])
		if (length(p)!=0){
			tt=t[p,1]
			position=apply(as.matrix(tt),1,function(x){which(drug==x)})
			st_label[i,position]=1
			}
}
write.table(file='st_label.txt',st_label,row.names = FALSE,col.names = FALSE,quote = FALSE)