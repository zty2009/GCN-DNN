library(stringr)
setwd('C:/Users/zty20/Documents/R/drug target/demo')
t=read.table('demodata.txt',stringsAsFactor=F)
target=unique(t[,2])
drug=unique(t[,1])
con <- file("drug_detail2.txt", "r")
line=readLines(con,n=60555)
close(con)
database=line[seq(1,60555,5)]

index=apply(as.matrix(drug),1,function(x){which(x==database)})

##############sd drug net################
drugnet=line[(index-1)*5+4]
net=str_split(drugnet,' ')
pt=unlist(net)
a=as.data.frame(table(pt))
alldrug=as.character(a[which(a[,2]>8),1])
alldrug=unique(c(alldrug,drug))


database2=line[seq(4,60555,5)]
dindex=apply(as.matrix(alldrug),1,function(x){grep(x,database)})


dindex=unlist(dindex)	
index=apply(as.matrix(alldrug),1,function(x){grep(x,database2)})
sd_drugnet=matrix(0,length(alldrug),length(alldrug))
for (i in 1:length(dindex)){
	a=intersect(index[[i]],dindex)
	position=apply(as.matrix(a),1,function(x){which(x==dindex)})
	sd_drugnet[i,position]=1
	}

write.table(file='sd_drugnet.txt',sd_drugnet,row.names = FALSE,col.names = FALSE,quote = FALSE)
###############st drug net####################
dindex=apply(as.matrix(drug),1,function(x){grep(x,database)})
index=apply(as.matrix(drug),1,function(x){grep(x,database2)})
st_drugnet=matrix(0,length(drug),length(drug))
for (i in 1:length(dindex)){
	a=intersect(index[[i]],dindex)
	position=apply(as.matrix(a),1,function(x){which(x==dindex)})
	st_drugnet[i,position]=1
	}

write.table(file='st_drugnet.txt',st_drugnet,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(file='st_drugname.txt',drug,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
########################################################
index=apply(as.matrix(drug),1,function(x){which(x==database)})
drugcate=line[(index-1)*5+2]
cate=str_split(drugcate,' ')
pt=unlist(cate)
dcate=as.matrix(unique(pt))
a=as.data.frame(table(pt))
catename=as.character(a[which(a[,2]>1),1])
###############st drug feature####################
st_drugfeature=matrix(0,length(drug),length(catename))
for (i in 1:length(drug)){
	name=cate[[i]]
	pl=apply(as.matrix(name),1,function(x){which(x==catename)})
	st_drugfeature[i,unlist(pl)]=1
}
write.table(file='st_drugfeature27.txt',st_drugfeature,row.names = FALSE,col.names = FALSE,quote = FALSE)

##############sd drug feature################

dindex=apply(as.matrix(alldrug),1,function(x){which(x==database)})

drugcate=line[(dindex-1)*5+2]
cate=str_split(drugcate,' ')
pt=unlist(cate)
dcate=as.matrix(unique(pt))
a=as.data.frame(table(pt))
catename=as.character(a[which(a[,2]>100),1])

sd_drugfeature=matrix(0,length(alldrug),length(catename))
for (i in 1:length(alldrug)){
	name=cate[[i]]
	pl=apply(as.matrix(name),1,function(x){which(x==catename)})
	sd_drugfeature[i,unlist(pl)]=1
}
write.table(file='sd_drugfeature46.txt',sd_drugfeature,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(file='sd_drug_name.txt',alldrug,row.names = FALSE,col.names = FALSE,quote = FALSE)

######################sd_label#####################
sd_label=matrix(0,length(alldrug),length(target))
drug=t[,1]
protein=t[,2]
for (i in 1:length(alldrug)){
		p=which(alldrug[i]==drug)
		if (length(p)!=0){
			tt=protein[p]
			position=apply(as.matrix(tt),1,function(x){which(target==x)})
			sd_label[i,position]=1
			}
}
write.table(file='sd_label.txt',sd_label,row.names = FALSE,col.names = FALSE,quote = FALSE)


