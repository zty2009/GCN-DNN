setwd('C:/Users/zty20/Documents/R/drug target/demo')
library(data.table)
library(keras)
library(sigmoid)
source('GCN.R')
############st##################
netedge=read.table('sd_drugnet.txt')
X=as.matrix(read.table('sd_drugfeature46.txt'))
pnetedge=read.table('sd_proteinnet.txt')
pX=as.matrix(read.table('sd_proteinfeature.txt'))

original=read.table('demodata.txt',stringsAsFactor=F)
alldrug=read.table('sd_drug_name.txt',stringsAsFactor=F)
allprotein=read.table('sd_proteinname.txt',stringsAsFactor=F)
drug_index=apply(as.matrix(original[,1]),1,function(x){which(x==alldrug)})
newdrug_index=setdiff(1:nrow(alldrug),drug_index)
newdrug=alldrug[newdrug_index,]




label=matrix(0,nrow(alldrug),nrow(allprotein))

for (i in 1:nrow(alldrug)){
		index1=which(original[,1]==alldrug[i,1])
		index2=apply(as.matrix(original[index1,2]),1,function(x){which(allprotein==x)})
		label[i,index2]=1
		}
allfeature=matrix(0,nrow(pX)*nrow(X),ncol(X)+ncol(pX))
for (i in 1:nrow(X)){
	for (j in 1:nrow(pX)){
		allfeature[((i-1)*nrow(pX)+j),]=c(pX[j,],X[i,])
		}
}
feature_gcn=GCN(netedge,pnetedge,label,allfeature)


label=as.vector(label)
label1=matrix(0,length(label),3)
for (i in 1:nrow(allprotein)){
	for (j in 1:nrow(alldrug)){
			label1[(i-1)*nrow(alldrug)+j,1]=alldrug[j,1]
			label1[(i-1)*nrow(alldrug)+j,2]=allprotein[i,1]
			label1[(i-1)*nrow(alldrug)+j,3]=label[(i-1)*nrow(alldrug)+j]
			}
}
index=which(label1[,1]==newdrug)
index=apply(as.matrix(newdrug),1,function(x){which(label1[,1]==x)})
posi=which(label==1)
negi=as.vector(index)
neg=feature_gcn[negi,]
pos=feature_gcn[posi,]
set.seed(100)
samplepos=sample(1:nrow(pos))
set.seed(100)
sampletest=sample(1:nrow(neg))
library(pROC)
library(PRROC)

posindex=cut(1:nrow(pos),10, labels=F)
negindex=cut(1:nrow(neg),10, labels=F)
for (ll in  1:10){
	cv.test1=samplepos[which(posindex==ll)]

	cv.test2=sampletest[which(negindex==ll)]

	cv.sample=rbind(pos[setdiff(samplepos,cv.test1),],neg[setdiff(sampletest,cv.test2),])
	cv.test=rbind(pos[cv.test1,],neg[cv.test2,])
	testindex=c(posi[cv.test1],negi[cv.test2])
	samplelabel=matrix(0,nrow(cv.sample),1)
	samplelabel[1:nrow(pos[setdiff(samplepos,cv.test1),])]=1
	testlabel=matrix(0,nrow(cv.test),1)
	testlabel[1:length(cv.test1)]=1
	samplelabel1=to_categorical(samplelabel,2)
	result=matrix(0,length(testlabel),5)
	for (kk in 1:5){ 
	model <- keras_model_sequential() 
	model %>%
	

  layer_dense(units = 256, activation = 'relu') %>% 
  layer_dropout(rate = 0.3) %>%   
  layer_dense(units = 128, activation = 'relu') %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units = 2, activation = 'sigmoid')
  
  model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

model %>% fit(
 as.matrix(cv.sample), samplelabel1, 
  epochs = 50, batch_size = 3, 
  validation_split = 0.1
)
result[,kk]=predict(model,as.matrix(cv.test))[,2]
}
print(ll)
plp=data.frame(label1[testindex,],rowSums(result)/5)
write.table(file='Sd_result.txt',plp,append=T,row.names = FALSE,col.names = FALSE,quote = FALSE)
}
plp=read.table('Sd_result.txt')
roc(plp[,3],plp[,4])
pr.curve(plp[plp[,3]==1,4],plp[plp[,3]==0,4])