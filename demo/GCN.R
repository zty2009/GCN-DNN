
GCN<-function(netedge,pnetedge,label,allfeature){
A=matrix(0,(nrow(netedge)*nrow(pnetedge)),(nrow(netedge)*nrow(pnetedge)))
for (i in 1:nrow(netedge)){
	for (j in 1:nrow(pnetedge)){
	now=matrix(0,nrow(netedge),nrow(pnetedge))
		index1=which(label[i,]==1)
		now[i,index1]=1
		index2=which(label[,j]==1)
		now[index2,j]=1
		index3=which(netedge[i,]==1)
		if (length(index3)!=0){
		index3.p=apply(as.matrix(index3),1,function(x){which(label[x,]==1)})
		for (k in 1:length(index3)){
			now[index3[k],index3.p[[k]]]=0.5
			}
			}
		index4=which(pnetedge[i,]>0)
		if (length(index4)!=0){
		index4.p=apply(as.matrix(index4),1,function(x){which(label[,x]==1)})
		for (k in 1:length(index4)){
			now[index4.p[[k]],index4[k]]=0.5
			}
			}
		now[i,j]=1
		A[((i-1)*nrow(pnetedge)+j),]=as.vector(now)
		}	
}
D=colSums(A)
D_hat=diag(D^(-0.5))
a=relu(D_hat%*%A%*%D_hat%*%allfeature)
return(a)
}
