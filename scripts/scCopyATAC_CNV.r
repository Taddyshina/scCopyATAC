library(pheatmap)
library("dendextend")
dim(sum_logfc_A7[,which(min_w>80)])
new_A7<-sum_logfc_A7[,c(-752,-753,-754,-755,-603,-604)]
sum_logfc_A7<-new_A7
ann_colors = list(
  ward.D2=c("1"="#000000",
        "2"="#AE433D",
        "3"="#97FFFF",
        "4"="#354E87"),
  subgroup=c("G3"="#F1DF82",
        "G34"="#EE3B3B",
        "G4"="#008500",
        "immune"="#354E87",
        "proj_like"="#66C2A5",
        "SHH"="#999999",
        "unkn"="#8DA0CB",
        "WNT"="#FFFFB3"))


hclust_avg<-dist(sum_logfc_A7)
tree<-hclust(hclust_avg,method="ward.D2")
cut<- cutree(tree,k=6,order_clusters_as_data=FALSE)
c<-as.data.frame(table(cut))
test_A7<-sum_logfc_A7[which(cut %in% c[which(c$Freq>=100),]$cut),]
dim(test_A7)
ward.D2_clusters_A7<-cut[which(cut %in% c[which(c$Freq>=100),]$cut)]
ward.D2_clusters_A7<-as.data.frame(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
sub_A7_cell_info<-A7_cell_info[rownames(test_A7),]
sub_A7_cell_info<-cbind(sub_A7_cell_info,ward.D2_clusters_A7$ward.D2_clusters_A7)
colnames(sub_A7_cell_info)<-c("Sample","Clusters","subgroup","ward.D2")
sub_A7_cell_info<-sub_A7_cell_info[order(sub_A7_cell_info$ward.D2),]
pheatmap(sum_logfc_A7[rownames(sub_cell_info_A7),],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 3,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = chromosome,
         main ='sub_870w_0.36 A7 cut 4 CNV',
         annotation_colors = ann_colors
)





hclust_avg<-dist(sum_logfc_A7)
tree<-hclust(hclust_avg,method="ward.D2")
plot(tree,hang=-0.01,cex=0.7,label=FALSE,main="A7 tree step1",ylab=NULL,xlab="cell barcodes")
abline(h=1000,col='blue')
abline(h=460,col='red')
abline(h=200,col='green')
abline(h=170,col='orange')
abline(h=130,col='navy')
abline(h=85,col='grey')
cut<- cutree(tree,h=1000,order_clusters_as_data=FALSE)
c<-as.data.frame(table(cut))
table(cut)
sub_A7_1<-sum_logfc_A7[which(cut==2),]
sub_cut_1<-cut[which(cut==2)]
#sub_cut_1[which(cut==2)]="one"
table(sub_cut_1)

cut<- cutree(tree,h=460,order_clusters_as_data=FALSE)
c<-as.data.frame(table(cut))
table(cut)
sub_A7_2<-sum_logfc_A7[which(cut==1),]
sub_cut_2<-cut[which(cut==1)]
table(sub_cut_2)

cut<- cutree(tree,h=200,order_clusters_as_data=FALSE)
c<-as.data.frame(table(cut))
table(cut)
sub_A7_3<-sum_logfc_A7[which(cut==26|cut==28),]
sub_cut_3<-cut[which(cut==26|cut==28)]
table(sub_cut_3)
cut<- cutree(tree,h=170,order_clusters_as_data=FALSE)
cut<- cutree(tree,h=130,order_clusters_as_data=FALSE)
cut<- cutree(tree,h=85,order_clusters_as_data=FALSE)
table(cut)
c<-as.data.frame(table(cut))
sub_A7_4<-sum_logfc_A7[which(cut %in% c[which(c$Freq>=40),]$cut & cut != 98&cut!=125),]
sub_cut_4<-cut[which(cut %in% c[which(c$Freq>=40),]$cut & cut != 98&cut!=125)]
table(sub_cut_4)

sub_cut_sum<-c(sub_cut_1,sub_cut_2,sub_cut_3,sub_cut_4)
sub_A7_sum<-rbind(sub_A7_1,sub_A7_2,sub_A7_3,sub_A7_4)

table(sub_cut_sum)
norm.mat.smooth<-t(sub_A7)
SDM <-NULL
SSD <-NULL
for(i in min(as.integer(sub_cell_info_A7$ward.D2)):max(as.integer(sub_cell_info_A7$ward.D2))){
  
  data.c <- apply(norm.mat.smooth[, which(sub_cell_info_A7$ward.D2==i)],1, median)
  sx <- max(c(0.05, 0.5*sd(data.c)))
  GM3 <- mixtools::normalmixEM(data.c, lambda = rep(1,3)/3, mu = c(-0.2, 0, 0.2), sigma = sx,arbvar=FALSE,ECM=FALSE,maxit=5000)
  SDM <- c(SDM, GM3$sigma[1])
  SSD <- c(SSD, sd(data.c))
  i <- i+1
}

sub_cluster2_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="2"),]),]
sub_cluster3_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="3"),]),]
sub_cluster4_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="4"),]),]



ward.D2_clusters_A7<-cut
ward.D2_clusters_A7[which(ward.D2_clusters_A7!=1&ward.D2_clusters_A7!=2&ward.D2_clusters_A7!=3)]<-4


ward.D2_clusters_A7<-as.data.frame(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
sub_cluster1_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="1"),]),]
sub_cluster2_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="2"),]),]
sub_cluster3_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="3"),]),]
sub_cluster4_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="4"),]),]

sub_cluster1_A7<-sub_cluster1_A7[sample(c(1:length(rownames((sub_cluster1_A7)))),40,replace = FALSE),]
sub_cluster2_A7<-sub_cluster2_A7[sample(c(1:length(rownames((sub_cluster2_A7)))),40,replace = FALSE),]
sub_cluster3_A7<-sub_cluster3_A7[sample(c(1:length(rownames((sub_cluster3_A7)))),40,replace = FALSE),]
sub_cluster4_A7<-sub_cluster4_A7[sample(c(1:length(rownames((sub_cluster4_A7)))),40,replace = FALSE),]
sub_A7<-rbind(sub_cluster3_A7,sub_cluster4_A7)
sub_cell_info_A7<-A7_cell_info[rownames(sub_A7),]
sub_cell_info_A7<-cbind(sub_cell_info_A7,ward.D2_clusters_A7[rownames(sub_A7),]$ward.D2_clusters_A7)
colnames(sub_cell_info_A7)<-c("Sample","Clusters","subgroup","ward.D2")

A7_barcode<-read.table("scATAC_A7.barcodeTranslate.tsv")
A7_barcode[,2]<-paste("A7",A7_barcode[,2],sep="#")
uniq_A7<-A7_barcode[!duplicated(A7_barcode$V2),]
rownames(uniq_A7)<-uniq_A7$V2
uniq_A7<-uniq_A7[rownames(cell_info)[which(cell_info$Sample=="A7")],]
barcode_A7<-uniq_A7$V1
A7_cell_info<-cell_info[rownames(uniq_A7),]
rownames(uniq_A7)<-uniq_A7$V1
uniq_A7<-uniq_A7[rownames(sum_logfc_A7),]
A7_cell_info<-cell_info[uniq_A7$V2,]
rownames(A7_cell_info)<-rownames(sum_logfc_A7)
sub_cell_info_A7


sub_cell_info_A7<-A7_cell_info[rownames(ward.D2_clusters_A7),]
sub_cell_info_A7<-cbind(A7_cell_info,as.character(ward.D2_clusters_A7))
colnames(sub_cell_info_A7)<-c("Sample","Clusters","subgroup","ward.D2")
chromosome<-sum_annotation_col_A6[c(-752,-753,-754,-755,-603,-604)]
colnames(chromosome)="Chr"
pheatmap(sum_logfc_A7[rownames(sub_cell_info_A7),],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 3,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = chromosome,
         main ='sub_870w_0.36 A7 cut 4 CNV',
         annotation_colors = ann_colors
)

pheatmap(sum_logfc_A7[rownames(sub_cell_info_A7),],
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 300,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = new_col_ann,
         main ='sub_870w_0.36 A7 cut 4 CNV',
         clustering_method="ward.D2",
         annotation_colors = ann_colors
)

pheatmap(sub_A7[rownames(sub_sub_cell_info_A7),],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = as.data.frame(sub_sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 3,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.4,
         breaks=seq(-2,2,0.04),
         annotation_col = new_col_ann,
         annotation_colors = ann_colors,
         main ='sub_870w_0.36 A7 sub cut 4 CNV'
         
)

library(vioplot)
vioplot(colMeans(sub_cluster1_A7), colMeans(sub_cluster2_A7),colMeans(sub_cluster3_A7),colMeans(sub_cluster4_A7)
        ,main = "Distribution of CNV signals on 804 windows", 
        col=c("#000000","#AE433D","#97FFFF","#354E87"), xlab="sample",ylab="window mean CNV signal") 
legend("topleft", legend=c("cluster1", "cluster2", "cluster3","cluster4"),
       fill=c("#000000","#AE433D","#97FFFF","#354E87"), cex = 0.8)



sub_A7<-rbind(sub_cluster2_A7,sub_cluster3_A7)[,which(min_w>80)]
sub_sub_cell_info_A7<-sub_cell_info_A7[rownames(sub_A7),]



sub_A7
hclust_avg<-dist(sub_A7)
tree<-hclust(hclust_avg,method="ward.D2")
plot(tree,hang=-0.01,cex=0.7,label=FALSE,main="A7 tree",ylab=NULL,xlab="cell barcodes")
abline(h=380,col='red')
abline(h=190,col='green')
abline(h=90,col='blue')
cut<- cutree(tree,h=380,order_clusters_as_data=FALSE)
table(cut)
cut1<- cutree(tree,h=190,order_clusters_as_data=FALSE)
table(cut1)

ward.D2_clusters_A7<-as.data.frame(cut1,rownames(cut1))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
colnames(ward.D2_clusters_A7)<-c("cut1","ward.D2")
sub_cluster1_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut1=="1"),]),]
sub_cluster3_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut1=="3"),]),]
sub_cluster4_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut1=="4"),]),]
sub_A7<-rbind(sub_cluster1_A7,sub_cluster3_A7,sub_cluster4_A7)

hclust_avg<-dist(sub_A7)
tree<-hclust(hclust_avg,method="ward.D2")
plot(tree,hang=-0.01,cex=0.7,label=FALSE,main="A7 tree",ylab=NULL,xlab="cell barcodes")
abline(h=170,col='green')
cut<- cutree(tree,h=170,order_clusters_as_data=FALSE)
table(cut)
abline(h=81,col='blue')
cut<- cutree(tree,h=81,order_clusters_as_data=FALSE)

ward.D2_clusters_A7<-as.data.frame(cut,rownames(cut))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
colnames(ward.D2_clusters_A7)<-c("cut","ward.D2")
sub_cluster1_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="1"),]),]
sub_cluster2_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="2"),]),]
sub_cluster3_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="3"),]),]
sub_cluster23_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="23"),]),]
sub_cluster24_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="24"),]),]
sub_cluster31_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="31"),]),]
sub_cluster32_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="32"),]),]
sub_cluster40_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="40"),]),]
sub_cluster41_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="41"),]),]
sub_cluster42_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="42"),]),]
sub_cluster50_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$cut=="50"),]),]

sub_A7<-rbind(sub_cluster1_A7,sub_cluster2_A7,sub_cluster3_A7,sub_cluster23_A7,sub_cluster24_A7,
              sub_cluster31_A7,sub_cluster32_A7,sub_cluster40_A7,sub_cluster41_A7,sub_cluster42_A7,sub_cluster50_A7)


sub_A7
hclust_avg<-dist(sub_A7)
tree<-hclust(hclust_avg,method="ward.D2")
plot(tree,hang=-0.01,cex=0.7,label=FALSE,main="A7 tree",ylab=NULL,xlab="cell barcodes")
abline(h=81,col='blue')
cut<- cutree(tree,h=81,order_clusters_as_data=FALSE)
table(cut)
ward.D2_clusters_A7<-as.data.frame(cut,rownames(cut))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
colnames(ward.D2_clusters_A7)<-c("cut1","ward.D2")
sub_cell_info_A7<-A7_cell_info[ward.D2_clusters_A7$ward.D2,]
sub_cell_info_A7<-cbind(sub_cell_info_A7,as.character(ward.D2_clusters_A7$cut1))
colnames(sub_cell_info_A7)<-c("Sample","Clusters","subgroup","ward.D2")

pheatmap(sum_logfc_A7[ward.D2_clusters_A7$ward.D2,],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 3,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = chromosome,
         main ='sub_870w_0.36 A7 cut 4 CNV',
         annotation_colors = ann_colors
)






ward.D2_clusters_A7<-as.data.frame(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
sub_cluster1_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="1"),]),]
sub_cluster2_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="2"),]),]
sub_cluster3_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="3"),]),]
sub_cluster4_A7<-sum_logfc_A7[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="4"),]),]




plot(tree,hang=-0.01,cex=0.7,label=FALSE,main="A7 tree",ylab=NULL,xlab="cell barcodes")
abline(h=130,col='red')
abline(h=70,col='green')
abline(h=90,col='blue')
cut<- cutree(tree,h=130,order_clusters_as_data=FALSE)
table(cut)
cut2<- cutree(tree,h=70,order_clusters_as_data=FALSE)
table(cut2)
cut3<- cutree(tree,h=90,order_clusters_as_data=FALSE)
table(cut3)

pheatmap(sub_A7[rownames(sub_sub_cell_info_A7),],
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_row = as.data.frame(sub_sub_cell_info_A7)[,c(1,3)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 300,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = new_col_ann,
         annotation_colors = ann_colors,
         main ='sub_809w A7 sub cut 3 CNV',
         clustering_method="ward.D2"
)

ward.D2_clusters_A7<-cut
ward.D2_clusters_A7[which(ward.D2_clusters_A7!=1&ward.D2_clusters_A7!=2&ward.D2_clusters_A7!=3)]<-4
ward.D2_clusters_A7[which(cut2==1)]<-"cut1/1"
ward.D2_clusters_A7[which(cut2==2)]<-"cut1/2"
ward.D2_clusters_A7[which(cut2==3)]<-"cut1/3"
ward.D2_clusters_A7[which(cut3==9)]<-"cut3/1"
ward.D2_clusters_A7[which(cut3==10)]<-"cut3/2"
ward.D2_clusters_A7[which(cut3==11)]<-"cut3/3"
ward.D2_clusters_A7[which(cut3==12)]<-"cut3/4"
ward.D2_clusters_A7[which(cut3==13)]<-"cut3/5"
table(ward.D2_clusters_A7)
head(ward.D2_clusters_A7)
ward.D2_clusters_A7<-as.data.frame(ward.D2_clusters_A7)
ward.D2_clusters_A7<-cbind(ward.D2_clusters_A7,rownames(ward.D2_clusters_A7))
sub_cell_info_A7<-A7_cell_info[rownames(ward.D2_clusters_A7),]
sub_cell_info_A7<-cbind(sub_cell_info_A7,as.character(ward.D2_clusters_A7$ward.D2_clusters_A7))
colnames(sub_cell_info_A7)<-c("Sample","Clusters","subgroup","ward.D2")
new_col_ann<-as.data.frame(sum_annotation_col_A7[which(min_w>80),])
rownames(new_col_ann)<-colnames(sum_logfc_A7[rownames(sub_cell_info_A7),])
colnames(new_col_ann)<-c("chr")


new_ann_colors = list(
  ward.D2=c("cut1/1"="#FFD92F",
            "cut1/2"="#FB8072",
            "cut1/3"="#8B7D68",
            "cut3/1"="#458B74",
            "cut3/2"="#8B312F",
            "cut3/3"="#AAED0E",
            "cut3/4"="#8EE5EE",
            "cut3/5"="#8B0000",
            "2"="#AE433D",
            "4"="#354E87"),
  subgroup=c("G3"="#F1DF82",
             "G34"="#EE3B3B",
             "G4"="#008500",
             "immune"="#354E87",
             "proj_like"="#66C2A5",
             "SHH"="#999999",
             "unkn"="#8DA0CB",
             "WNT"="#FFFFB3"))

pheatmap(sub_A7,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(3,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 300,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = new_col_ann,
         annotation_colors = new_ann_colors,
         main ='sub_809w A7 sub cut 4 CNV',
         clustering_method="ward.D2"
)

pheatmap(sub_A7,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_row = as.data.frame(sub_cell_info_A7)[,c(1,4)],
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 300,
         color = mycolors,
         cellwidth = 0.4,
         cellheight = 0.15,
         breaks=seq(-2,2,0.04),
         annotation_col = new_col_ann,
         annotation_colors = new_ann_colors,
         main ='sub_809w A7 sub cut 4 CNV',
         clustering_method="ward.D2"
)



sub_cluster1_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut1/1"),]),]
sub_cluster2_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut1/2"),]),]
sub_cluster3_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut1/3"),]),]
sub_cluster4_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut3/1"),]),]
sub_cluster5_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut3/2"),]),]
sub_cluster6_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut3/3"),]),]
sub_cluster7_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut3/4"),]),]
sub_cluster8_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="cut3/3"),]),]
sub_cluster9_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="2"),]),]
sub_cluster10_A7<-rbind(sub_logfc_A7,sub_logfc_A72)[rownames(ward.D2_clusters_A7[which(ward.D2_clusters_A7$ward.D2_clusters_A7=="4"),]),]

vioplot(colMeans(sub_cluster1_A7), colMeans(sub_cluster2_A7),colMeans(sub_cluster3_A7),colMeans(sub_cluster9_A7)
        ,colMeans(sub_cluster4_A7), colMeans(sub_cluster5_A7),colMeans(sub_cluster6_A7),colMeans(sub_cluster7_A7)
        ,colMeans(sub_cluster8_A7),colMeans(sub_cluster10_A7)
        ,main = "Distribution of CNV signals on 804 windows", 
        col=c("#FFD92F","#FB8072","#8B7D68","#AE433D","#458B74","#8B312F","#AAED0E","#8EE5EE","#8B0000","#354E87"), xlab="sample",ylab="window mean CNV signal") 
legend("bottomleft", legend=c("cut1/1", "cut1/2", "cut1/3","2","cut3/1","cut3/2","cut3/3","cut3/4","cut3/5","4"),
       fill=c("#FFD92F","#FB8072","#8B7D68","#AE433D","#458B74","#8B312F","#AAED0E","#8EE5EE","#8B0000","#354E87"), cex = 0.8)
abline(h=0,col='red')

library(ggplot2)
library(viridis)
library(hrbrthemes)
data <- as.data.frame(sub_cell_info_A7)[,c(3,4)]
data<-cbind(data,1)
# Stacked + percent
ggplot(data, aes(fill=subgroup, x=ward.D2,y="1")) + 
  geom_bar(position="fill", stat="identity")


ggplot(data, aes(fill=subgroup, y="1", x=ward.D2)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T) +
  ggtitle("ATAC subgroup distrubute in tree clusters ") +
  xlab("")+
  ylab("")

colors<-c(rep("#FFD92F",sum(sub_cell_info_A7$ward.D2=="cut1/1")),
          rep("#FB8072",sum(sub_cell_info_A7$ward.D2=="cut1/2")),
          rep("#8B7D68",sum(sub_cell_info_A7$ward.D2=="cut1/3")),
          rep("#AE433D",sum(sub_cell_info_A7$ward.D2=="2")),
          rep("#458B74",sum(sub_cell_info_A7$ward.D2=="cut3/1")),
          rep("#8B312F",sum(sub_cell_info_A7$ward.D2=="cut3/2")),
          rep("#AAED0E",sum(sub_cell_info_A7$ward.D2=="cut3/3")),
          rep("#8EE5EE",sum(sub_cell_info_A7$ward.D2=="cut3/4")),
          rep("#8B0000",sum(sub_cell_info_A7$ward.D2=="cut3/5")),
          rep("#354E87",sum(sub_cell_info_A7$ward.D2=="4"))
          )

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)
setwd("C:/Users/Ted/Desktop/infer_scCNV_from_ATAC/A7/")
counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc) <- annotations


pbmc
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments


pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc<-pbmc[,rownames(sub_cell_info_A7)]

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc$clusters<-sub_cell_info_A7[colnames(pbmc),]$ward.D2
  
  
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc,pt.size =2,group.by="clusters",cols=c("cut1/1"="#FFD92F",
                                                            "cut1/2"="#FB8072",
                                                            "cut1/3"="#8B7D68",
                                                            "cut3/1"="#458B74",
                                                            "cut3/2"="#8B312F",
                                                            "cut3/3"="#AAED0E",
                                                            "cut3/4"="#8EE5EE",
                                                            "cut3/5"="#8B0000",
                                                            "2"="#AE433D",
                                                            "4"="#354E87")) 

DimPlot(object = pbmc,pt.size =2,label = TRUE)




library(ggplot2)
library(viridis)
library(hrbrthemes)
data <- cbind(as.character(pbmc$seurat_clusters),sub_cell_info_A7[colnames(pbmc),]$Clusters)
data<-as.data.frame(cbind(data,1))
colnames(data)<-c("clusters","seurat","1")
# Stacked + percent
ggplot(data, aes(fill=clusters, x=seurat,y="1")) + 
  geom_bar(position="fill", stat="identity")


ggplot(data, aes(fill=clusters, y="1", x=seurat)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T) +
  ggtitle("sum subgroup distrubute in  A7 clusters ") +
  xlab("")+
  ylab("")


data <- cbind(as.character(pbmc$seurat_clusters),sub_cell_info_A7[colnames(pbmc),]$ward.D2)
data<-as.data.frame(cbind(data,1))
colnames(data)<-c("clusters","seurat","1")
# Stacked + percent
ggplot(data, aes(fill=clusters, x=seurat,y="1")) + 
  geom_bar(position="fill", stat="identity")


ggplot(data, aes(fill=clusters, y="1", x=seurat)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T) +
  ggtitle("A7 clusters distrubute in tree clusters ") +
  xlab("")+
  ylab("")

ggplot(data, aes(x=clusters, y="1", fill=seurat)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T) +
  ggtitle("tree clusters distrubute in A7 clusters ") +
  xlab("")+
  ylab("")


gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
DefaultAssay(pbmc) <- 'RNA'

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

FeaturePlot(pbmc, features = c("PCBP3", "EDEM2", "NPAS4", "NR4A1", "CTD-3088G3.8", "SIPA1L3"),pt.size = 2)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(1, 2,9), min.pct = 0.25)
head(cluster5.markers, n = 5)

FeaturePlot(pbmc, features = c("TTYH1", "TST", "LCNL1", "RGMA", "XYLT1"),pt.size = 1.5)

CNV_pbmc <- CreateSeuratObject(counts = t(sub_A7),
                               assay = "peaks",
                               meta.data = as.data.frame(sub_cell_info_A7))
CNV_pbmc <- RunTFIDF(CNV_pbmc)
CNV_pbmc <- FindTopFeatures(CNV_pbmc, min.cutoff = 'q0')
CNV_pbmc<- RunSVD(CNV_pbmc)
CNV_pbmc <- RunUMAP(object = CNV_pbmc,reduction = 'lsi',dims = 1:50)**
CNV_pbmc <- FindNeighbors(object = CNV_pbmc, reduction = 'lsi', dims = 1:50)
CNV_pbmc <- FindClusters(object = CNV_pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = CNV_pbmc, label = TRUE) + NoLegend()
CNV_pbmc$clusters<-sub_cell_info_A7[colnames(CNV_pbmc),]$ward.D2
DimPlot(object = CNV_pbmc,pt.size =2,group.by="clusters",cols=c("cut1/1"="#FFD92F",
                                                            "cut1/2"="#FB8072",
                                                            "cut1/3"="#8B7D68",
                                                            "cut3/1"="#458B74",
                                                            "cut3/2"="#8B312F",
                                                            "cut3/3"="#AAED0E",
                                                            "cut3/4"="#8EE5EE",
                                                            "cut3/5"="#8B0000",
                                                            "2"="#AE433D",
                                                            "4"="#354E87")) 
CNV_pbmc<-CNV_pbmc[,colnames(gene.activities)]
CNV_pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
CNV_pbmc <- NormalizeData(
  object = CNV_pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(CNV_pbmc$nCount_RNA)
)
DefaultAssay(CNV_pbmc) <- 'RNA'
cluster5.markers <- FindMarkers(CNV_pbmc, ident.1 = 'cut3/4', min.pct = 0.25,group.by = "clusters")
head(cluster5.markers, n = 6)
FeaturePlot(CNV_pbmc, features = c("TRIM55", "PPP2R2B", "ZNRF3", "DCUN1D2", "PPP2R2C"),pt.size = 2)
FeaturePlot(CNV_pbmc, features = c( "PPP2R2B", "TRIM5", "CD163", "CSF1R", "CDK5","IDO1","MICA","ULBP2"),pt.size = 2)
write.csv(rownames(cluster5.markers),"genelists.txt",quote = F,row.names = F)


library(celldex)
ref <- HumanPrimaryCellAtlasData()
library(SingleR)
comomgenes<-intersect(rownames(ref),rownames(gene.activities))
genes<-gene.activities[comomgenes,]
com_ref<-ref[comomgenes,]
pred2 <- SingleR(test=genes, ref=com_ref, clusters = CNV_pbmc$clusters, labels=com_ref$label.main)
table(pred2$labels)
plotScoreHeatmap(pred2)


