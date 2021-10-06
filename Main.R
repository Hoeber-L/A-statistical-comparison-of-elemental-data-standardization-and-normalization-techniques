#Author: Lars Höber
#Last updated: 04.10.21


library("readr")
library("here")
library("ggplot2")
library("ggplotify")
library("ggrepel")
library("gridExtra")
library("ggcorrplot")
library("ggfortify")
library("psych")
library("GGally")

#--------------------------------------------------------------------------------
#-------------------------data formatting---------------------------------------------
temp="[]"
temp=temp[-1]

#compare for relevant data/ return irrelevant column indices
comparecsv=function(data,temp) (
  for(i in 1:ncol(data))
  {
    if(data[1,i]=="d13C" | data[1,i]=="d18O"| 
       data[1,i]=="Ca"       |
       data[1,i]=="Mg"| data[1,i]=="Fe"| data[1,i]=="Mn"| data[1,i]=="Sr")
    {}
    else
    {
      temp=append(temp,i)
    }
    if(i==ncol(data))
    {
      return(temp)
    }
  }
)

#trimming data/ remove irrelevant columns
remex=function(data)(
  for(i in length(as.double(comparecsv(data,temp))):1)
  {
    data=data[-as.double(comparecsv(data,temp))[i]]
    if(i==1)
    {
      colnames(data)=data[1,]
      data=data[-1,]
      return(data)
    }
  }
  
)

#creating output matrix filled with zeroes
matex0=function(data)(
  return(matrix(0,nrow(data),ncol(data)))
)

#repalcing zeroes with data
mx0=function(data,mat)(
  for(i in 1:ncol(data))
  {
    for(j in 1:nrow(data))
    {
      mat[j,i]=as.double(data[j,i])
      if(i==ncol(data))
      {
        if(j==nrow(data))
        {
          colnames(mat)=colnames(data)
          return(mat)
        }
      }
    }
  }
)

#everything in one command
mx=function(data)(
  mx0(remex(data),matex0(remex(data)))
)

#shortcut for standard path
loads=function(data)
{
  mx(read.csv(here("data",data), sep=";",header = FALSE))
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Normalisierung CA
corm1=function(x)(
  return((x)/(x[,3]))
)
c1x0=function(data,mat){
  for (k in 1:2) {
    mat[,k]=(data)[,k]
    if(k==ncol(data))
    {
      colnames(mat)=colnames(data)
    }
    
  }
  for(i in 3:ncol(data))
  {
    mat[,i]=corm1(data)[,i]
    if(i==ncol(data))
    {
      colnames(mat)=colnames(data)
      return(mat)
    }
  }
}
#Normalisierung einer Matrix durch Ca
c1x=function(data)(
  c1x0(data,matex0(data))
)
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#-Normalisierung 0-1
#function die einzelne spalten auf 0-1 normiert
norm1=function(x)(
  return((x-min(x))/(max(x)-min(x)))
)

n1x0=function(data,mat)(
  for(i in 1:ncol(data))
  {
    mat[,i]=norm1(data[,i])
    if(i==ncol(data))
    {
      colnames(mat)=colnames(data)
      return(mat)
    }
  }
)
#normalisierung einer Matrix von 0-1
n1x=function(data)(
  n1x0(data,matex0(data))
)
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#normalisierung Standard-Score
morm1=function(x)(
  return((x-mean(x))/(sd(x)))
)
m1x0=function(data,mat)(
  for(i in 1:ncol(data))
  {
    mat[,i]=morm1(data[,i])
    if(i==ncol(data))
    {
      colnames(mat)=colnames(data)
      return(mat)
    }
  }
)
#normalisierung einer Matrix Standart-Score
m1x=function(data)(
  m1x0(data,matex0(data))
)
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#formatting

fmt_dcimals = function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Principal component analysis

#Scree Plot scaled
pcascree0=function(data)(
  barplot(round(prcomp((data), scale=TRUE)$sdev^2/sum(prcomp((data), scale=TRUE)$sdev^2)*100, 1),
          main="Scree Plot",
          xlab="Principal Component", ylab="Percent Variation")
)

#Scree Plot unscaled
pcascree1=function(data)(
  barplot(prcomp((data), scale=FALSE)$sdev^2/sum(prcomp((data), scale=FALSE)$sdev^2)*100,
          
          main="Scree Plot",
          
          xlab="Principal Component", ylab="Percent Variation")
)
#biplot scaled 
bip0=function(data)(biplot(prcomp((data), scale=TRUE), scale = 0) )
#biplot unscaled 
bip1=function(data)(biplot(prcomp((data), scale=FALSE)) )

#PCA
PPCA=function(data,decimals=0)
{
  #theme_set(theme_classic())

  raw    =ggplot(cbind(loads(data)[,-3],prcomp(loads(data)[,-3],scale = T)$x[,1:2]), #generating raw data plot
                 aes(PC1,PC2))+
 
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+ #drawing confidence ellipse 
 
    geom_segment(data = as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L:2L]* #drawing Loadingarrows 
                   min(max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L])),
                   max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 2L])))
                   *0.8,
                 mapping = ggplot2::aes_string(x = 0,
                                               y = 0,
                                               xend = colnames(as.data.frame(prcomp(loads(data)[,-3])$rotation[,1:5]))[1L],
                                               yend = colnames(as.data.frame(prcomp(loads(data)[,-3])$rotation[,1:5]))[2L]),
                 arrow = grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+

    geom_label_repel(data = as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L:2L]* #labeling loadings
                  min(max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L])),
                  max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 2L])))
                  *0.8,
                  mapping = aes(x=(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L])
                                  *min(max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                  max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 2L])))
                                  *0.8,
                                y = (as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 2L])
                                  *min(max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                  max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(loads(data)[,-3],scale=T)$rotation[,1:5])[, 2L])))
                                  *0.8),
                  label=rownames(prcomp(loads(data)[,-3],scale=T)$rotation))+

    geom_point(shape=16, col="black", alpha=0.7)+ #adding datapoints
    
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+

    labs(x=paste("PC 1 (",round(summary(prcomp(loads(data)[,-3],scale=T ))$importance[2]*100, digits=1),"%)",sep=""),#adding graph labels
         y=paste("PC 2 (",round(summary(prcomp(loads(data)[,-3],scale=T ))$importance[5]*100, digits=1),"%)",sep=""),
         title= "PCA-biplot initial data",
         tag="A")
  
 #repeat with other data data 
  calc   =ggplot(cbind(c1x(loads(data))[,-3],
                       prcomp(c1x(loads(data))[,-3],
                              scale = T)$x[,1:2]),
                 aes(PC1,PC2))+
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+
    geom_segment(data = as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L:2L]*
                   min(max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                       max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                   0.8,
                 mapping = ggplot2::aes_string(x = 0, y = 0,xend = colnames(as.data.frame(prcomp(c1x(loads(data))[,-3])$rotation[,1:5]))[1L],yend = colnames(as.data.frame(prcomp(c1x(loads(data))[,-3])$rotation[,1:5]))[2L]),
                 arrow = grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    geom_label_repel(data = as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L:2L]*
                       min(max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                           max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))
                     *0.8,
                     mapping = aes(x=(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])
                                   *min(max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8,
                                   y = (as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])
                                   *min(max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(c1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8),
                     label=rownames(prcomp(c1x(loads(data))[,-3],scale=T)$rotation))+
    
    geom_point(shape=16, col="black", alpha=0.7) +
    
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    
    labs(x=paste("PC 1 (",round(summary(prcomp(c1x(loads(data))[,-3],scale=T))$importance[2]*100, digits=1),"%)",sep=""),
                                                        y=paste("PC 2 (",round(summary(prcomp(c1x(loads(data))[,-3],scale=T))$importance[5]*100, digits=1),"%)",sep=""),
                                                        title= "PCA-biplot Ca-normalized data", tag="B")
  
  
  minma  =ggplot(cbind(n1x(loads(data))[,-3],
                       prcomp(n1x(loads(data))[,-3],
                              scale = T)$x[,1:2]),
                 aes(PC1,PC2))+
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+
    geom_segment(data = as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L:2L]*
                   min(max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                       max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                   0.8,
                 mapping = ggplot2::aes_string(x = 0, y = 0,xend = colnames(as.data.frame(prcomp(n1x(loads(data))[,-3])$rotation[,1:5]))[1L],yend = colnames(as.data.frame(prcomp(n1x(loads(data))[,-3])$rotation[,1:5]))[2L]),
                 arrow = grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    geom_label_repel(data = as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L:2L]*
                       min(max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                           max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))
                     *0.8,
                     mapping = aes(x=(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])
                                   *min(max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8,
                                   y = (as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])
                                   *min(max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(n1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8),
                     label=rownames(prcomp(n1x(loads(data))[,-3],scale=T)$rotation))+
    
    geom_point(shape=16, col="black", alpha=0.7) +
    
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    
    labs(x=paste("PC 1 (",round(summary(prcomp(n1x(loads(data))[,-3],scale=T))$importance[2]*100, digits=1),"%)",sep=""),
                                                        y=paste("PC 2 (",round(summary(prcomp(n1x(loads(data))[,-3],scale=T))$importance[5]*100, digits=1),"%)",sep=""),
                                                        title= "PCA-biplot Min_Max data",
                                                        tag="C")
  
  sscore =ggplot(cbind(m1x(loads(data))[,-3],
                       prcomp(m1x(loads(data))[,-3],
                              scale = T)$x[,1:2]),
                 aes(PC1,PC2))+
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+
    geom_segment(data = as.data.frame(prcomp(m1x(loads(data))[,-3])$rotation[,1:5])[, 1L:2L]*
                   min(max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3])$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3])$rotation[,1:5])[, 1L])),
                       max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3])$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3])$rotation[,1:5])[, 2L])))*
                   0.8,
                 mapping = ggplot2::aes_string(x = 0, y = 0,xend = colnames(as.data.frame(prcomp(m1x(loads(data))[,-3])$rotation[,1:5]))[1L],yend = colnames(as.data.frame(prcomp(m1x(loads(data))[,-3])$rotation[,1:5]))[2L]),
                 arrow = grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    geom_label_repel(data = as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L:2L]*
                       min(max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                           max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))
                     *0.8,
                     mapping = aes(x=(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])
                                   *min(max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8,
                                   y = (as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])
                                   *min(max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 1L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 1L])),
                                        max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$x[,1:2])[, 2L])) /max(abs(as.data.frame(prcomp(m1x(loads(data))[,-3],scale=T)$rotation[,1:5])[, 2L])))*
                                     0.8),
                     label=rownames(prcomp(m1x(loads(data))[,-3],scale=T)$rotation))+
    
    geom_point(shape=16, col="black", alpha=0.7) +
    
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    
    labs(x=paste("PC 1 (",round(summary(prcomp(m1x(loads(data))[,-3],scale=T))$importance[2]*100, digits=1),"%)",sep=""),
         y=paste("PC 2 (",round(summary(prcomp(m1x(loads(data))[,-3],scale=T))$importance[5]*100, digits=1),"%)",sep=""),
         title= "PCA-biplot z-score data",
         tag="D")
  
  ggsave(paste(data,"_PCA_biplot_scaled.jpg",sep=""),arrangeGrob(raw,calc,minma,sscore,ncol=2),width = 8, height = 8,dpi=600)  
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Factoanalysis
FA=function(data,decimals=0)
{
  #generating raw data plot
  raw=ggplot(cbind(loads(data)[,-3],factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores), 
             aes(Factor1,Factor2))+
    #drawing confidence ellipse 
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+ 
    #drawing arrows
    geom_segment(
      data =    factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 2L]))))
      *0.8,
      mapping = ggplot2::aes_string(
        x = 0,
        y = 0,
        xend = rownames(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[1,1])),
        yend = rownames(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[1,2]))),
      arrow =   grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    #labeling
    geom_label_repel(
      data =    factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 2L]))))
      *0.8,
      mapping = aes(
        x=factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L]*  
          min(max(abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 2L]))))
        *0.8,
        y = factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 2L]* 
          min(max(abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(loads(data)[,-3],factanal(loads(data)[,-3],scale=T,2,method="varimax"))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings[, 2L]))))
        *0.8),
      label=rownames(prcomp(loads(data)[,-3],scale=T)$rotation))+
    #adding datapoints   
    geom_point(shape=16, col="black", alpha=0.7)+ 
    #adjusting decimal places
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    #generating axislabels with Variance information    
    labs(x=paste("Factor 1 (",round(colSums(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings^2)/nrow(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings),3)[1]*100," %)",sep=""),#adding graph labels
         y=paste("Factor 2 (",round(colSums(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings^2)/nrow(factanal(loads(data)[,-3],scale=T,2,method="varimax")$loadings),3)[2]*100," %)",sep=""),
         title= "FA-biplot initial data",
         tag="A")
  
  #generating ca normalized data plot
  calc=ggplot(cbind(c1x(loads(data))[,-3],factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores), 
              aes(Factor1,Factor2))+
    #drawing confidence ellipse 
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+ 
    #drawing arrows
    geom_segment(
      data =    factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = ggplot2::aes_string(
        x = 0,
        y = 0,
        xend = rownames(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[1,1])),
        yend = rownames(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[1,2]))),
      arrow =   grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    #labeling
    geom_label_repel(
      data =    factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = aes(
        x=factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]*  
          min(max(abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8,
        y = factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]* 
          min(max(abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(c1x(loads(data))[,-3],factanal(c1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(c1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8),
      label=rownames(prcomp(c1x(loads(data))[,-3],scale=T)$rotation))+
    #adding datapoints   
    geom_point(shape=16, col="black", alpha=0.7)+ 
    #adjusting decimal places
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    #generating axislabels with Variance information    
    labs(x=paste("Factor 1 (",round(colSums(factanal(c1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(c1x(loads(data))[,-3],2)$loadings),3)[1]*100," %)",sep=""),#adding graph labels
         y=paste("Factor 2 (",round(colSums(factanal(c1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(c1x(loads(data))[,-3],2)$loadings),3)[2]*100," %)",sep=""),
         title= "FA-biplot Ca-normalized data",
         tag="B")
  
  #generating Min-Max data plot
  minma=ggplot(cbind(n1x(loads(data))[,-3],factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores), 
               aes(Factor1,Factor2))+
    #drawing confidence ellipse 
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+ 
    #drawing arrows
    geom_segment(
      data =    factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = ggplot2::aes_string(
        x = 0,
        y = 0,
        xend = rownames(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[1,1])),
        yend = rownames(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[1,2]))),
      arrow =   grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    #labeling
    geom_label_repel(
      data =    factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = aes(
        x=factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]*  
          min(max(abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8,
        y = factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]* 
          min(max(abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(n1x(loads(data))[,-3],factanal(n1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(n1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8),
      label=rownames(prcomp(n1x(loads(data))[,-3],scale=T)$rotation))+
    #adding datapoints   
    geom_point(shape=16, col="black", alpha=0.7)+ 
    #adjusting decimal places
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    #generating axislabels with Variance information    
    labs(x=paste("Factor 1 (",round(colSums(factanal(n1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(n1x(loads(data))[,-3],2)$loadings),3)[1]*100," %)",sep=""),#adding graph labels
         y=paste("Factor 2 (",round(colSums(factanal(n1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(n1x(loads(data))[,-3],2)$loadings),3)[2]*100," %)",sep=""),
         title= "FA-biplot Min_Max data",
         tag="C")
  
  #generating z-score data plot
  sscore=ggplot(cbind(m1x(loads(data))[,-3],factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores), 
                aes(Factor1,Factor2))+
    #drawing confidence ellipse 
    stat_ellipse(geom="polygon", fill=5,col="black", alpha=0.5, level=0.95)+ 
    #drawing arrows
    geom_segment(
      data =    factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = ggplot2::aes_string(
        x = 0,
        y = 0,
        xend = rownames(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[1,1])),
        yend = rownames(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[1,2]))),
      arrow =   grid::arrow(length = grid::unit(8, 'points')), colour = "Red")+
    #labeling
    geom_label_repel(
      data =    factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L:2L]* 
        min(max(abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
            max(    abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
      *0.8,
      mapping = aes(
        x=factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]*  
          min(max(abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8,
        y = factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]* 
          min(max(abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 1L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 1L]))),
              max(    abs(as.data.frame(factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],2))$scores[,1:2])[, 2L])) /max(abs(as.data.frame(factanal(m1x(loads(data))[,-3],scale=T,2)$loadings[, 2L]))))
        *0.8),
      label=rownames(prcomp(m1x(loads(data))[,-3],scale=T)$rotation))+
    #adding datapoints   
    geom_point(shape=16, col="black", alpha=0.7)+ 
    #adjusting decimal places
    scale_x_continuous(labels = fmt_dcimals(1))+
    scale_y_continuous(labels = fmt_dcimals(1))+
    #generating axislabels with Variance information    
    labs(x=paste("Factor 1 (",round(colSums(factanal(m1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(m1x(loads(data))[,-3],2)$loadings),3)[1]*100," %)",sep=""),#adding graph labels
         y=paste("Factor 2 (",round(colSums(factanal(m1x(loads(data))[,-3],2)$loadings^2)/nrow(factanal(m1x(loads(data))[,-3],2)$loadings),3)[2]*100," %)",sep=""),
         title= "FA-biplot z-score data",
         tag="D") 
  ggsave(paste(data,"_FA_biplot_scaled.jpg",sep=""),arrangeGrob(raw,calc,minma,sscore,ncol=2),width = 8, height = 8,dpi=600)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Correlation Matrix ///http://www.sthda.com/english/wiki/correlation-matrix-an-r-function-to-do-all-you-need
Correlations=function(data)
{
  theme_set(theme_classic())
  raw    =ggcorrplot(cor(loads(data)[,-3]), hc.order = F, type = "lower",lab = TRUE)+
    labs(title= "Correlations initial-data",tag="A")
  calc   =ggcorrplot(cor(c1x(loads(data))[,-3]), hc.order = F, type = "lower",lab = TRUE)+
    labs(title= "Correlations calcium_normalized-data",tag="B")
  minma  =ggcorrplot(cor(n1x(loads(data))[,-3]), hc.order = F, type = "lower",lab = TRUE)+
    labs(title= "Correlations Min_Max-data",tag="C")
  sscore =ggcorrplot(cor(m1x(loads(data))[,-3]), hc.order = F, type = "lower",lab = TRUE)+
    labs(title= "Correlations standard_score-data",tag="D")
  ggsave(paste(data,"_Correlations.jpg",sep=""),arrangeGrob(raw,calc,minma,sscore,ncol=2),width = 8, height = 8,dpi=300)  
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#measure of the difference -> Euclidean distance and ca range seems best

#ca value range
ma=function(data){
  dd=max(loads(data)[,3])-min(loads(data)[,3])
  return(dd)
}
#difference in FA
mc=function(data){
  d=abs(factor.scores(m1x(c1x(loads(data)))[,-3],factanal(m1x(c1x(loads(data)))[,-3],scale=T,2,method="varimax"))$scores-factor.scores(m1x(loads(data))[,-3],factanal(m1x(loads(data))[,-3],scale=T,2,method="varimax"))$scores)
  dd=1
  dd=dd[-1]
  for (i in 1:nrow(loads(data))) 
  {
    dd[i]=sqrt(d[i,1]^2+d[i,2])
  }
  dd=mean(dd)
  return(dd)
}

#for factor loadings
nn=function(data){
  d=abs(sqrt(factanal(m1x(c1x(loads(data)))[,-3],2,method="varimax")$loading[,]^2)-sqrt(factanal(m1x(loads(data))[,-3],2,method="varimax")$loading[,]^2))
  dd=1
  dd=dd[-1]
  for (i in 1:ncol(loads(data)[,-3])) 
  {
    dd[i]=sqrt(d[i,1]^2+d[i,2])
  }
  dd=mean(dd)
  return(dd)
}

#differerence in original data
md=function(data)
{
  d=abs(m1x(c1x(loads(data)))[,-3]-m1x(loads(data))[,-3])
  dd=1
  dd=dd[-1]
  for (i in 1:nrow(loads(data))) {
    dd[i]=sqrt(+d[i,3]^2+d[i,4]^2+d[i,5]^2+d[i,6]^2)
  }
  dd=mean(dd)
  return(dd)
}
#difference in pca
mmd=function(data)
{
  
  d=abs(prcomp(m1x(c1x(loads(data))[,-3]))$x-prcomp(m1x(loads(data))[,-3])$x)
  dd=1
  dd=dd[-1]
  for (i in 1:nrow(loads(data))) {
    dd[i]=sqrt(d[i,1]^2+d[i,2]^2+d[i,3]^2+d[i,4]^2+d[i,5]^2+d[i,6]^2)
  }
  dd=mean(dd)
  return(dd)
}
#comparing the measurements
mm=function(dstart,dend){
  #print(noquote(c("range","normal","PCA","FA-scores","FA-loadings")))
  z=matrix(0,10,5)
  for (i in dstart:dend)
  {
    a=paste("Dataset",i,".csv",sep="")
    for (j in 1:5) {
      if (j==1){z[i,j]=ma(a)}
      if (j==2){z[i,j]=md(a)}
      if (j==3){z[i,j]=mmd(a)}
      if (j==4){z[i,j]=mc(a)}
      if (j==5){z[i,j]=nn(a)}
    }
  }
  
  return(round(z,3))
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
  #generaiting matrix representation of data
  lower = function(data,clust,inp){
    
    if (inp == "A")
    {
      cpairs=
        {
          set.seed(1234)
          d=loads(data)[,-3]
          nvar=c(colnames(d),"Cluster")
          d=cbind(d,kmeans(d,clust)$cluster)
          colnames(d)=nvar
          
          
          ggpairs(d[,cbind("d13C","d18O","Mg","Fe","Mn","Sr")],aes(color="Cluster"),diag="blank",upper="blank")+
            geom_point(colour=c("Black","Gray75","Grey50")[kmeans(d,clust)$cluster])+
            labs(title=paste("Dataset",gsub("([Dataset.csv])", "", data)," initial \"raw\" data", sep=""),tag="A") 
        }
      
      cpairs$plots = cpairs$plots[-(1:cpairs$nrow)]
      cpairs$yAxisLabels = cpairs$yAxisLabels[-1]
      cpairs$nrow = cpairs$nrow -1
      
      cpairs$plots = cpairs$plots[-(seq(cpairs$ncol, length(cpairs$plots), by = cpairs$ncol))]
      cpairs$xAxisLabels = cpairs$xAxisLabels[-cpairs$ncol]
      cpairs$ncol = cpairs$ncol - 1
    }
    if (inp == "B")
    {
      cpairs=
        {
          set.seed(1234)
          d=c1x(loads(data))[,-3]
          nvar=c(colnames(d),"Cluster")
          d=cbind(d,kmeans(d,clust)$cluster)
          colnames(d)=nvar
          
          
          ggpairs(d[,cbind("d13C","d18O","Mg","Fe","Mn","Sr")],aes(color="Cluster"),diag="blank",upper="blank")+
            geom_point(colour=c("Black","Gray75","Grey50")[kmeans(d,clust)$cluster])+
            labs(title=paste("Dataset",gsub("([Dataset.csv])", "", data)," initial Ca-normalized data", sep=""),tag="B") 
        }
      
      cpairs$plots = cpairs$plots[-(1:cpairs$nrow)]
      cpairs$yAxisLabels = cpairs$yAxisLabels[-1]
      cpairs$nrow = cpairs$nrow -1
      
      cpairs$plots = cpairs$plots[-(seq(cpairs$ncol, length(cpairs$plots), by = cpairs$ncol))]
      cpairs$xAxisLabels = cpairs$xAxisLabels[-cpairs$ncol]
      cpairs$ncol = cpairs$ncol - 1
    }
    if (inp == "C")
    {
      cpairs=
        {
          set.seed(1234)
          d=n1x(loads(data)[,-3])
          nvar=c(colnames(d),"Cluster")
          d=cbind(d,kmeans(d,clust)$cluster)
          colnames(d)=nvar
          
          
          ggpairs(d[,cbind("d13C","d18O","Mg","Fe","Mn","Sr")],aes(color="Cluster",shape="Cluster"),diag="blank",upper="blank")+
            geom_point(colour=c("Black","Gray75","Grey50")[kmeans(d,clust)$cluster])+
            labs(title= paste("Dataset",gsub("([Dataset.csv])", "", data)," initial Min-Max data", sep=""),tag="C") 
        }
      
      cpairs$plots = cpairs$plots[-(1:cpairs$nrow)]
      cpairs$yAxisLabels = cpairs$yAxisLabels[-1]
      cpairs$nrow = cpairs$nrow -1
      
      cpairs$plots = cpairs$plots[-(seq(cpairs$ncol, length(cpairs$plots), by = cpairs$ncol))]
      cpairs$xAxisLabels = cpairs$xAxisLabels[-cpairs$ncol]
      cpairs$ncol = cpairs$ncol - 1
    }
    if (inp == "D")
    {
      cpairs=
        {
          set.seed(1234)
          d=m1x(loads(data)[,-3])
          nvar=c(colnames(d),"Cluster")
          d=cbind(d,kmeans(d,clust)$cluster)
          colnames(d)=nvar
          
          
          ggpairs(d[,cbind("d13C","d18O","Mg","Fe","Mn","Sr")],diag="blank",upper="blank")+
            geom_point(colour=c("Black","Gray75","Grey50")[kmeans(d,clust)$cluster])+
            labs(title= paste("Dataset",gsub("([Dataset.csv])", "", data)," initial z-score data", sep=""),tag="D") 
        }
      
      cpairs$plots = cpairs$plots[-(1:cpairs$nrow)]
      cpairs$yAxisLabels = cpairs$yAxisLabels[-1]
      cpairs$nrow = cpairs$nrow -1
      
      cpairs$plots = cpairs$plots[-(seq(cpairs$ncol, length(cpairs$plots), by = cpairs$ncol))]
      cpairs$xAxisLabels = cpairs$xAxisLabels[-cpairs$ncol]
      cpairs$ncol = cpairs$ncol - 1
    }
    ###plot option
    cpairs
    
    ###saving option
    #ggsave(paste(data,"_Pairs",inp,".jpg",sep=""),cpairs,width = 14, height = 8,dpi=800)
  }
  
  #automatic plotting pairs
  ownpairs=function(dstart,dend)
  {
    for (j in dstart:dend)
    {
      for (i in 1:4)
      {
        if (i==1)
        {
          lower(paste("Dataset",j,".csv",sep=""),3,"A")
        }
        if (i==2)
        {
          lower(paste("Dataset",j,".csv",sep=""),3,"B")
        }
        if (i==3)
        {
          lower(paste("Dataset",j,".csv",sep=""),3,"C")
        }
        if (i==4)
        {
          lower(paste("Dataset",j,".csv",sep=""),3,"D")
        }
      }
    }
    
  }