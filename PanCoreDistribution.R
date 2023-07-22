###########
###author: Chengyu Gao

if (!require(iterpc))
  install.packages("iterpc",repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require(tidyverse))
  install.packages("tidyverse",repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

library(iterpc)
library(tidyverse)


args<-commandArgs(TRUE)
inputfile<-args[1]
step<-as.integer(args[3])
startNum<-as.integer(args[2])
frequen<-as.integer(args[4])
outdir<-args[5]

setwd(outdir)

df<-read.csv(inputfile,row.names = 1,header = T)


sampleNum<-ncol(df)
for (i in 1:sampleNum){
  df[which(df[,i]!=0),i]=1
}

Sample<-c()
Core<-c()
Pan<-c()



while (startNum<=sampleNum){
  combNum<-choose(sampleNum,startNum)
  combination<-iterpc(sampleNum,startNum)
  
  if (combNum <frequen){
    for (i in 1:sampleNum){
      result<-getnext(combination)
      aRowSum<-apply(df[result], 1, sum)
      cores<-length(grep(startNum,aRowSum,fixed = T))
      pans<-length(aRowSum)-length(grep(0,aRowSum,fixed = T))
      Core<-append(Core,cores)
      Pan<-append(Pan,pans)
      Sample<-append(Sample,startNum)
    }
  }else{
    randomSample<-sample(sampleNum,frequen)
    for (i in 1:sampleNum){
      if (i %in% randomSample){
        result<-getnext(combination)
        aRowSum<-apply(df[result], 1, sum)
        cores<-length(grep(startNum,aRowSum,fixed = T))
        pans<-length(aRowSum)-length(grep(0,aRowSum,fixed = T))
        Core<-append(Core,cores)
        Pan<-append(Pan,pans)
        Sample<-append(Sample,startNum)
      }
    }
  }
  startNum=startNum+step
}

outdf<-as.data.frame(cbind(Sample,Core,Pan))
write.table(outdf,"results.txt",sep = "\t",quote = F,row.names = F)
#绘图

longdf<-pivot_longer(outdf,cols = c("Core","Pan"),names_to = "group",values_to = "value")
pdf("plot.pdf",width = 8,height = 5)
ggplot(longdf,aes(x=Sample,y=value,color=group))+
  geom_point(shape=4)+
  geom_smooth(aes(x=Sample,y=value),
              method = "lm")+
  geom_smooth(aes(x=Sample,y=value),
              method = "lm")+
  scale_color_manual(values=c("#D17431","#39729C"))+
  labs(x="Number of genomes",y="Number of gene families",color="")+
  theme_classic()
dev.off()

coredf<-longdf%>%filter(group=="Core")
pandf<-longdf%>%filter(group=="Pan")
corem<-summary(lm(value~Sample,data = coredf))
panm<-summary(lm(value~Sample,data = pandf))

corem_df<-data.frame(Slope=corem$coefficients[2,1],Intercept=corem$coefficients[1,1],R=corem$adj.r.squared)

panm_df<-data.frame(Slope=panm$coefficients[2,1],Intercept=panm$coefficients[1,1],R=panm$adj.r.squared)


write.table(corem_df,"core_formula.txt",sep = "\t",quote = F,row.names = F)
write.table(panm_df,"pan_formula.txt",sep = "\t",quote = F,row.names = F)