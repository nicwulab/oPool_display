library(ggplot2)

TB = read.csv("result/20240305_Fragmented.csv")

LEN = c()
for(i in c(1:nrow(TB))){
  LEN = c(LEN, length(strsplit(TB$Seg1[i], "")[[1]]))
}
TB["len_Seg1"] = LEN

ggplot(data = TB, aes(x= len_Seg1)) + geom_density() + theme_bw()


TB$VL_AA[TB$len_Seg1>=316]
