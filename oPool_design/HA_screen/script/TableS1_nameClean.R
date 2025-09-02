
# Update the antibodies' name
TbS1 <- read.csv('data/TableS1.csv',  header = T, skip = 1)
TbS1$OpName <- gsub("_", "-", TbS1$Name)
TbS1$OpName <- gsub("-Heavy", "", TbS1$OpName)

write.csv(TbS1, 'data/TableS1_nameClean.csv', row.names = F)

