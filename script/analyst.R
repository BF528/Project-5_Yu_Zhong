## set wd
setwd("/projectnb/bf528/users/dreadlocks/project_5_yu_zhong")

## load datasets
cuffdiff_table <- read.table("result/cuffdiff/gene_exp.diff", header = T)
# sort datasets
sorted_cuffdiff <- cuffdiff_table[order(cuffdiff_table$q_value, decreasing = F),]
top_ten <- sorted_cuffdiff[, c("gene","value_1","value_2",
                               "log2.fold_change.","p_value","q_value")] %>% head(10)
top_ten

## examine the expression level
par(mfrow = c(1,2))
hist(cuffdiff_table$log2.fold_change., main = "",
     xlab = "log2FC", breaks = 20, xlim = c(-10,10), probability = T, 
     col = "red")
title(main = "A", adj = 0)

# remove not significant ones
sign_cuffdiff <- subset(cuffdiff_table, significant == "yes")

hist(sign_cuffdiff$log2.fold_change., main = "",
     xlab = "log2FC", breaks = 20, xlim = c(-10,10), probability = T,
     col = "blue")
title(main = "B", adj = 0)

## select up and down degs
up_gene <- subset(sign_cuffdiff, log2.fold_change. > 0 & q_value <= 0.01)
down_gene <- subset(sign_cuffdiff, log2.fold_change. < 0 & q_value <= 0.01)

## export the results
write.csv(up_gene$gene, file = "result/up_gene.csv", 
          row.names = F, quote = F)
write.csv(down_gene$gene, file = "result/down_gene.csv",
          row.names = F, quote = F)
