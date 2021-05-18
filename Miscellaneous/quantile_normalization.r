## Quantile normalization

```r
#Load data
df <- read.table("toquantile.txt",header=T)

head(data)

rownames(df)<-df$ID

df_rank <- apply(df,2,rank,ties.method="min")
df_sorted <- data.frame(apply(df, 2, sort))

# Create quantile normalization function
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

new_data <- quantile_normalisation(df[,2:37]) #Set data columns
new_data

write.table(new_data,"quantiled_data.txt",sep="\t",na="",quote=F)
```

```r
# Alternatively, you can use preprocessCore package on Bioconductor

# Instalation
source('http://bioconductor.org/biocLite.R')
biocLite('preprocessCore')
# Load package
library(preprocessCore)
# Function expects a matrix
normalize.quantiles(mat)
```