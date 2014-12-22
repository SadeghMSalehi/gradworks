suppressMessages(require(gtools))
suppressMessages(require(ggplot2))

combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1
     || n < 1 || (n %% 1) != 0) stop("bad value of n")
  if(mode(r) != "numeric" || length(r) != 1
     || r < 1 || (r %% 1) != 0) stop("bad value of r")
  if(!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
    {
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(n == 1) matrix(v, 1, r) else
            rbind( cbind(v[1], Recall(n, r-1, v)),
                   Recall(n-1, r, v[-1]))
    }
  else
    sub <- function(n, r, v)
    {
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(r == n) matrix(v, 1, n) else
            rbind(cbind(v[1], Recall(n-1, r-1, v[-1])),
                  Recall(n-1, r, v[-1]))
    }
  sub(n, r, v[1:n])
}

permhist = function(o2, pvalues, observed, ppvalue, count, total) {
  pdf(o2,width=8,height=6)
  hist(pvalues, main=sprintf("Permutation Histogram (p=%5.4f, n=%d, N=%d)", ppvalue, count, total), xlab = "mean", xlim=c(0,.02))
  par(new=TRUE)   
  abline(v=observed, col="blue", lty=2)
  par(new=TRUE)   
  plot(density(pvalues), col=2, yaxt="n", xaxt="n", bty='n', xlab="", ylab="", main='', xlim=c(0,.02))
  dev.off()
}

permtest = function(f1,f2,o1,o2=NULL) {
  data1.input_df = read.table(f1)[,-1]
  data2.input_df = read.table(f2)[,-1]
  
  data1 = as.matrix(data1.input_df)
  data2 = as.matrix(data2.input_df)
  
  combined = cbind(data1,data2)
  diff.observed = mean(c(data1)) - mean(c(data2))
  diff.abs_observed = abs(mean(c(data1)) - mean(c(data2)))
  diff.pvalue_observed = t.test(c(data1), c(data2))$p.value
  
  diff.random = NULL
  diff.ttest = NULL
  diff.pvalues = NULL
  
  len1 = ncol(data1)
  len2 = ncol(combined)
  
  indexes = t(combinations(len2,len1))
  indexes = indexes[,1:(ncol(indexes)/2)]
  #indexes.all = indexes.all[,(1:ncol(indexes.all)/2)]
  
  for (i in 1:ncol(indexes)) {
    data1.random = combined[,indexes[,i]]
    data2.random = combined[,-indexes[,i]]
    
    diff.random[i] = mean(c(data2.random)) - mean(c(data1.random))
    diff.pvalues[i] = t.test(c(data2.random), c(data1.random))$p.value
  }
  write.table(diff.pvalues, o1)
  
  diff.abs_random = abs(diff.random)
  diff.mean_max = max(diff.abs_random)
  
  count = length(diff.pvalues[diff.pvalues <= diff.pvalue_observed])
  total = ncol(indexes)
  pvalue = count / total
  
  permhist(o2, diff.pvalues, diff.pvalue_observed, pvalue, count, total)
  return (diff)
}


# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) >= 4) {
#   dir = args[1]
#   data1_prefix = args[2]
#   data2_prefix = args[3]
#   output_prefix = args[4]
#   
#   if (!file.exists(dir)) {
#     stop(sprintf("'%s' not exists", dir))
#   }
# } else {
#   stop(sprintf("permtest.R root-dir data1-prefix data2-prefix output-prefix"))
# }
#dir = "/Prime/Rodent-Thickness-Data/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right"

# for (j in 1:20) {
#   region_number = j
#   data1_file = sprintf("%s/Statistical/%s%02d.txt", dir, data1_prefix, region_number)
#   data2_file = sprintf("%s/Statistical/%s%02d.txt", dir, data2_prefix, region_number)
#   data_out = sprintf("%s/Statistical/%s%02d.txt", dir, output_prefix, region_number)
#   pdf_out = sprintf("%s/Statistical/%s%02d.pdf", dir, output_prefix, region_number)
#   if (!file.exists(data1_file)) {
#     cat(sprintf("'%s' not exists", data1_file))
#     break
#   } 
#   if (!file.exists(data2_file)) {
#     cat(sprintf("'%s' not exists", data2_file))
#     break
#   }
#   diff = permtest(data1_file, data2_file, data_out, pdf_out)
# }