for (k in 1:462) {
  f1 = sprintf("/Prime/Thesis-Data/Rodent-Thickness/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right/Statistical/Permtest/dataA_%03d.txt", k)
  
  f2 = sprintf("/Prime/Thesis-Data/Rodent-Thickness/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right/Statistical/Permtest/dataB_%03d.txt", k)
  
  d1 = read.table(f1);
  d2 = read.table(f2);
  
  pvalues = NULL
  for (j in 1:nrow(d1)) {
    pvalues[j] = t.test(d1[j,],d2[j,])$p.value
  }
  
  write.table(pvalues, sprintf("/Prime/Thesis-Data/Rodent-Thickness/RPV_Thickness/RPV1+3/RPV3_AIE_vs_RPV3_Control_Right/Statistical/Permtest/pvalues_%03d.txt", k),row.names=FALSE,col.names=FALSE)
}