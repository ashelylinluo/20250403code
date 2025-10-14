=== 路径 ===
DST <- "/home/data/t210549/luolin/WGCNA_WRKY_Figure6"
setwd(DST)

=== 加载现有 Rdata（稳健处理不同文件名） ===
step1: datExpr, traitData/datTraits
load(file.path(DST, "power_8/step1_input.Rdata"))

step3: moduleColors（注意你保存的是 *_50.Rdata）
if (file.exists(file.path(DST, "power_8/step3_genes_modules_50.Rdata"))) {
  load(file.path(DST, "power_8/step3_genes_modules_50.Rdata"))  # 提供 moduleColors/net
} else {
  load(file.path(DST, "power_8/step3_genes_modules.Rdata"))
}

step4: design（师弟保存为 *_50.Rdata）
if (file.exists(file.path(DST, "power_8/step4_design_50.Rdata"))) {
  load(file.path(DST, "power_8/step4_design_50.Rdata"))         # 提供 design
} else if (file.exists(file.path(DST, "power_8/step4_design.Rdata"))) {
  load(file.path(DST, "power_8/step4_design.Rdata"))
} else {
  stop("找不到 step4_design*.Rdata，请检查 power_8/ 下的文件。")
}

=== 对齐样本（极重要） ===
if (exists("datTraits")) traitData <- datTraits
stopifnot(all(rownames(datExpr) %in% rownames(traitData)))
traitData <- traitData[rownames(datExpr), , drop = FALSE]

=== 不用 WGCNA：手动计算模块特征向量 (ME = 每个模块基因矩阵的PC1) ===
compute_MEs <- function(datExpr, moduleColors) {
  mods <- unique(moduleColors); mods <- mods[mods != "grey"]
  ME_df <- sapply(mods, function(mod){
    idx <- moduleColors == mod
    X   <- as.matrix(datExpr[, idx, drop = FALSE])     # 样本×基因
    Xs  <- scale(X, center = TRUE, scale = TRUE)
    pc  <- prcomp(Xs, center = FALSE, scale. = FALSE)
    me  <- pc$x[,1]
    sgn <- sign(cor(me, rowMeans(Xs), use = "pairwise.complete.obs"))
    if (is.na(sgn) || sgn == 0) sgn <- 1
    me * sgn
  })
  ME_df <- as.data.frame(ME_df)
  colnames(ME_df) <- paste0("ME", mods)
  rownames(ME_df) <- rownames(datExpr)
模块间相关性聚类后重排（近似 orderMEs）
  cmat <- cor(ME_df, use = "pairwise.complete.obs")
  ME_df[, hclust(as.dist(1 - cmat))$order, drop = FALSE]
}

MEs <- compute_MEs(datExpr, moduleColors)

=== 用 design（CK/drought/recovery 的0/1矩阵）做相关性与显著性 ===
如果 design 列名不是 control/drought/recovery 请按热图对应修改：
X <- design   # 列名应该是 CK / drought_14d / rewater_24h（或 control/drought/recovery）

相关性矩阵与 p 值
cor_mat <- cor(MEs, X, use = "p")    # 行=ME模块，列=性状
n <- nrow(MEs)
p_mat  <- matrix(NA, nrow = ncol(MEs), ncol = ncol(X),
                 dimnames = list(colnames(MEs), colnames(X)))
for (j in seq_len(ncol(X))) {
  p_mat[, j] <- corPvalueStudent(cor_mat[, j], n)      # 不依赖 WGCNA 的版本
}

FDR（按列或全局都可以；这里给“全局更保守”的一份）
p_vec <- as.vector(p_mat)
fdr_vec <- p.adjust(p_vec, method = "BH")
FDR_global <- matrix(fdr_vec, nrow = nrow(p_mat), ncol = ncol(p_mat),
                     dimnames = dimnames(p_mat), byrow = FALSE)

dir.create(file.path(DST, "power_8/exports"), showWarnings = FALSE, recursive = TRUE)
write.csv(cor_mat,      file.path(DST, "power_8/exports", "ModuleTraitCor.csv"))
write.csv(p_mat,        file.path(DST, "power_8/exports", "ModuleTraitPvalue.csv"))
write.csv(FDR_global,   file.path(DST, "power_8/exports", "ModuleTraitFDR_global.csv"))

=== 按阈值筛选“与胁迫显著相关的模块” ===
根据你的热图标签，选择胁迫那一列（drought 或 drought_14d）
stress_col <- intersect(colnames(X), c("drought", "drought_14d"))[1]
recovery_col <- intersect(colnames(X), c("recovery", "rewater_24h"))[1]
control_col  <- intersect(colnames(X), c("control", "CK"))[1]

stopifnot(!is.na(stress_col))

effect_cut <- 0.25   # 可调到 0.30 更严格
res_stress <- data.frame(
  Module    = rownames(cor_mat),
  Cor       = cor_mat[, stress_col],
  Pvalue    = p_mat[, stress_col],
  FDR_global= FDR_global[, stress_col],
  row.names = NULL
) |>
  dplyr::arrange(FDR_global)

sig_stress <- res_stress |>
  dplyr::filter(Pvalue < 0.05 & abs(Cor) >= effect_cut)

write.csv(res_stress, file.path(DST, "power_8/exports", "Stress_cor_full.csv"), row.names = FALSE)
write.csv(sig_stress, file.path(DST, "power_8/exports", "Stress_cor_significant.csv"), row.names = FALSE)

需要恢复期的也来一份（如有）
if (!is.na(recovery_col)) {
  res_recovery <- data.frame(
    Module    = rownames(cor_mat),
    Cor       = cor_mat[, recovery_col],
    Pvalue    = p_mat[, recovery_col],
    FDR_global= FDR_global[, recovery_col],
    row.names = NULL
  ) |>
    dplyr::arrange(FDR_global)
  sig_recovery <- res_recovery |>
    dplyr::filter(Pvalue < 0.05 & abs(Cor) >= effect_cut)
  write.csv(res_recovery, file.path(DST, "power_8/exports", "Recovery_cor_full.csv"), row.names = FALSE)
  write.csv(sig_recovery, file.path(DST, "power_8/exports", "Recovery_cor_significant.csv"), row.names = FALSE)
}

cat("完成：请在 power_8/exports/ 查看 ModuleTraitCor.csv / ModuleTraitPvalue.csv / Stress_cor_significant.csv 等文件。\n")
