# ===============================================
# ALGORITMA GWMR UNTUK 3 RESPONS DAN 4 PREDIKTOR
# ===============================================
library(MASS)   # Untuk fungsi ginv() (Generalized Inverse)
library(readxl) # Untuk membaca file Excel
df <- read_excel("C:/Users/LENOVO/Downloads/Dataset IKP Besar.xlsx", sheet = "Dataset")

# Mapping variabel (WAJIB dijalankan dulu sebelum scale())
df$Y1 <- df$`Indeks Ketersediaan`
df$Y2 <- df$`Indeks Keterjangkauan`
df$Y3 <- df$`Indeks Pemanfaatan`
df$X1 <- df$`Persentase Penduduk Miskin`
df$X2 <- df$`Luas Panen Padi (Ha)`
df$X3 <- df$`Persentase Sanitasi Layak`
df$X4 <- df$`Produksi Beras Per Kapita (Ton/Ribu Jiwa)`
df$X5 <- df$`Nilai Tuker Petani`
df$X6 <- df$`PDRB atas Harga Berlaku`
df$X7 <- df$`Persentase Air Minum Layak`
df$X8 <- df$`PoU`

n <- nrow(df)
m <- 3

# Standardisasi variabel X — kombinasi 
X1_s <- scale(df$X4)
X2_s <- scale(df$X6)
X3_s <- scale(df$X7)
X4_s <- scale(df$X)

# Menyusun Matriks Desain (X) dengan tambahan kolom Intersep
X <- cbind(1, X1_s, X2_s, X3_s, X4_s)
colnames(X) <- c("Intercept", "X1", "X2", "X3", "X4")

# Menyusun Matriks Respons (Y)
Y <- as.matrix(df[, c("Y1", "Y2", "Y3")])

# Membentuk Matriks Koordinat & Matriks Jarak Euclidean (u_i, v_i)
coords <- as.matrix(df[, c("Longitude", "Latitude")]) 
dist_matrix <- as.matrix(dist(coords))

# Fungsi Kernel Adaptive Bi-square
adaptive_bisquare <- function(distance_vec, h) {
  weight_vec <- ifelse(distance_vec <= h, (1 - (distance_vec / h)^2)^2, 0)
  return(diag(weight_vec))
}

# ------------------------------------------------------------------------------
# 3. OPTIMASI BANDWIDTH (MEMINIMALKAN AICc)
# ------------------------------------------------------------------------------
calculate_aicc <- function(h) {
  tr_S <- 0
  E_mat <- matrix(0, nrow = n, ncol = m)
  
  for (i in 1:n) {
    W_i <- adaptive_bisquare(dist_matrix[i, ], h)
    
    if (sum(diag(W_i)) <= 4) return(Inf) 
    
    X_Wi_X_inv <- tryCatch(ginv(t(X) %*% W_i %*% X), error = function(e) return(NULL))
    if (is.null(X_Wi_X_inv)) return(Inf)
    
    S_ii <- (X[i, , drop=FALSE] %*% X_Wi_X_inv %*% t(X[i, , drop=FALSE])) * W_i[i, i]
    tr_S <- tr_S + as.numeric(S_ii)
    
    B_hat_i <- X_Wi_X_inv %*% t(X) %*% W_i %*% Y
    Y_hat_i <- X[i, , drop=FALSE] %*% B_hat_i
    E_mat[i, ] <- Y[i, ] - Y_hat_i
  }
  
  Sigma_hat <- (t(E_mat) %*% E_mat) / n
  det_Sigma <- det(Sigma_hat)
  
  if (det_Sigma <= 0) return(Inf)
  
  AICc <- (n * log(det_Sigma)) + (n * m * log(2 * pi)) + 
    (n * (n + tr_S) * m) / (n - tr_S - m - 1)
  return(list(AICc = AICc, tr_S = tr_S))
}

cat("Mencari bandwidth optimal...\n")
h_candidates <- seq(max(dist_matrix) * 0.1, max(dist_matrix), length.out = 30)
best_h <- NULL
best_AICc <- Inf
global_eta <- NULL

for (h in h_candidates) {
  res <- calculate_aicc(h)
  if (is.list(res) && res$AICc < best_AICc) {
    best_AICc <- res$AICc
    best_h <- h
    global_eta <- res$tr_S
  }
}

cat("Bandwidth Optimal (h):", round(best_h, 4), "\n")
cat("AICc Terkecil:", round(best_AICc, 4), "\n")
cat("Effective Parameters (eta):", round(global_eta, 4), "\n\n")

# ------------------------------------------------------------------------------
# 4. FUNGSI UJI HIPOTESIS MGLH LOKAL
# ------------------------------------------------------------------------------
gwmr_mglh_local <- function(X, Y, Wi, L, M, C, eta) {
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  
  L_B_M_minus_C <- (L %*% B_hat %*% M) - C
  H_mat <- t(L_B_M_minus_C) %*% ginv(L %*% inv_X_Wi_X %*% t(L)) %*% L_B_M_minus_C
  
  E_mat <- t(M) %*% (t(Y) %*% Wi %*% Y - t(B_hat) %*% X_Wi_X %*% B_hat) %*% M
  
  V <- sum(diag(H_mat %*% ginv(H_mat + E_mat)))
  
  a <- nrow(L); b <- ncol(M); q <- min(a, b)
  theta <- (abs(a - b) - 1) / 2
  lambda <- (eta - b - 1) / 2
  
  num <- 2 * lambda + q + 1
  den <- 2 * theta + q + 1
  
  if(abs(q - V) < 1e-10) {
    p_val <- NA
  } else {
    F_star <- (num / den) * (V / (q - V))
    df1 <- a * b; df2 <- q * num
    p_val <- pf(F_star, df1, df2, lower.tail = FALSE)
  }
  
  return(list(B_hat = B_hat, p_value = p_val))
}

# ------------------------------------------------------------------------------
# 5. EKSEKUSI PENGUJIAN MGLH (PARSIAL UNTUK SELURUH KOMBINASI X DAN Y)
# ------------------------------------------------------------------------------
cat("\n=== TAHAP 5: PENGUJIAN HIPOTESIS PARSIAL MGLH ===\n")

# Menyiapkan matriks kosong untuk menyimpan 12 koefisien dan 12 p-value tiap lokasi
n_combinations <- 4 * 3
coef_matrix <- matrix(0, nrow = n, ncol = n_combinations)
pval_matrix <- matrix(0, nrow = n, ncol = n_combinations)

col_names <- c()
for(p in 1:4) {
  for(r in 1:3) {
    col_names <- c(col_names, paste0("X", p, "_Y", r))
  }
}
colnames(coef_matrix) <- paste0("Coef_", col_names)
colnames(pval_matrix) <- paste0("Pval_", col_names)

cat("Memproses pengujian untuk", n, "lokasi...\n")

for (i in 1:n) {
  Wi <- adaptive_bisquare(dist_matrix[i, ], best_h)
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  Y_Wi_Y <- t(Y) %*% Wi %*% Y
  B_X_Wi_X_B <- t(B_hat) %*% X_Wi_X %*% B_hat
  
  col_idx <- 1
  
  # 2. Iterasi untuk 4 Prediktor dan 3 Respons
  for (p in 1:4) {
    for (r in 1:3) {
      # H0: beta_pr = 0
      L_mat <- matrix(0, nrow = 1, ncol = 5)
      L_mat[1, p + 1] <- 1  # Index p+1 karena kolom 1 adalah Intersep
      
      M_mat <- matrix(0, nrow = 3, ncol = 1)
      M_mat[r, 1] <- 1
      
      C_mat <- matrix(0, nrow = 1, ncol = 1)
      
      L_B_M_minus_C <- (L_mat %*% B_hat %*% M_mat) - C_mat
      H_mat <- t(L_B_M_minus_C) %*% ginv(L_mat %*% inv_X_Wi_X %*% t(L_mat)) %*% L_B_M_minus_C
      E_mat <- t(M_mat) %*% (Y_Wi_Y - B_X_Wi_X_B) %*% M_mat
      
      V <- sum(diag(H_mat %*% ginv(H_mat + E_mat)))
      
      a <- 1; b <- 1; q <- 1; theta <- 0
      lambda <- (global_eta - b - 1) / 2
      num <- 2 * lambda + q + 1
      den <- 2 * theta + q + 1
      
      if (abs(q - V) < 1e-10) {
        p_val <- NA
      } else {
        F_star <- (num / den) * (V / (q - V))
        p_val <- pf(F_star, df1 = a * b, df2 = q * num, lower.tail = FALSE)
      }
      
      coef_matrix[i, col_idx] <- B_hat[p + 1, r]
      pval_matrix[i, col_idx] <- p_val
      
      col_idx <- col_idx + 1
    }
  }
}

# ------------------------------------------------------------------------------
# 6. OUTPUT HASIL AKHIR
# ------------------------------------------------------------------------------
output_final <- data.frame(Provinsi = df$Provinsi)
output_final <- cbind(output_final, round(coef_matrix, 4), round(pval_matrix, 5))

print(output_final)

# write.csv(output_final, "Hasil_Estimasi_MGLH_GWMR.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. PERBANDINGAN DENGAN MODEL GLOBAL (MMR)
# ------------------------------------------------------------------------------
df_model <- data.frame(Y1=Y[,1], Y2=Y[,2], Y3=Y[,3], 
                       X1=X1_s, X2=X2_s, X3=X3_s, X4=X4_s)

fit_MMR <- lm(cbind(Y1, Y2, Y3) ~ X1 + X2 + X3 + X4, data = df_model)

E_MMR <- residuals(fit_MMR)
RSS_MMR <- sum(E_MMR^2)

Sigma_MMR <- (t(E_MMR) %*% E_MMR) / n
det_Sigma_MMR <- det(Sigma_MMR)
p_params <- 5 # Jumlah parameter per respons (Intercept + 4 Prediktor)

AICc_MMR <- (n * log(det_Sigma_MMR)) + (n * m * log(2 * pi)) + 
  (n * (n + p_params) * m) / (n - p_params - m - 1)

# Menghitung RSS untuk GWMR (Berdasarkan bandwidth optimal)
E_GWMR <- matrix(0, n, m)
B_GWMR_array <- array(0, dim=c(n, p_params, m))

for(i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  B_hat_i <- ginv(t(X) %*% W_i %*% X) %*% t(X) %*% W_i %*% Y
  
  B_GWMR_array[i, , ] <- B_hat_i
  
  Y_hat_i <- X[i, , drop=FALSE] %*% B_hat_i
  E_GWMR[i, ] <- Y[i, ] - Y_hat_i
}

RSS_GWMR <- sum(E_GWMR^2)

cat("RSS MMR Global   :", round(RSS_MMR, 4), "\n")
cat("RSS GWMR Lokal   :", round(RSS_GWMR, 4), "\n")
cat("AICc MMR Global  :", round(AICc_MMR, 4), "\n")
cat("AICc GWMR Lokal  :", round(best_AICc, 4), "\n")

# -----------------------------------------------------------------
# UJI HETEROGENITAS SPASIAL (MONTE CARLO TEST UNTUK SEMUA VARIABEL)
# -----------------------------------------------------------------
var_orig <- matrix(0, nrow = 4, ncol = 3)
rownames(var_orig) <- c("X1", "X2", "X3", "X4")
colnames(var_orig) <- c("Y1", "Y2", "Y3")

for(p in 1:4) {
  for(r in 1:3) {
    var_orig[p, r] <- var(B_GWMR_array[, p + 1, r]) 
  }
}

B_sim <- 999
var_perm_array <- array(0, dim = c(B_sim, 4, 3))

set.seed(2026)
cat("Menjalankan", B_sim, "permutasi spasial untuk", 4*3, "kombinasi variabel...\n")

for(b in 1:B_sim) {
  idx_perm <- sample(1:n)
  coords_perm <- coords[idx_perm, ]
  dist_perm <- as.matrix(dist(coords_perm))
  
  coef_perm_tmp <- array(0, dim = c(n, 4, 3)) 
  
  for(i in 1:n) {
    W_i_perm <- adaptive_bisquare(dist_perm[i, ], best_h)
    XWX_inv <- tryCatch(ginv(t(X) %*% W_i_perm %*% X), error=function(e) NULL)
    
    if(!is.null(XWX_inv)) {
      B_hat_perm <- XWX_inv %*% t(X) %*% W_i_perm %*% Y
      # Ekstrak baris prediktor (2 s.d 5) untuk semua kolom respons (1 s.d 3)
      coef_perm_tmp[i, , ] <- B_hat_perm[2:5, ] 
    }
  }
  
  for(p in 1:4) {
    for(r in 1:3) {
      var_perm_array[b, p, r] <- var(coef_perm_tmp[, p, r], na.rm = TRUE)
    }
  }
  
  if(b %% 10 == 0) cat("Iterasi", b, "/", B_sim, "selesai...\n")
}

result_mc_df <- data.frame(
  Predictor = character(),
  Response = character(),
  Original_Var = numeric(),
  P_Value = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

for(p in 1:4) {
  for(r in 1:3) {
    p_val <- (1 + sum(var_perm_array[, p, r] >= var_orig[p, r])) / (1 + B_sim)
    
    status <- ifelse(p_val < 0.05, 
                     "Signifikan (Heterogen Spasial)", 
                     "Tidak Signifikan (Stasioner/Global)")
    
    result_mc_df <- rbind(result_mc_df, data.frame(
      Predictor = paste0("X", p),
      Response = paste0("Y", r),
      Original_Var = round(var_orig[p, r], 4),
      P_Value = round(p_val, 4),
      Status = status
    ))
  }
}
cat("\n=== HASIL UJI HETEROGENITAS SPASIAL ===\n")
print(result_mc_df)


cat("\n=== EVALUASI KEBAIKAN MODEL (R2, MSE, MAPE) ===\n")
local_R2 <- numeric(n)
local_MSE <- numeric(n)
local_MAPE <- numeric(n)

Y_hat_global <- matrix(0, nrow = n, ncol = m)
colnames(Y_hat_global) <- c("Y1_hat", "Y2_hat", "Y3_hat")

Y_bar <- colMeans(Y)

cat("Menghitung metrik evaluasi untuk", n, "lokasi...\n")

for (i in 1:n) {
  Wi <- adaptive_bisquare(dist_matrix[i, ], best_h)
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat_i <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  
  Y_hat_all_s <- X %*% B_hat_i
  
  Y_hat_i <- X[i, , drop = FALSE] %*% B_hat_i
  Y_hat_global[i, ] <- Y_hat_i
  
  num_R2 <- 0
  den_R2 <- 0
  for(s in 1:n) {
    if(Wi[s, s] > 0) {
      num_R2 <- num_R2 + Wi[s, s] * sum((Y[s, ] - Y_hat_all_s[s, ])^2)
      den_R2 <- den_R2 + Wi[s, s] * sum((Y[s, ] - Y_bar)^2)
    }
  }
  local_R2[i] <- 1 - (num_R2 / den_R2)
  
  local_MSE[i] <- mean((Y[i, ] - Y_hat_i)^2)
  
  local_MAPE[i] <- mean(abs((Y[i, ] - Y_hat_i) / Y[i, ])) * 100
}

global_MSE <- mean((Y - Y_hat_global)^2)
global_MAPE <- mean(abs((Y - Y_hat_global) / Y)) * 100

SST_global <- sum(sweep(Y, 2, Y_bar)^2)
SSE_global <- sum((Y - Y_hat_global)^2)
global_R2 <- 1 - (SSE_global / SST_global)

cat("\n--- RINGKASAN METRIK GLOBAL ---\n")
cat("Global R-Squared :", round(global_R2, 4), "\n")
cat("Global MSE       :", round(global_MSE, 4), "\n")
cat("Global MAPE      :", round(global_MAPE, 2), "%\n")

output_final$Local_R2 <- round(local_R2, 4)
output_final$Local_MSE <- round(local_MSE, 4)
output_final$Local_MAPE <- round(local_MAPE, 2)

print(output_final[, c("Provinsi", "Local_R2", "Local_MSE", "Local_MAPE")])