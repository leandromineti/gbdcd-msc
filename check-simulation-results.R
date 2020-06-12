rm(list=ls())

library(dplyr)

df_cluster_1_c_0001 <- as.data.frame(get(load(file = './results/result_cluster_one_t_100_c_1.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 1) & (1 <= hpd_2))) %>% 
  mutate(cluster_real_number = 1,
         priori_param = 0.001)
df_cluster_1_c_0333 <- as.data.frame(get(load(file = './results/result_cluster_one_t_100_c_2.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 1) & (1 <= hpd_2))) %>% 
  mutate(cluster_real_number = 1,
         priori_param = 0.333)
df_cluster_2_c_0001 <- as.data.frame(get(load(file = './results/result_cluster_diagonal_t_100_c_1.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 2) & (2 <= hpd_2))) %>% 
  mutate(cluster_real_number = 2,
         priori_param = 0.001)
df_cluster_2_c_0333 <- as.data.frame(get(load(file = './results/result_cluster_diagonal_t_100_c_2.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 2) & (2 <= hpd_2))) %>% 
  mutate(cluster_real_number = 2,
         priori_param = 0.333)
df_cluster_4_c_0001 <- as.data.frame(get(load(file = './results/result_cluster_four_t_100_c_1.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 4) & (4 <= hpd_2))) %>% 
  mutate(cluster_real_number = 4,
         priori_param = 0.001)
df_cluster_4_c_0333 <- as.data.frame(get(load(file = './results/result_cluster_four_t_100_c_2.RData'))) %>% 
  mutate(hpd_includes_true_k = ((hpd_1 <= 4) & (4 <= hpd_2))) %>% 
  mutate(cluster_real_number = 4,
         priori_param = 0.333)

df_results <- df_cluster_1_c_0001 %>% 
  bind_rows(df_cluster_1_c_0333) %>% 
  bind_rows(df_cluster_2_c_0001) %>% 
  bind_rows(df_cluster_2_c_0333) %>% 
  bind_rows(df_cluster_4_c_0001) %>% 
  bind_rows(df_cluster_4_c_0333)

df_final <- df_results %>% 
  group_by(cluster_real_number, priori_param) %>% 
  summarise(mean_k_hat = mean(mode_k), 
            mean_hpd_size = mean(hpd_size), 
            hpd_includes_true_k = mean(hpd_includes_true_k)) %>% 
  ungroup()
