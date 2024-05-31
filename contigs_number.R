library(openxlsx)
library(tidyverse)

data <- read.xlsx("total-600.final.confidence_table.xlsx")
colnames(data)[1] <- "X1"

data_split <- data %>%
  separate(X1, into = c("SRA", "other_cols"), sep = "_", extra = "merge")

data_split_selected <- select(data_split, 1, 2, 5, 6, 7, 8, 9,10,17, 18,19)

data_unique <- data_split_selected %>% 
  distinct(cluster,  SRA, other_cols, NR_blastx) %>% 
  count(SRA, cluster, name = "unique_cols") %>% 
  mutate(SRA_cluster = paste(SRA, cluster, sep = "_")) %>% 
  select(-c(SRA, cluster))

data_split_star <- filter(data_split, grepl("\\*", cor))%>% select(1, 9, 17, 18,19) %>%
                   mutate(SRA_cluster = paste(SRA, cluster, sep = "_")) %>% 
                   select(-c(SRA, cluster,cor)) %>% 
                   extract(NR_blastx, into = "virus_name", "\\[(.*)\\]", remove = FALSE)%>% 
                   select(-c(NR_blastx))

merged_data <- left_join(data_unique, data_split_star, by = "SRA_cluster")

merged_data_split <- merged_data %>%
  separate(SRA_cluster, into = c("SRA", "other_cols"), sep = "_", extra = "merge") %>%  select(-c(other_cols))

wide_data <- merged_data_split %>%
  pivot_wider(names_from = virus_name, values_from = unique_cols)

write.xlsx(wide_data, "total-600.final.wide_data.xlsx")

