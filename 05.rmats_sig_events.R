### DAS events
# |ΔPSI| > 0.1 & FDR < 0.05
# ============================================
setwd("/workpath/AS_out/")
library(dplyr)
library(tidyr)
library(stringr)

# file list
event_files <- c(
  "A3SS.MATS.JCEC.txt",
  "A5SS.MATS.JCEC.txt", 
  "MXE.MATS.JCEC.txt",
  "RI.MATS.JCEC.txt",
  "SE.MATS.JCEC.txt"
)

event_types <- c("A3SS", "A5SS", "MXE", "RI", "SE")

all_results <- list()

for (i in seq_along(event_files)) {
  file_path <- event_files[i]
  event_type <- event_types[i]
  
  # read files
  df <- read.table(file_path, header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE)
  # filter: |ΔPSI| > 0.1 & FDR < 0.05
  df_filtered <- df %>%
    filter(abs(IncLevelDifference) > 0.1, FDR < 0.05)
  result <- df_filtered %>%
        select(
          event_id = ID,
          gene_id = GeneID,
          gene_symbol = geneSymbol,
          psi_diff = IncLevelDifference,
          fdr = FDR
        ) %>%
        mutate(event_type = event_type, .before = 1)
  
    all_results[[event_type]] <- result
}

# bind rows
final_result <- bind_rows(all_results)

# print result

summary_stats <- final_result %>%
  group_by(event_type) %>%
  summarise(count = n(), .groups = "drop")
print(summary_stats)

cat("\nall DAS:", nrow(final_result), "\n")

# save result
write.table(final_result, "significant_splicing_events.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write.csv(final_result, "significant_splicing_events.csv", 
          row.names = FALSE, quote = FALSE)

write.csv(summary_stats, "sig_event_counts.csv", 
          row.names = FALSE, quote = FALSE)
