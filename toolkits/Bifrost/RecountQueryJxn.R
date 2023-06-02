library("recount3")
library(dplyr)
#library("recount")
#library(biomaRt)

ref <- read.csv("/home/rmaruzani/storage/hydra2/vus/toolkits/bifrost/reference.csv")
alt <- read.csv("/home/rmaruzani/storage/hydra2/vus/toolkits/bifrost/alternate.csv")

human_projects <- available_projects()

process_data <- function(df, human_projects) {
  # Group the dataframe by Variant.ID
  grouped_df <- df %>% group_by(Variant.ID)
  
  # Loop through each variant
  for (variant in unique(df$Variant.ID)) {
    # Create a file name based on the dataframe and variant
    file_name <- paste0(deparse(substitute(df)),"_", variant, ".output.txt")
    
    # Open the file in write mode
    file_conn <- file(file_name, "w")
    
    # Filter rows for the current variant
    variant_df <- grouped_df %>% filter(Variant.ID == variant)
    # Get unique studies for the variant
    unique_studies <- unique(variant_df$Study)
    
    # Loop through the unique studies (get all junctions from recount3)
    for (study in unique_studies) {
      proj_info <- subset(human_projects, project == study & project_type == "data_sources")
      dwnload <- create_rse(proj_info, type = "jxn")
      rse_jxn <- create_rse(proj_info, type = "jxn")
      table_data <- as_tibble(rowRanges(rse_jxn))
      
      # Filter rows for the current study
      study_df <- variant_df %>% filter(Study == study)
      # Loop through the filtered rows- should only be one - and 
      #return the junctions that fall within the gene coordinates 
      for (i in 1:nrow(study_df)) {
        curr_row <- study_df[i, ]
        chr_value <- curr_row$Chr
        
        if (!is.na(chr_value)) {
          for (i in 1:nrow(table_data)) {
            if (!is.na(table_data$seqnames[i]) && table_data$seqnames[i] == chr_value &&
                !is.na(table_data$start[i]) && !is.na(table_data$end[i]) &&
                table_data$start[i] >= curr_row$`Gene.Start` && table_data$end[i] <= curr_row$`Gene.Stop`) {
              output_string <- paste(table_data$seqnames[i], table_data$start[i], table_data$end[i], sep = "-") # nolint: line_length_linter.
              
              # Write the output string (junctions) to the file
              writeLines(output_string, file_conn)
              
              # Print some details to the console
              print(paste("variant:", variant, "study:", study, "gene:", curr_row$Gene, "jxn:", output_string))
            }
          }
        } else {
          print("Missing value in Chr vector.")
        }
      }
    }
    
    # Close the file connection
    close(file_conn)
  }
}






process_data(ref, human_projects)
process_data(alt, human_projects)