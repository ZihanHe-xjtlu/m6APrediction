#' Encode DNA 5-mer sequences into factor-based features
#'
#' @description
#' This function converts a vector of DNA 5-mer strings into a feature
#' data frame where each nucleotide position is treated as a categorical
#' variable. The encoded features are used as input for the m6A prediction model.
#'
#' @param dna_strings A character vector of DNA strings (e.g., 5-mer sequences).
#'
#' @return A data frame where each column represents a nucleotide position
#' (nt_pos1–nt_pos5), and each value is a factor with levels A/T/C/G.
#'
#' @import randomForest
#' @export
#'
#' @examples
#' dna_encoding(c("ATCGA", "TGGCA"))
#'
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A methylation for multiple input sites
#'
#' @description
#' This function generates m6A prediction probability and predicted class
#' ("Positive"/"Negative") for a batch of input RNA sites. It automatically
#' encodes the DNA_5mer sequence feature using \code{dna_encoding()},
#' combines the encoded features with the original feature table,
#' and runs the random forest model.
#'
#' @param ml_fit A fitted randomForest model object (loaded from rf_fit.rds).
#' @param feature_df A data frame containing input features:
#' \code{gc_content}, \code{RNA_type}, \code{RNA_region},
#' \code{exon_length}, \code{distance_to_junction},
#' \code{evolutionary_conservation}, and \code{DNA_5mer}.
#' @param positive_threshold A numeric cutoff for classifying output
#' probability into Positive/Negative. Default = 0.5.
#'
#' @return A data frame containing the original input columns plus:
#' \itemize{
#'   \item \code{predicted_m6A_prob} — predicted probability for m6A modification
#'   \item \code{predicted_m6A_status} — predicted class ("Positive"/"Negative")
#' }
#'
#'@importFrom stats predict
#' @export
#'
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                            package = "m6APrediction"))
#' pred_result <- prediction_multiple(model, df, positive_threshold = 0.5)
#' head(pred_result[, c("gc_content", "RNA_type", "predicted_m6A_prob", "predicted_m6A_status")])
#'

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_encoded <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_encoded$RNA_type <- factor(feature_encoded$RNA_type,
                                     levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_encoded$RNA_region <- factor(feature_encoded$RNA_region,
                                       levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  pred_prob <- predict(ml_fit, newdata = feature_encoded, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- ifelse(pred_prob > positive_threshold, "Positive", "Negative")

  return(feature_df)
}

#' Predict m6A methylation for a single input site
#'
#' @description
#' This function performs m6A prediction for a single RNA site. It internally
#' constructs a one-row data frame from user-provided values and calls
#' \code{prediction_multiple()} to generate prediction probability and class.
#'
#' @param ml_fit A fitted randomForest model object (loaded from rf_fit.rds).
#' @param gc_content Numeric value of GC content.
#' @param RNA_type RNA type as a string; must be one of:
#'   "mRNA", "lincRNA", "lncRNA", "pseudogene".
#' @param RNA_region RNA region as a string; must be one of:
#'   "CDS", "intron", "3'UTR", "5'UTR".
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to nearest exon–intron junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer A 5-mer DNA sequence string.
#' @param positive_threshold Numeric probability threshold (default = 0.5).
#'
#' @return A named vector containing:
#' \itemize{
#'   \item \code{predicted_m6A_prob}
#'   \item \code{predicted_m6A_status}
#' }
#'
#' @export
#'
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' pred_single_result <- prediction_single(
#'     ml_fit = model,
#'     gc_content = 0.55,
#'     RNA_type = "mRNA",
#'     RNA_region = "CDS",
#'     exon_length = 1500,
#'     distance_to_junction = 120,
#'     evolutionary_conservation = 0.32,
#'     DNA_5mer = "ATCGA"
#' )
#' print(pred_single_result)

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  single_df$RNA_type <- factor(single_df$RNA_type,
                               levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  single_df$RNA_region <- factor(single_df$RNA_region,
                                 levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  pred_result <- prediction_multiple(ml_fit, feature_df = single_df, positive_threshold = positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = pred_result$predicted_m6A_prob,
    predicted_m6A_status = pred_result$predicted_m6A_status
  )
  return(returned_vector)
}
