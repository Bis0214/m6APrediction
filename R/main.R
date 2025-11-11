#' Encode DNA Sequence
#' This function encodes a DNA sequence into a suitable format for machine learning.
#' @param dna_strings A character string representing the DNA sequence to be encoded.
#' @return A numeric vector representing the encoded sequence.
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict Multiple m6A Sites
#' This function predicts multiple m6A sites based on input features using a trained model.
#' @param ml_fit A trained model (e.g., a random forest model) used for making predictions.
#' @param feature_df A data frame containing the features for prediction. The data frame must include:
#'   - `gc_content`: The GC content of the RNA sequence (numeric).
#'   - `RNA_type`: The type of RNA (character, levels: "mRNA", "lincRNA", "lncRNA", "pseudogene").
#'   - `RNA_region`: The region of the RNA sequence (character, levels: "CDS", "intron", "3'UTR", "5'UTR").
#'   - `exon_length`: The length of the exon (numeric).
#'   - `distance_to_junction`: The distance from the m6A site to the nearest splice junction (numeric).
#'   - `evolutionary_conservation`: A measure of evolutionary conservation at the m6A site (numeric).
#'   - `DNA_5mer`: The 5-mer surrounding the m6A site (character).
#' @param positive_threshold A threshold for positive predictions (numeric, default is 0.5).
#' @return A vector of predictions based on the input data.
#' @import randomForest
#' @export
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' input_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(ml_fit = model, feature_df = input_data)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content","RNA_type","RNA_region","exon_length","distance_to_junction","evolutionary_conservation","DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA","lincRNA","lncRNA","pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS","intron","3'UTR","5'UTR"))

  seq_df <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df, seq_df)

  probs <- predict(ml_fit, newdata = feature_df, type = "prob")[,"Positive"]
  feature_df$predicted_m6A_prob <- probs
  feature_df$predicted_m6A_status <- ifelse(probs > positive_threshold, "Positive", "Negative")

  return(feature_df)
}

#' Predict a Single m6A Site
#' This function predicts a single m6A site based on input features using a trained model.
#' @param ml_fit A trained model (e.g., a random forest model) used for making predictions.
#' @param gc_content The GC content of the RNA sequence as a numeric value.
#' @param RNA_type The type of RNA (e.g., mRNA, lncRNA) as a character string.
#' @param RNA_region The region of the RNA sequence (e.g., coding region, UTR) as a character string.
#' @param exon_length The length of the exon containing the m6A site as a numeric value.
#' @param distance_to_junction The distance from the m6A site to the nearest splice junction as a numeric value.
#' @param evolutionary_conservation A measure of the evolutionary conservation at the m6A site as a numeric value.
#' @param DNA_5mer The 5-mer DNA sequence surrounding the m6A site as a character string.
#' @param positive_threshold The threshold for positive predictions, default is 0.5.
#' @return A single prediction value based on the input vector.
#' @importFrom stats predict
#' @export
#' @examples
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' input_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' gc_content <- input_data$gc_content
#' RNA_type <- input_data$RNA_type
#' RNA_region <- input_data$RNA_region
#' exon_length <- input_data$exon_length
#' distance_to_junction <- input_data$distance_to_junction
#' evolutionary_conservation <- input_data$evolutionary_conservation
#' DNA_5mer <- input_data$DNA_5mer
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = gc_content,
#'   RNA_type = RNA_type,
#'   RNA_region = RNA_region,
#'   exon_length = exon_length,
#'   distance_to_junction = distance_to_junction,
#'   evolutionary_conservation = evolutionary_conservation,
#'   DNA_5mer = DNA_5mer,
#'   positive_threshold = 0.5
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = c("mRNA","lincRNA","lncRNA","pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS","intron","3'UTR","5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  seq_df <- dna_encoding(single_df$DNA_5mer)
  single_df <- cbind(single_df, seq_df)

  probs <- predict(ml_fit, single_df, type = "prob")[,"Positive"]
  c(predicted_m6A_prob = probs, predicted_m6A_status = ifelse(probs > positive_threshold, "Positive", "Negative"))
}
