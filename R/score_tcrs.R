#' score TCR sequences
#'
#' @param data TCR sequence data
#' @param chain whether the TCR data is paired a-b ("ab"), alpha only ("a") or beta only ("b")
#' @export
score_tcrs <- function(data, chain, MAIT_NKT = FALSE){
  score_names = c("TCR-innate", "TCR-CD8", "TCR-reg", "TCR-mem")
  ftz = featurize_tcrs(data, chain)
  print("TCRs featurized!")
  if (chain=="ab"){
    weights = weightsAB
    score_means = ABscore_mns
    score_sds = ABscore_sds
  } else if (chain=="a"){
    weights = weightsA
    score_names = gsub("TCR", "TCRalpha", score_names)
    score_means = Ascore_mns
    score_sds = Ascore_sds
  } else if (chain=="b"){
    weights = weightsB
    score_names = gsub("TCR", "TCRbeta", score_names)
    if (MAIT_NKT == TRUE){
      weights = weightsMAITNKT
      score_names = c("TCRbeta-MAIT", "TCRbeta-NKT")
      score_means = MNKTscore_mns
      score_sds = MNKTscore_sds
    } else {
      score_means = Bscore_mns
      score_sds = Bscore_sds
    }
  } else {
    print("please specify the 'chain' argument (a, b, or ab)")
  }
  m = mns[as.character(rownames(weights)),]
  s = sds[as.character(rownames(weights)),]
  print("scoring TCRs...")
  rownames(ftz) = as.character(ftz$id)
  ftz = scale_variables(ftz, m, s)
  scores = as.matrix(ftz) %*% as.matrix(weights)
  rownames(scores) = rownames(ftz)
  colnames(scores) = score_names
  names(score_mns) = score_names
  names(score_sds) = score_names
  scores = scale_variables(scores, score_mns, score_sds)
  ##scaling factors
  print("all done!")
  return(data.frame(scores))
}
