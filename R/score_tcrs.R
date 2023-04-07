#' @export
score_tcrs <- function(data, chain){
  ftz = featurize_tcrs(data, chain)
  if (chain=="ab"){
    ldgs = ABldgs
    mns_x = ABmns_x
    sds_x = ABsds_x
  } else if (chain=="a"){
    ldgs = Aldgs
    mns_x = Amns_x
    sds_x = Asds_x
  } else if (chain=="b"){
    ldgs = Bldgs
    mns_x = Bmns_x
    sds_x = Bsds_x
  } else {
    print("please specificy the 'chain' argument (a, b, or ab)")
  }
  rownames(ftz) = as.character(ftz$id)
  ftz = scale_variables(ftz, mns_x, sds_x)
  scores = as.matrix(ftz) %*% as.matrix(ldgs)
  rownames(scores) = rownames(ftz)
  scores[,1] = -scores[,1]
  scores[,3] = -scores[,3]
  colnames(scores) = c("TCRinnate", "TCR-8", "TCRmem", "TCRreg")
  return(data.frame(scores))
}
