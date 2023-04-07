scale_variables <- function(data, mns, sds){
  data = data[,as.character(names(mns))]
  data = sweep(data, 2, mns)
  data = sweep(data, MARGIN=2, FUN="/", sds)
  data[is.na(data)] <- 0
  return(data)
}

add_adjacent_ints <- function(x, prefix){
  factors = paste("AF", seq(1, 5), sep="")
  factor_grid = expand.grid(factors, factors)
  print(factor_grid)
  for (i in 2:18){
    pos2 = paste(prefix, i, sep="")
    if (!(paste(pos2, "AF1", sep="_") %in% colnames(x))) { next }
    pos1 = paste(prefix, (i-1), sep="")
    for (j in 1:nrow(factor_grid)){
      ind2 = which(colnames(x)==paste(pos2, factor_grid$Var1[j], sep="_"))
      ind1 = which(colnames(x)==paste(pos1, factor_grid$Var2[j], sep="_"))
      new = x[,ind1] * x[,ind2]
      x = cbind(x, new)
      colnames(x)[ncol(x)] = paste(colnames(x)[ind1], colnames(x)[ind2], sep="_by_")
    }
  }
  return(x)
}

reformat_gene <- function(string, gene, brds){
  string = gsub("\\*", "-", string)
  pref = strsplit(string, gene)[[1]][1]
  info = strsplit(string, gene)[[1]][2]
  fam = strsplit(info, "-")[[1]][1]
  mem = strsplit(info, "-")[[1]][2]
  if (substr(fam, 1, 1)=="0"){
    fam = substr(fam, 2, nchar(fam))
  }
  if (!(is.na(mem))){
    if (substr(mem, 1, 1)=="0"){
      mem = substr(mem, 2, nchar(mem))
    }
  }
  base = paste(paste(gsub("C", "", pref), gene, sep=""), fam, sep="")
  if (!(base %in% brds) & !(is.na(mem))){
    return(paste(paste(gsub("C", "", pref), gene, sep=""), paste(fam, mem, sep="-"), sep=""))
  } else {
    return(base)
  }
}

create_hashmap <- function(key_vector, val_vector){
  h = hash()
  for (i in 1:length(key_vector)){
    h[[key_vector[i]]] = val_vector[i]
  }
  return(h)
}

str_index <- function(str_length, pos, scheme, max_length){
  if (scheme=="LR" & pos<=ceiling(str_length/2)){
    return(pos)
  }
  if (scheme=="LR" & pos>ceiling(str_length/2)){
    diff = max_length - str_length
    if ((pos-diff) > ceiling(str_length/2)){
      return(pos - diff)
    } else {
      return(".")
    }
  }
  if (scheme=="mid"){
    tab = floor((max_length - str_length)/2)
    if (tab>=pos){
      return(".")
    } else {
      if ((pos-tab)<=str_length){
        return(pos - tab)
      } else {
        return(".")
      }
    }
  }
}

get_feat_score <- function(x, amap){
  if (is.na(x)) { return(NA) }
  if (x=="."){
    return(NA)
  }
  sum = 0
  for (i in 1:nchar(x)){
    sum = sum + amap[[substr(x, i, i)]]
  }
  return(sum/nchar(x))
}

featurize_tcrs <- function(data, cdr3_align="mid", add_ints52 = FALSE, return_seq_grid=FALSE, beta_only=FALSE, do_jgenes = TRUE, restrict_length=TRUE){

  library(dplyr)
  library(stringr)
  library(hash)

  TRBVref$TRB_line <- NULL
  TRBVref$TRB_type <- NULL
  TRAVref$TRA_line <- NULL
  TRAVref$TRA_type <- NULL
  TRAJref = TRAJref[,c("TRA_Gene", "TRA_Jtail")]
  TRBJref = TRBJref[,c("TRB_Gene", "TRB_Jtail")]

  brd_BV = TRBVgrid$gene[!(grepl("-", TRBVgrid$gene))]
  brd_BJ = TRBJgrid$gene[!(grepl("-", TRBJgrid$gene))]
  brd_AV = TRAVgrid$gene[!(grepl("-", TRAVgrid$gene))]
  brd_AJ = TRAJgrid$gene[!(grepl("-", TRAJgrid$gene))]
  brds <- c(brd_BV, brd_BJ, brd_AV, brd_AJ)
  ##possibly reformat V and J genes
  print(nrow(data))
  print(head(data))
  data$TCRB_vgene = as.character(data$TCRB_vgene)
  data$TCRB_jgene = as.character(data$TCRB_jgene)
  data = data[data$TCRB_vgene!="unresolved",]
  if (do_jgenes){ data = data[data$TCRB_jgene!="unresolved",]; data = data[!(is.na(data$TCRB_jgene)),] }
  print(nrow(data))
  data = data[data$TCRB_vgene!="",]
  data = data[!(is.na(data$TCRB_vgene)),]
  print(nrow(data))
  print(table(data$TCRB_vgene))
  data$TCRB_vgene = sapply(as.character(data$TCRB_vgene), function(x) reformat_gene(x, gene="V", brds))
  print(table(data$TCRB_jgene))
  if (do_jgenes){
    data$TCRB_jgene = sapply(as.character(data$TCRB_jgene), function(x) reformat_gene(x, gene="J", brds))
    data = data[data$TCRB_jgene %in% TRBJgrid$gene,]
  }
  data = data[data$TCRB_vgene %in% TRBVgrid$gene,]
  print(nrow(data))
  data$TCRB_cdr3aa = as.character(data$TCRB_cdr3aa)
  if (restrict_length){
    nc = sapply(data$TCRB_cdr3aa, function(x) nchar(x))
    data = data[nc >=11 & nc <= 18,]
  }
  if (!(beta_only)){
    data$TCRA_vgene = as.character(data$TCRA_vgene)
    data$TCRA_jgene = as.character(data$TCRA_jgene)
    data = data[data$TCRA_vgene!="unresolved",]
    data = data[data$TCRA_vgene!="",]
    if (do_jgenes) { data = data[data$TCRA_jgene!="unresolved",] }
    data$TCRA_vgene = sapply(as.character(data$TCRA_vgene), function(x) reformat_gene(x, gene="V", brds))
    data = data[data$TCRA_vgene %in% TRAVgrid$gene,]
    if (do_jgenes){
      data$TCRA_jgene = sapply(as.character(data$TCRA_jgene), function(x) reformat_gene(x, gene="J", brds))
      data = data[data$TCRA_jgene %in% TRAJgrid$gene,]
    }
    data$TCRA_cdr3aa = as.character(data$TCRA_cdr3aa)
    if (restrict_length){
      nc = sapply(data$TCRA_cdr3aa, function(x) nchar(x))
      data = data[nc >=10 & nc <= 17,]
    }
  }

  if (!(beta_only)){
    pos = left_join(data, TRAVgrid, by=c("TCRA_vgene"="gene"))
    pos = left_join(pos, TRAJgrid, by=c("TCRA_jgene"="gene"))
    pos = left_join(pos, TRBVgrid, by=c("TCRB_vgene"="gene"))
  } else {
    pos = left_join(data, TRBVgrid, by=c("TCRB_vgene"="gene"))
  }
  pos = left_join(pos, TRBJgrid, by=c("TCRB_jgene"="gene"))
  ##now need to think about LR or mid alignment
  print(head(data))
  nc = sapply(data$TCRB_cdr3aa, function(x) nchar(x))
  print(min(nc))
  print(max(nc))
  print("any NAs?")
  print(table(is.na(data$TCRB_cdr3aa)))
  print(table(is.na(nc)))
  if (!(beta_only)){
    for (i in 1:17){
      pos$new = sapply(data$TCRA_cdr3aa, function(y) substr(y, str_index(nchar(y), i, cdr3_align, 17), str_index(nchar(y), i, cdr3_align, 17)))
      pos$new[is.na(pos$new)] = "."
      colnames(pos)[ncol(pos)] = paste("TRAcdr3_p", i, sep="")
    }
  }
  for (i in 1:18){
    pos$new = sapply(data$TCRB_cdr3aa, function(y) substr(y, str_index(nchar(y), i, cdr3_align, 18), str_index(nchar(y), i, cdr3_align, 18)))
    pos$new[is.na(pos$new)] = "."
    colnames(pos)[ncol(pos)] = paste("TRBcdr3_p", i, sep="")
  }
  if (return_seq_grid){
    ind = which(colnames(pos)=="TRA_p1")
    pos = pos[,ind:ncol(pos)]
    rownames(pos) = as.character(data[,1])
    return(pos)
  }

  if (!(beta_only)){
    loops = left_join(data, TRAVref, by=c("TCRA_vgene"="gene"))
    loops = left_join(loops, TRAJref, by=c("TCRA_jgene"="TRA_Gene"))
    loops = left_join(loops, TRBVref, by=c("TCRB_vgene"="gene"))
    loops = left_join(loops, TRBJref, by=c("TCRB_jgene"="TRB_Gene"))
  } else {
    loops = left_join(data, TRBVref, by=c("TCRB_vgene"="gene"))
    loops = left_join(loops, TRBJref, by=c("TCRB_jgene"="TRB_Gene"))
  }
  prc_fts = data.frame(id = as.character(data[,1]))
  ##we don't calculate G because this is our reference
  aminos <- c("A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  if (!(beta_only)){
    ind = which(colnames(loops)=="TRA_FR1")
    segs = c(colnames(loops)[ind:ncol(loops)], "TCRA_cdr3aa", "TCRB_cdr3aa")
  } else {
    ind = which(colnames(loops)=="TRB_FR1")
    segs = c(colnames(loops)[ind:ncol(loops)], "TCRB_cdr3aa")
  }
  for (i in 1:length(segs)){
    ind = which(colnames(loops)==segs[i])
    loops[,ind] = as.character(loops[,ind])
    prc_fts$new = sapply(loops[,ind], function(x) nchar(x) - str_count(x, "\\."))
    colnames(prc_fts)[ncol(prc_fts)]  = paste(segs[i], "length", sep="_")
    for (j in 1:length(aminos)){
      prc_fts$new = sapply(loops[,ind], function(x) str_count(x, aminos[j])/(nchar(x)-str_count(x, "\\.")))
      colnames(prc_fts)[ncol(prc_fts)] = paste(segs[i], paste("prc", aminos[j], sep=""), sep="_")
    }
  }
  atch_maps <- list()
  for (i in 2:ncol(atch)){
    atch_maps[[i-1]] <- create_hashmap(c(as.character(atch$symbol), "."), c(as.numeric(as.character(atch[,i])), NA))
  }

  atch_fts = data.frame(id = as.character(data[,1]))
  if (!(beta_only)){
    ind = which(colnames(pos)=="TRA_p1")
  } else {
    ind = which(colnames(pos)=="TRB_p1")
  }

  for (i in ind:ncol(pos)){
    for (j in 1:5){
      atch_fts$new = sapply(pos[,i], function(x) get_feat_score(x, atch_maps[[j]]))
      colnames(atch_fts)[ncol(atch_fts)] = paste(colnames(pos)[i], paste("AF", j, sep=""), sep="_")
    }
  }
  if (!(beta_only)){
    res = data.frame(cbind(prc_fts[,1:21], atch_fts[,2:which(colnames(atch_fts)=="TRA_p26_AF5")],
                           prc_fts[,22:41], atch_fts[which(colnames(atch_fts)=="TRA_p27_AF1"):which(colnames(atch_fts)=="TRA_p38_AF5")],
                           prc_fts[,42:61], atch_fts[which(colnames(atch_fts)=="TRA_p39_AF1"):which(colnames(atch_fts)=="TRA_p55_AF5")],
                           prc_fts[,62:81], atch_fts[which(colnames(atch_fts)=="TRA_p56_AF1"):which(colnames(atch_fts)=="TRA_p65_AF5")],
                           prc_fts[,82:101], atch_fts[which(colnames(atch_fts)=="TRA_p66_AF1"):which(colnames(atch_fts)=="TRA_p103_AF5")],
                           prc_fts[,242:261], atch_fts[which(colnames(atch_fts)=="TRAcdr3_p1_AF1"):which(colnames(atch_fts)=="TRAcdr3_p17_AF5")],
                           prc_fts[,102:121], atch_fts[which(colnames(atch_fts)=="TRA_p119_AF1"):which(colnames(atch_fts)=="TRA_p128_AF5")],

                           prc_fts[,122:141], atch_fts[,which(colnames(atch_fts)=="TRB_p1_AF1"):which(colnames(atch_fts)=="TRB_p26_AF5")],
                           prc_fts[,142:161], atch_fts[which(colnames(atch_fts)=="TRB_p27_AF1"):which(colnames(atch_fts)=="TRB_p38_AF5")],
                           prc_fts[,162:181], atch_fts[which(colnames(atch_fts)=="TRB_p39_AF1"):which(colnames(atch_fts)=="TRB_p55_AF5")],
                           prc_fts[,182:201], atch_fts[which(colnames(atch_fts)=="TRB_p56_AF1"):which(colnames(atch_fts)=="TRB_p65_AF5")],
                           prc_fts[,202:221], atch_fts[which(colnames(atch_fts)=="TRB_p66_AF1"):which(colnames(atch_fts)=="TRB_p103_AF5")],
                           prc_fts[,262:281], atch_fts[which(colnames(atch_fts)=="TRBcdr3_p1_AF1"):which(colnames(atch_fts)=="TRBcdr3_p18_AF5")],
                           prc_fts[,222:241], atch_fts[which(colnames(atch_fts)=="TRB_p119_AF1"):which(colnames(atch_fts)=="TRB_p127_AF5")]))
  } else {
    res = data.frame(cbind(prc_fts[,1:21], atch_fts[,2:which(colnames(atch_fts)=="TRB_p26_AF5")],
                           prc_fts[,22:41], atch_fts[which(colnames(atch_fts)=="TRB_p27_AF1"):which(colnames(atch_fts)=="TRB_p38_AF5")],
                           prc_fts[,42:61], atch_fts[which(colnames(atch_fts)=="TRB_p39_AF1"):which(colnames(atch_fts)=="TRB_p55_AF5")],
                           prc_fts[,62:81], atch_fts[which(colnames(atch_fts)=="TRB_p56_AF1"):which(colnames(atch_fts)=="TRB_p65_AF5")],
                           prc_fts[,82:101], atch_fts[which(colnames(atch_fts)=="TRB_p66_AF1"):which(colnames(atch_fts)=="TRB_p103_AF5")],
                           prc_fts[,122:141], atch_fts[which(colnames(atch_fts)=="TRBcdr3_p1_AF1"):which(colnames(atch_fts)=="TRBcdr3_p18_AF5")],
                           prc_fts[,102:121], atch_fts[which(colnames(atch_fts)=="TRB_p119_AF1"):which(colnames(atch_fts)=="TRB_p127_AF5")]))
  }
  if (add_ints52){
    res = add_adjacent_ints(res, "TRAcdr3_p")
    res = add_adjacent_ints(res, "TRBcdr3_p")
  }
  return(res)
}

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
  ftz = scale_variables(ftz, mns_x, sds_x)
  scores = as.matrix(ftz) %*% as.matrix(ldgs)
  rownames(scores) = rownames(ftz)
  scores[,1] = -scores[,1]
  scores[,3] = -scores[,3]
  colnames(scores) = c("TCRinnnate", "TCR-8", "TCRmem", "TCRreg")
  return(scores)
}


