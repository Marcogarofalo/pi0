args <- commandArgs(trailingOnly = TRUE)

# if (length(args)<1){
#   stop("usage: Rscript convert.R  infile", call.=FALSE)
# }
# args[1] <- "data/B64/merging_info.R"
args[1] <- "data/B64/merging_info_conn.R"
source(args[1])
# L<-as.integer(args[2])

# out<-args[1]

# out<-"data/B64/B64.txt"
# L<- as.integer("64")
# args[3]="data/B64/tri_1_light_defl_nev400_S1.txt"
# args[4]="data/B64/tri_1_light_defl_nev400_S2.txt"


T <- L * 2
l <- list()
Nin <- length(infiles)

for (i in seq_along(infiles)) {
  print(infiles[i])
  df <- read.table(infiles[i], comment.char = "", fill = TRUE)
  l[[i]] <- df
}

Nconf <- list()
confs <- list()
# data<-list()

for (i in c(1:Nin)) {
  Nconf[[i]] <- length(l[[i]][, 1]) / (L + 1 + 1)
  print(Nconf[[i]])
  lc <- which(l[[i]][, 1] == "#")
  confs[[i]] <- l[[i]][lc, 2]
  # data[[i]]<-l[[i]][-lc,2]
}

size <- ncorr * 2 * T

conf_all_the_same <- TRUE
for (i in seq_along(confs)) {
  if (!all(confs[[1]] == confs[[i]])) {
    cat("conf of file ", i, " differs\n")
    conf_all_the_same <- FALSE
  }
}

# fileConn<-file(out)
# options(scipen = 1)
# options(digits=12)




#
# if (conf_all_the_same) {
print("merging")
Njack <- 0
for (i in seq_along(Nconf)) {
  Njack <- Njack + Nconf[[1]]
}
con <- file(out, "wb")
Sint <- 4
writeBin(as.integer(Njack), con, size = Sint, endian = "little")
writeBin(as.integer(T), con, size = Sint, endian = "little")
writeBin(as.integer(L), con, size = Sint, endian = "little")
writeBin(as.integer(ncorr), con, size = Sint, endian = "little")

Sdoub <- 8
writeBin(as.numeric(beta), con, size = Sdoub, endian = "little")
writeBin(as.numeric(kappa), con, size = Sdoub, endian = "little")

writeBin(length(mus), con, size = Sint, endian = "little")
for (mu in mus) {
  writeBin(as.numeric(mu), con, size = Sdoub, endian = "little")
}
writeBin(length(rs), con, size = Sint, endian = "little")
for (r in rs) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(length(thetas), con, size = Sint, endian = "little")
for (r in thetas) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(length(gammas), con, size = Sint, endian = "little")
for (r in gammas) {
  writeChar(r, con, nchars = nchar(r, type = "chars"), eos = "", useBytes = FALSE)
}
writeBin(length(smearing), con, size = Sint, endian = "little")
for (r in smearing) {
  writeChar(r, con, nchars = nchar(r, type = "chars"), eos = "", useBytes = FALSE)
}
writeBin(length(bananas), con, size = Sint, endian = "little")
for (r in bananas) {
  writeBin(as.integer(r), con, size = Sint, endian = "little")
}
writeBin(length(oranges), con, size = Sint, endian = "little")
for (r in oranges) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(as.integer(size), con, size = Sint, endian = "little")

#   sink(out)
count <- 0
for (i in seq_along(confs[[1]])) {
  #     cat("\n")
  #     cat(paste("#",confs[[1]][i]),"\n")
  #     cat("\n")
  for (j in c(1:Nin)) {
    writeBin(as.integer(count), con, size = Sint, endian = "little")
    for (t in c(1:(T / 2 ))) {
      
      r <-c( as.numeric(l[[j]][t + i + (T / 2 + 1) * (i - 1), 1]), 0.0)
      writeBin(r, con, size = Sdoub, endian = "little")
    }
    for (t in c((T / 2 + 1):T)) {
      
      r <- c(as.numeric(l[[j]][ (T  +1- t) + i + (T / 2 + 1) * (i - 1), 1]),0.0)
      writeBin(r, con, size = Sdoub, endian = "little")
      # cat(t, "  ", r, "  ",  i, " \n")
    }
    # r <- r / Nin
    count <- count + 1
    #       cat(formatC(r, format = "e", digits = 12),"\n")
  }
  #
}
#   sink()
close(con)
# } else {
#   stop("conf are not all the same, nothing impremented here", call. = FALSE)
# }
