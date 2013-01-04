### This file contains functions for "seq-gen".

seqgen <- function(opts = NULL, rooted.tree = NULL, newick.tree = NULL,
    input = NULL){
  argv <- "seq-gen"
  temp.file.seqgen <- tempfile("seqgen.")
  temp.file.ms <- tempfile("ms.")

  if(! is.null(opts)){
    if((! is.null(rooted.tree)) && (class(rooted.tree) == "phylo")){
      newick.tree <- write.tree(rooted.tree, digits = 12)
    }
    if((!is.null(newick.tree)) && (!is.null(input))){
      stop("rooted.tree/newick.tree and input can not work at the same time.")
    }
    if(! is.null(newick.tree)){
      write(newick.tree, file = temp.file.ms, sep = "")
    } else if(! is.null(input)){
      write(input, file = temp.file.ms, sep = "\n")
    } else{
      stop("A newick or rooted/phylo tree is required.")
    }

    argv <- c(argv, unlist(strsplit(opts, " ")), temp.file.ms)
    .Call("R_seq_gen_main", argv, temp.file.seqgen, PACKAGE = "phyclust")

    ret <- scan(file = temp.file.seqgen,
                what = "character", sep = "\n", quiet = TRUE)
    class(ret) <- "seqgen"

    unlink(temp.file.seqgen)
    unlink(temp.file.ms)
    return(ret)
  }

  argv <- c(argv, "-h")
  try(.Call("R_seq_gen_main", argv, temp.file.ms, PACKAGE = "phyclust"),
      silent = TRUE)
  unlink(temp.file.seqgen)
  unlink(temp.file.ms)
  invisible()
} # End of seqgen().

print.seqgen <- function(x, ...){
  seqgen <- x
  cat(seqgen, sep = "\n")
} # Enf of print.seqgen().

