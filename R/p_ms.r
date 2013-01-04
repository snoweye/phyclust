### This file contains functions for "ms".

ms <- function(nsam = NULL, nreps = 1, opts = NULL, temp.file = NULL){
  if(is.null(temp.file)){
    temp.file.ms <- tempfile("ms.")
  } else{
    temp.file.ms <- temp.file
  }

  if(! is.null(opts) && ! is.null(nsam)){
    if(nsam >= 2){
      nsam <- as.character(nsam)
      nreps <- as.character(nreps)
      argv <- c("ms", nsam, nreps, unlist(strsplit(opts, " ")))

      .Call("R_ms_main", argv, temp.file.ms, PACKAGE = "phyclust")
      ret <- scan(file = temp.file.ms,
                  what = "character", sep = "\n", quiet = TRUE)
      class(ret) <- "ms"

      if(is.null(temp.file)){
        unlink(temp.file.ms)
      }
      return(ret)
    }
  }

  argv <- c("ms", "-h")
  try(.Call("R_ms_main", argv, temp.file.ms, PACKAGE = "phyclust"),
      silent = TRUE)
  unlink(temp.file.ms)
  invisible()
} # End of ms().

print.ms <- function(x, ...){
  ms <- x
  cat(ms, sep = "\n")
} # End of print.ms().
