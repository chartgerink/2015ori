#' Digit analysis
#'
#' @param x A vector or matrix of numeric values.
#' @param type Type of digit analysis ('benford' or 'terminal')
#'
#' @return
#' @export
#'
#' @examples digit_analysis(c(1.234, 65.4321, 53.222), type = 'terminal')

digit_analysis <- function(x, type = 'terminal') {
  if(!(is.matrix(x) | is.vector(x))) {
    stop('Please only specify a vector or matrix. 
         If specifying a matrix, ensure results from one set of digits are in the 
         rows!')
  }
  
  if(!(type == 'terminal' | type == 'benford')) {
    stop("Only 'benford' and 'terminal' allowed as types.")
  }
  
  if (type == 'terminal') {
    df <- 10 - 1 
    
    # x <- decimator(x)
    
    matrix_of_digits <- apply(
      t(x),
      1,
      function(y)
      {
        splitted <- strsplit(as.character(y), "")
        digits <- lapply(splitted, function(z) tail(z, 1))
        return(unlist(digits))
      })
  } else if (type == 'benford')
  { 
    df <- 9 - 1
    
    matrix_of_digits <- apply(
      t(x),
      1,
      function(y)
      {
        splitted <- strsplit(as.character(y), "")
        digits <- lapply(splitted, function(z) head(z, 1))
        return(unlist(digits))
      })
  } else
  {
    stop("Something went awry.")
  }
  
  mode(matrix_of_digits) <- "numeric"
  
  chi <- apply(matrix_of_digits,
    2, function(q) {
      obs <- digit_counter(q, type = type)
      exp <- expected_digit_counter(q, type = type)
      testval <- sum((obs-exp)^2/exp)
      return(testval)
    })
  
  pval <- pchisq(q = chi, df = df,
    lower.tail = FALSE)
  res <- data.frame(chi, df = df, 
                    pval)
  
  return(res)
}

digit_counter <- function(x, type) {
  if(!is.vector(x)) stop("Currently only works with vectors.")
  
  x <- decimator(x)
  if (type == 'terminal') {
    counts <- rep(0, 11)
    names(counts) <- seq(0, 10)
    digits <- table(regmatches(x, regexpr("\\d$", x)))
    counts[names(digits)] <- digits
  } else if (type == 'benford') {
    counts <- rep(0, 10)
    names(counts) <- seq(1, 10)
    digits <- table(regmatches(x, regexpr("\\d", x)))
    counts[names(digits)] <- digits
  } else {
    stop("ERROR")
  }
  
  return(counts)
}

expected_digit_counter <- function(x, type = type) {
  if (!is.vector(x)) stop("Currently only works with vectors.")
  if (!(type == 'terminal' | type == 'benford')) {
    stop("Only benford and terminal allowed as types.")
  }
  
  if (type == 'terminal') {
    cell <- length(x) / 10
    counts <- rep(cell, 10)
  } else if (type == 'benford') {
    d <- 1:9
    counts <- length(x) * log(x = ((d + 1) / d), base = 10)
  } else {
    stop("Something went awry.")
  }
  
  return(counts)
}

decimator <- function(x) {
  decimal <- regexpr('\\.', x)  
  length <- nchar(x)

  decimated <- 10 ^ (length - decimal) * x
  # Make sure integers remain
  decimated[decimal == -1] <- x[decimal == -1]
  
  return(decimated)
}