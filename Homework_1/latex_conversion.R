# ============================================================
# latex_export.R
# Formatting utility only.
# Converts already-computed Problem 7 results into LaTeX.
# ============================================================

outfile <- "problem7_results.tex"


# ------------------------------------------------------------
# Helper: format a number
# ------------------------------------------------------------
latex_num <- function(x, digits = 10) {

  if (abs(x) < 1e-4 && x != 0) {
    return(format(x, scientific = TRUE, digits = digits))
  }

  format(x, scientific = FALSE, digits = digits, trim = TRUE)
}


# ------------------------------------------------------------
# Helper: convert vector to LaTeX column vector
# ------------------------------------------------------------
latex_vector <- function(x, digits = 10) {

  values <- sapply(x, latex_num, digits = digits)

  paste0(
    "\\begin{bmatrix}\n",
    paste(values, collapse = " \\\\\n"),
    "\n\\end{bmatrix}"
  )
}


# ------------------------------------------------------------
# Helper: convert matrix to LaTeX matrix
# ------------------------------------------------------------
latex_matrix <- function(A, digits = 10) {

  rows <- apply(A, 1, function(row) {

    values <- sapply(row, latex_num, digits = digits)

    paste(values, collapse = " & ")
  })

  paste0(
    "\\begin{bmatrix}\n",
    paste(rows, collapse = " \\\\\n"),
    "\n\\end{bmatrix}"
  )
}


# ------------------------------------------------------------
# Helper: format singular values
# ------------------------------------------------------------
latex_singular_values <- function(x, digits = 10) {

  values <- sapply(x, latex_num, digits = digits)

  paste(values, collapse = ",\\; ")
}

# ------------------------------------------------------------
# Helper: format scientific notation
# ------------------------------------------------------------
latex_num <- function(x, digits = 10) {

    if(x == 0){
        return("0")
    }

    if (abs(x) < 1e-4 || abs(x) > 1e5) {
        s <- format(x, scientific = TRUE, digits = digits)

        parts <- strsplit(s, "e")[[1]]

        coefficient <- parts[1]
        exponent <- as.integer(parts[2])

        return(
            paste0(coefficient, "\\times10^{", exponent, "}")
        )
    }

    format(x, scientific = FALSE, digits = digits, trim = TRUE)
}

# ------------------------------------------------------------
# Write one case
# ------------------------------------------------------------
write_case <- function(file, case) {

  cat(
    paste0("\\subsubsection*{Case ", case$name, "}\n\n"),
    file = file,
    append = TRUE
  )


  cat("\\begin{itemize}\n", file = file, append = TRUE)


  cat(
    paste0(
      "\\item Dimensions: $",
      case$m,
      " \\times ",
      case$n,
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item System type: ",
      case$system_type,
      ".\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Singular values: $",
      latex_singular_values(case$singular_values),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Numerical-rank tolerance: $\\tau_G = ",
      latex_num(case$tau),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Numerical rank: $r_G = ",
      case$rank,
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    "\\item Range projector:\n",
    file = file,
    append = TRUE
  )

  cat(
    "\\[\n",
    "P_{\\mathrm{Range}(G)} = ",
    latex_matrix(case$P),
    "\n\\]\n",
    file = file,
    append = TRUE,
    sep = ""
  )


  cat(
    paste0(
      "\\item Symmetry check: $\\|P-P^T\\|_F = ",
      latex_num(case$symmetry_error),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Idempotence check: $\\|P^2-P\\|_F = ",
      latex_num(case$idempotence_error),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    "\\item Projected vector:\n",
    file = file,
    append = TRUE
  )

  cat(
    "\\[\n",
    "P_{\\mathrm{Range}(G)}f = ",
    latex_vector(case$Pf),
    "\n\\]\n",
    file = file,
    append = TRUE,
    sep = ""
  )


  cat(
    "\\item Residual vector:\n",
    file = file,
    append = TRUE
  )

  cat(
    "\\[\n",
    "r_f = ",
    latex_vector(case$rf),
    "\n\\]\n",
    file = file,
    append = TRUE,
    sep = ""
  )


  cat(
    paste0(
      "\\item Residual norm: $\\|r_f\\|_2 = ",
      latex_num(case$residual_norm),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Range-membership threshold: $",
      latex_num(case$threshold),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Range membership: ",
      case$range_statement,
      ".\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Solution type: ",
      case$solution_type,
      ".\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    "\\item Computed solution:\n",
    file = file,
    append = TRUE
  )

  cat(
    "\\[\n",
    "w = ",
    latex_vector(case$solution),
    "\n\\]\n",
    file = file,
    append = TRUE,
    sep = ""
  )


  cat(
    paste0(
      "\\item Solution norm: $\\|w\\|_2 = ",
      latex_num(case$solution_norm),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat(
    paste0(
      "\\item Final solution residual norm: $",
      latex_num(case$solution_residual),
      "$.\n"
    ),
    file = file,
    append = TRUE
  )


  cat("\\end{itemize}\n\n", file = file, append = TRUE)
}


# ============================================================
# RESULTS
#
# Paste values already computed by your Problem 7 script here.
# This section does NOT perform any linear algebra.
# ============================================================

case_1 <- list(

  name = "1",

  m = 3,
  n = 4,

  system_type = "Underdetermined",

  singular_values = c(
    4.421502343,
    2.020442439,
    0.5973737652
  ),

  tau = 4.421502343e-7,

  rank = 3,

  P = matrix(
    c(
      1,                 2.775557562e-17, -6.938893904e-18,
      2.775557562e-17,   1,                1.804112415e-16,
     -6.938893904e-18,   1.804112415e-16,  1
    ),
    nrow = 3,
    byrow = TRUE
  ),

  symmetry_error = 0,

  idempotence_error = 1.161306352e-15,

  Pf = c(
     0.61339829,
     0.11012282,
    -0.06426754
  ),

  rf = c(
    -4.440892099e-16,
     9.714451465e-17,
     0
  ),

  residual_norm = 4.545902364e-16,

  threshold = 1e-7,

  range_statement = "$f \\in \\mathrm{Range}(G)$",

  solution_type = "Infinitely many exact solutions",

  solution = c(
    -0.07025444832,
     0.1984662065,
     0.01998729011,
     0.07207220944
  ),

  solution_norm = 0.2234242999,

  solution_residual = 2.407706197e-16
)


# ============================================================
# WRITE FILE
# ============================================================

# Clear the file ONCE
cat("", file = outfile)


# Write cases
write_case(outfile, case_1)


# ============================================================
# CONFIRM
# ============================================================

cat("\nLaTeX export complete.\n")
cat("File:", normalizePath(outfile), "\n\n")

cat("Contents written:\n")
cat("----------------------------------------\n")
cat(readLines(outfile), sep = "\n")