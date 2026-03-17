if (!exists("HamiltonianMultiplet", inherits = TRUE)) {
  local_class_path <- c(
    "R/HamiltonianMultiplet.R",
    file.path("..", "..", "R", "HamiltonianMultiplet.R")
  )
  local_class_path <- local_class_path[file.exists(local_class_path)][1]

  if (!is.na(local_class_path)) {
    if (!exists("R6Class", inherits = TRUE)) {
      assign("R6Class", get("R6Class", envir = asNamespace("R6")), envir = globalenv())
    }
    source(local_class_path, local = globalenv())
  } else if (requireNamespace("fitnmr", quietly = TRUE)) {
    assign(
      "HamiltonianMultiplet",
      get("HamiltonianMultiplet", envir = asNamespace("fitnmr")),
      envir = globalenv()
    )
  }
}
