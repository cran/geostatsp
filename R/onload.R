.onLoad <- function(libname, pkgname) {
  # During R CMD check, _R_CHECK_LIMIT_CORES_ is typically set.
  # Cap OpenMP threads to avoid example/test failures on constrained builders.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))) {
    Sys.setenv(
      OMP_NUM_THREADS = "2",
      OMP_THREAD_LIMIT = "2",
      KMP_DEVICE_THREAD_LIMIT = "2",
      KMP_TEAMS_THREAD_LIMIT = "2"
    )
  }
}
