julia_pkg_OptimalSPCADepot <- function() {
  dir.create('Optimal-SPCA-Depot')
  with_envvar(
    c(
      JULIA_DEPOT_PATH=paste0(getwd(), '/Optimal-SPCA-Depot'),
      LD_LIBRARY_PATH=paste0(obtainEnvironmentPath(julia_env), '/lib')
    ),
    stopifnot(
    system2(
      paste0(
        obtainEnvironmentPath(julia_env),
        '/bin/julia'
      ),
      c(
        '--project=.',
        '-e',
        '"import Pkg"',
        '-e',
        '\'Pkg.add(["Arpack","CSV","DataFrames","StatsBase","Tables"])\''
      )
    )
    == 0)
  )
}