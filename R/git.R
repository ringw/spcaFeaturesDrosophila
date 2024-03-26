git_clone_OptimalSPCA <- function() {
  stopifnot(
    system2('git', c('clone', 'https://github.com/ringw/Optimal-SPCA')) == 0)
  with_dir(
    'Optimal-SPCA',
    stopifnot(
      # Branch name "deterministic" - includes SPCA.jl.
      system2('git', c('checkout', 'e0a13c3ee559108b0d867ac5523e5c57471f88b3')) == 0))
  'Optimal-SPCA'
}