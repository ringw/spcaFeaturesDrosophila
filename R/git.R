git_clone_OptimalSPCA <- function() {
  stopifnot(
    system2('git', c('clone', 'https://github.com/ringw/Optimal-SPCA')) == 0)
  with_dir(
    'Optimal-SPCA',
    stopifnot(
      # Branch name "deterministic" - includes SPCA.jl.
      system2('git', c('checkout', 'eb3f784cfdce82a6e03680584d7d0e331f2155af')) == 0))
  'Optimal-SPCA'
}