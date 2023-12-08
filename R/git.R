git_clone_OptimalSPCA <- function() {
  stopifnot(
    system2('git', c('clone', 'https://github.com/ringw/Optimal-SPCA')) == 0)
  with_dir(
    'Optimal-SPCA',
    stopifnot(
      # TODO: Replace this branch name with a commit hash
      system2('git', c('checkout', 'deterministic')) == 0))
  'Optimal-SPCA'
}