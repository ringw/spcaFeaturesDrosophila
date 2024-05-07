julia_env <- basilisk::BasiliskEnvironment(
  envname='julia_env',
  pkgname='spcaFeaturesDrosophila',
  packages='julia=1.9.4',
  channels='conda-forge'
)

scikit_learn_env <- basilisk::BasiliskEnvironment(
  envname='scikit_learn_env',
  pkgname='spcaFeaturesDrosophila',
  packages=c('python=3.10.10', 'scikit-learn=1.3.0')
)