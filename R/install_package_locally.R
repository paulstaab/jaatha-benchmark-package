#' @export
installPackageLocally = function(package.path, install.path=tempdir()) {
  install.path <- tempfile('tmp_repo_', tmpdir=install.path)
  print(install.path)
  dir.create(install.path, recursive=TRUE, showWarnings=FALSE)
  install.packages(package.path, dependencies=TRUE, 
                   lib=install.path, repos = NULL) 
  for (package in list.dirs(install.path, full.names=FALSE, recursive=FALSE)) {
    library(package, lib.loc=install.path, character.only=TRUE)
  }
}