
.onAttach <- function(lib, pkg)
{
  unlockBinding(".sFPCA", asNamespace("sFPCA")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    { # > figlet -f slant MCLUST
      packageStartupMessage(
"    __  ___________    __  _____________
   /  |/  / ____/ /   / / / / ___/_  __/
  / /|_/ / /   / /   / / / /\\__ \\ / /   
 / /  / / /___/ /___/ /_/ /___/ // /    
/_/  /_/\\____/_____/\\____//____//_/    version ", version)
}
else
  { packageStartupMessage("Package 'sFPCA' version ", version) } 

  invisible()
}

