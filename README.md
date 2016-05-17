# idfgen
Import, Manipulate, and Export Genotype Data

`idfgen` is an R package to import biological and genotype data from a Progeny database or www.fishgen.net. Once data is imported
into `idfgen` it can then be manipulated using a variety of functions prior to exporting the data in a useful format (e.g., for 
Structure, genepop, GenAlEx, gsi_sim, snppit, Colony).

To install `idfgen` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:
```
install.packages("devtools")
library(devtools)
```
Once `devtools` is successfully installed, use the following to install `idfgen`:
```
devtools::install_github("mackerman44/idfgen")
```
Finally, to access all of the functions available within `idfgen`, please run:
```
devtools::load_all()
```
If you are interested in making contributions, consider getting a GitHub account, fork this repository, modify, and send a pull request.

Enjoy!
