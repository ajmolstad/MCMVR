# MCMVR
An R package implementing the method proposed in [An explicit mean-covariance parameterization for multivariate response linear regression](https://www.tandfonline.com/doi/abs/10.1080/10618600.2020.1853551?casa_token=dQzCJAFc1ZoAAAAA%3AUaq0GRdBijyS7kavHT9njRKCFqCvnE-XBddXiI_w8BAEf0ZCllJVy_ALwrcXpGxSJSKcdS4i7P_q&journalCode=ucgs20). Please contact [amolstad@umn.edu](mailto:amolstad@umn.edu) with any comments or questions. 

### Installation
MCMVR can be loaded directly into R through the the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/MCMVR")
```
### Citation instructions
Please cite the most recent version of the article mentioned above. As of June 2023, this was the following (in bibtex): 
```
@article{Molstad2021Explicit,
  author = {Aaron J. Molstad and Guangwei Weng and Charles R. Doss and Adam J. Rothman},
  title = {An Explicit Mean-Covariance Parameterization for Multivariate Response Linear Regression},
  journal = {Journal of Computational and Graphical Statistics},
  volume = {30},
  number = {3},
  pages = {612-621},
  year  = {2021},
  publisher = {Taylor & Francis},
  doi = {10.1080/10618600.2020.1853551}
}
```
### Usage directions
Please visit [this example page](https://ajmolstad.github.io/docs/MCMVR_Example.html) for details on implementation and usage. 
