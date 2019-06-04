# Two studies on how researchers fabricate data

This project contains two studies where we asked researchers to fabricate data in a controlled setting. This is part of [Chris Hartgerink's dissertation](https://github.com/chartgerink/dissertation) and was funded by the [Office of Research Integrity](https://ori.hhs.gov). 

In study 2 we also interviewed the participating researchers and transcribed their interviews. The transcripts can be found [here](https://doi.org/10.5281/zenodo.832490) and the project we're working on [here](https://github.com/chartgerink/2015ori-2).

## R-pkgs

```R
pROC
latex2exp

```

## Generate

```R
rmarkdown::render('redraft.Rmd', 'bookdown::pdf_document2')
```