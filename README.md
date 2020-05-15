
title: “Improving DIABLO” author: “Amrit Singh” date: “15 May, 2020”

## What should this new method include?

  - n \<\< p (DIABLO )
  - missing values (DIABLO X)
  - integration of multiple datasets/blocks (DIABLO )
  - model the hierarchy between datasets (DIABLO and X, no
    directionality)
  - incorporate observational weights (DIABLO )
  - regression and discriminant analysis model (DIABLO and X, not
    implemented for regression)

# Multi-block Partial Least Squares (MB-PLS)

``` r
library(ade4)
data(chickenk)
Mortality <- chickenk[[1]]
dudiY.chick <- dudi.pca(Mortality, center = TRUE, scale = TRUE, scannf =
FALSE)
ktabX.chick <- ktab.list.df(chickenk[2:5])
resmbpls.chick <- mbpls(dudiY.chick, ktabX.chick, scale = TRUE,
option = "uniform", scannf = FALSE)
```

## Block importance

``` r
setNames(resmbpls.chick$bip, c("t1", "t2")) %>% 
  as.data.frame() %>% 
  mutate(block_name = rownames(.)) %>% 
  gather(Comp, value, -block_name) %>% 
  ggplot(aes(x = Comp, y = value)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~block_name)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Global component plot

### Xs and Ys

``` r
global_comps <- cbind(resmbpls.chick$lY[, 1, drop = FALSE], resmbpls.chick$lX[, 1, drop = FALSE], resmbpls.chick$lY[, 2, drop = FALSE], resmbpls.chick$lX[, 2, drop = FALSE])
colnames(global_comps) <- c("u1", "t1", "u2", "t2")

pairs(global_comps)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Individual component plots

``` r
indx <- resmbpls.chick$TlX %>% rownames(.) %>% as.numeric(.)
start <- which(indx == 1)
end <- which(indx == nrow(Mortality))

indComps <- mapply(function(start, end){
  resmbpls.chick$TlX[start:end, ]
}, start = start, end = end, SIMPLIFY = FALSE)

indComps <- lapply(seq_len(length(indComps)), function(i){
  x <- indComps[[i]]
  colnames(x) <- paste(c("t1", "t2"), i, sep = "_")
  x
})

par(mfrow = c(2, 2))
sapply(indComps, function(comps){ plot(comps)})
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ## NULL
    ## 
    ## [[3]]
    ## NULL
    ## 
    ## [[4]]
    ## NULL

## Variable importance

``` r
resmbpls.chick$vip %>% 
  as.data.frame() %>% 
  mutate(variable_name = rownames(.),
    Dataset = rep(paste0("Dataset", 1:(length(chickenk) - 1)), sapply(chickenk[2:5], ncol))) %>% 
  rename(t1 = Ax1, t2 = Ax2) %>% 
  gather(comp, value, c("t1", "t2")) %>% 
  ggplot(aes(x = comp, y = value,fill = Dataset )) +
  geom_bar(stat = "identity") +
  facet_wrap(Dataset~variable_name)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## compare to JIVE and MOFA
