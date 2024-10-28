user_data should link to /data2/platform/cell_type_workspace
```bash
ln -s /data2/platform/cell_type_workspace/user_data user_data
```

Some R packages prep here:
```r
install.packages('https://mirror-hk.koddos.net/CRAN/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz', repos = NULL)
install.packages("fastmap", dependencies = TRUE, repos = "https://mirror-hk.koddos.net/CRAN/")

```