# Looking into the data

Reading the data

``` r
library(data.table)
library(ggplot2)

simstats <- fread("data/02-dataprep-network-stats.csv.gz")
```

## Visual inspection

``` r
ggplot(simstats[final_preval > 0]) +
  geom_histogram(aes(x = final_preval)) +
  scale_y_log10() +
  labs(
    x = "Final prevalence",
    y = "Outbreak size\n(log10)",
  )
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    Warning: Transformation introduced infinite values in continuous y-axis

    Warning: Removed 9 rows containing missing values (`geom_bar()`).

![](03-report_files/figure-commonmark/unnamed-chunk-2-1.png)

``` r
ggplot(simstats) +
  geom_point(aes(x = rt, y = ergm_triangle)) 
```

    Warning: Removed 7 rows containing missing values (`geom_point()`).

![](03-report_files/figure-commonmark/unnamed-chunk-2-2.png)

``` r
ggplot(simstats) +
  geom_point(aes(x = rt, y = final_preval)) +
  scale_y_log10() +
  labs(
    x = "Reproduction number",
    y = "Outbreak size\n(log10)",
  )
```

    Warning: Removed 7 rows containing missing values (`geom_point()`).

![](03-report_files/figure-commonmark/unnamed-chunk-2-3.png)

``` r
# A group of figures showing the outbreak size as 
# a function of the igraph_* variables in the dataset simstats
igraph_vars <- grep("^igraph", names(simstats), value = TRUE)

igraph_vars <- setdiff(
  igraph_vars, c(
    "igraph_diameter",
    "igraph_components"
    ))

graphs <- list()
for (igraph_var in igraph_vars) {
  graphs[[igraph_var]] <-
    ggplot(simstats) +
      geom_point(aes_string(x = "final_preval", y = igraph_var)) +
      scale_y_log10() +
      labs(
        x = "Outbreak size (log10)",
        y = igraph_var,
      )
}
```

    Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ℹ Please use tidy evaluation idioms with `aes()`.
    ℹ See also `vignette("ggplot2-in-packages")` for more information.

``` r
# Build the same graphs for the ergm_* variables
ergm_vars <- grep("^ergm", names(simstats), value = TRUE)

graphs_ergm <- list()
for (ergm_var in ergm_vars) {
  graphs_ergm[[ergm_var]] <-
    ggplot(simstats) +
      geom_point(aes_string(x = "final_preval", y = ergm_var)) +
      scale_y_log10() +
      labs(
        x = "Outbreak size (log10)",
        y = ergm_var,
      )
}

# Use the patchwork R package to combine the graphs
library(patchwork)
wrap_plots(graphs, ncol = 4)
```

    Warning: Removed 7 rows containing missing values (`geom_point()`).

    Warning: Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).

    Warning: Removed 15 rows containing missing values (`geom_point()`).

    Warning: Removed 7 rows containing missing values (`geom_point()`).

![](03-report_files/figure-commonmark/unnamed-chunk-2-4.png)

``` r
wrap_plots(graphs_ergm, ncol = 4)
```

    Warning: Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).
    Removed 7 rows containing missing values (`geom_point()`).

    Warning: Transformation introduced infinite values in continuous y-axis

    Warning: Removed 7 rows containing missing values (`geom_point()`).

![](03-report_files/figure-commonmark/unnamed-chunk-2-5.png)

## Model fitting

``` r
# Scale all the variables starting with igraph in simstats
igraph_vars <- grep("^igraph", names(simstats), value = TRUE)
simstats[, (igraph_vars) := lapply(.SD, scale), .SDcols = igraph_vars]

# Do the same with the vars starting with ergm
ergm_vars <- grep("^ergm", names(simstats), value = TRUE)
simstats[, (ergm_vars) := lapply(.SD, scale), .SDcols = ergm_vars]
```

``` r
library(texreg)
```

    Version:  1.38.6
    Date:     2022-04-06
    Author:   Philip Leifeld (University of Essex)

    Consider submitting praise using the praise or praise_interactive functions.
    Please cite the JSS article in your publications -- see citation("texreg").

``` r
vars <- c(ergm_vars, igraph_vars)

# Main parameters
vars <- c(vars, "infectiousness", "inc_days", "recovery_rate")

baseline <- paste("~", paste(vars, collapse = " + ")) |>
  as.formula()
```

``` r
# Using the update function, run a series of models
# replacing the lhs of the formula with rt, peak_time, and final_preval

depvars <- c(
  "peak_time", "peak_preval", "rt", "dispersion", "final_preval")

model_peak_time <- update(baseline, paste(depvars[1], "~."))
model_peak_preval <- update(baseline, paste(depvars[2], "~."))
model_rt <- update(baseline, paste(depvars[3], "~."))
model_dispersion <- update(baseline, paste(depvars[4], "~."))
model_final_preval <- update(baseline, paste(depvars[5], "~."))
```

``` r
# Run the models
model_peak_preval <- glm(model_peak_preval, data = simstats)
model_peak_time <- glm(model_peak_time, data = simstats)
model_rt <- glm(model_rt, data = simstats)
model_dispersion <- glm(model_dispersion, data = simstats)
model_final_preval <- glm(model_final_preval, data = simstats, family = poisson())
```

``` r
texreg::knitreg(
  list(
    peak_preval = model_peak_preval,
    peak_time = model_peak_time,
    rt = model_rt,
    dispersion = model_dispersion,
    final_preval = model_final_preval
  ), single.row = TRUE)
```

<table class="texreg" style="margin: 10px auto;border-collapse: collapse;border-spacing: 0px;caption-side: bottom;color: #000000;border-top: 2px solid #000000;">
<caption>
Statistical models
</caption>
<thead>
<tr>
<th style="padding-left: 5px;padding-right: 5px;">
 
</th>
<th style="padding-left: 5px;padding-right: 5px;">
peak_preval
</th>
<th style="padding-left: 5px;padding-right: 5px;">
peak_time
</th>
<th style="padding-left: 5px;padding-right: 5px;">
rt
</th>
<th style="padding-left: 5px;padding-right: 5px;">
dispersion
</th>
<th style="padding-left: 5px;padding-right: 5px;">
final_preval
</th>
</tr>
</thead>
<tbody>
<tr style="border-top: 1px solid #000000;">
<td style="padding-left: 5px;padding-right: 5px;">
(Intercept)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
84.09 (1.14)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
5.31 (0.03)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.59 (0.02)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
1.97 (0.11)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
6.27 (0.00)<sup>\*\*\*</sup>
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_edges
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-68.49 (43.12)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
1.36 (1.28)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.03 (0.70)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-2.29 (4.08)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.10 (0.14)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_nodematch.gender
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.12 (0.29)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.01 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.01 (0.03)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_nodematch.grade
</td>
<td style="padding-left: 5px;padding-right: 5px;">
3.74 (1.91)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.12 (0.06)<sup>\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.04 (0.03)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.24 (0.18)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.01)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_triangle
</td>
<td style="padding-left: 5px;padding-right: 5px;">
3.32 (8.45)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.29 (0.25)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.07 (0.14)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.17 (0.80)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.03 (0.03)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_balance
</td>
<td style="padding-left: 5px;padding-right: 5px;">
63.18 (39.21)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-1.10 (1.16)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.05 (0.64)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
2.48 (3.71)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.08 (0.13)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
ergm_gwdeg.fixed.0.25
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.31 (0.32)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.01 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.04 (0.03)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_modularity
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.28 (0.42)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.01 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.04 (0.04)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_transitivity
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-2.54 (6.44)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.22 (0.19)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.05 (0.11)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.13 (0.61)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.02 (0.02)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_diameter
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.05 (0.20)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.01 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.00)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.02)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_avg_path_length
</td>
<td style="padding-left: 5px;padding-right: 5px;">
31.48 (12.32)<sup>\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.54 (0.37)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.08 (0.20)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.74 (1.16)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.03 (0.04)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_avg_closeness
</td>
<td style="padding-left: 5px;padding-right: 5px;">
33.74 (12.68)<sup>\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.60 (0.38)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.09 (0.21)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.79 (1.20)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.02 (0.04)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
igraph_avg_eigenvector
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.07 (0.21)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.02)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
infectiousness
</td>
<td style="padding-left: 5px;padding-right: 5px;">
388.73 (3.38)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-6.66 (0.10)<sup>\*\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.16 (0.06)<sup>\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.67 (0.32)<sup>\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.01)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
inc_days
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.06 (0.07)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.00 (0.00)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.01 (0.01)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.00 (0.00)
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
recovery_rate
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-1.55 (1.16)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.09 (0.03)<sup>\*\*</sup>
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-0.02 (0.02)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.08 (0.11)
</td>
<td style="padding-left: 5px;padding-right: 5px;">
0.03 (0.00)<sup>\*\*\*</sup>
</td>
</tr>
<tr style="border-top: 1px solid #000000;">
<td style="padding-left: 5px;padding-right: 5px;">
AIC
</td>
<td style="padding-left: 5px;padding-right: 5px;">
39726.94
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4607.36
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-1358.82
</td>
<td style="padding-left: 5px;padding-right: 5px;">
16180.73
</td>
<td style="padding-left: 5px;padding-right: 5px;">
43378.08
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
BIC
</td>
<td style="padding-left: 5px;padding-right: 5px;">
39837.70
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4718.13
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-1248.06
</td>
<td style="padding-left: 5px;padding-right: 5px;">
16291.50
</td>
<td style="padding-left: 5px;padding-right: 5px;">
43482.33
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
Log Likelihood
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-19846.47
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-2286.68
</td>
<td style="padding-left: 5px;padding-right: 5px;">
696.41
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-8073.37
</td>
<td style="padding-left: 5px;padding-right: 5px;">
-21673.04
</td>
</tr>
<tr>
<td style="padding-left: 5px;padding-right: 5px;">
Deviance
</td>
<td style="padding-left: 5px;padding-right: 5px;">
829870.04
</td>
<td style="padding-left: 5px;padding-right: 5px;">
730.59
</td>
<td style="padding-left: 5px;padding-right: 5px;">
221.12
</td>
<td style="padding-left: 5px;padding-right: 5px;">
7422.24
</td>
<td style="padding-left: 5px;padding-right: 5px;">
2838.32
</td>
</tr>
<tr style="border-bottom: 2px solid #000000;">
<td style="padding-left: 5px;padding-right: 5px;">
Num. obs.
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4992
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4992
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4992
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4992
</td>
<td style="padding-left: 5px;padding-right: 5px;">
4992
</td>
</tr>
</tbody>
<tfoot>
<tr>
<td style="font-size: 0.8em;" colspan="6">
<sup>\*\*\*</sup>p \< 0.001; <sup>\*\*</sup>p \< 0.01; <sup>\*</sup>p \<
0.05
</td>
</tr>
</tfoot>
</table>
