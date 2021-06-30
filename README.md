 
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EpiCompare <a href = "https://skgallagher.github.io/EpiCompare/"><img src="man/figures/logo.png" align="right" alt="EpiCompare" width="120" /></a>

<!-- badges: start -->

[![R-CMD-check, code
coverage](https://github.com/skgallagher/EpiCompare/workflows/R-CMD-check%20and%20code%20coverage/badge.svg)](https://github.com/skgallagher/EpiCompare/actions)
<!-- [![Travis build status](https://travis-ci.org/skgallagher/EpiCompare.svg?branch=master)](https://travis-ci.org/skgallagher/EpiCompare) -->
[![codecov](https://codecov.io/gh/skgallagher/EpiCompare/branch/master/graph/badge.svg)](https://codecov.io/gh/skgallagher/EpiCompare)

<!-- badges: end -->

The goal of `EpiCompare` is to provide the epidemiology community with
easy-to-use tools to encourage comparing and assessing epidemics and
epidemiology models in a time-free manner. All tools attempt to adhere
to `tidyverse`/`ggplot2` style to enhance ease of use.

Additionally, `EpiCompare` provides tools for time invariant analysis.
This allows for ‘fairer’ comparison of epidemics and model based
simulations by avoiding different scaling and shifting of time that can
confound time-based comparisons.

To achieve these goals, the package contains:

  - **Visualization tools** to visualize SIR epidemics and simulations
    from SIR models in a time-free manner using `ggtern`’s ternary plots
    and prediction bands. For agent-based SIR models we also provide
    visualization tools to let the user easily explore how different
    characteristics of the agents relate to different experiences in the
    epidemic.
  - General **comparison tools** to compare epidemics and epidemic
    models that have higher numbers of states (again in a time-free
    manner), allowing for the user to examine the differences between
    models through simulations, and if an epidemic is similar to a model
    through simulations and prediction bands.
  - **Conversion tools** to:
      - Convert and then compare models from *standard epidemic
        packages* like `EpiModels`, `pomp`, as well as internal
        agent-based models, and epidemics in a common framework.
      - Convert *agent-based* information into *aggregate* to compare in
        the aggregate framework described above.

## Installation

You can install the developmental version of EpiCompare from github
using:

``` r
# install.packages("devtools")
devtools::install_github("skgallagher/EpiCompare")
```

## Data

Description of data including in this package can be found in the data
section of the
[reference](https://skgallagher.github.io/EpiCompare/reference/index.html#section-data)
page of the documentation website.

## Example

``` r
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggtern)
library(EpiCompare)
```

The following example comes from a Measles outbreak in Hagelloch,
Germany in 1861. We have data on each child (agent) in the town.

``` r
hagelloch_raw %>%
  ggplot(aes(y = tI, z = tR)) +
    geom_aggregate() + 
    coord_tern() +
    labs(x = "S", y = "I", z = "R", title = "Town Analysis") +
    theme_sir()
#> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Previous work has suggested that the class (`CL`) the student was in
effected how the experienced the outbreak. The below figure shows
differences in the outbreak relative to this grouping.

``` r
hagelloch_raw %>% 
  rename(`school grade` = CL) %>%
  group_by(`school grade`) %>%
  summarize(`number of students` = n())
#> `summarise()` ungrouping output (override with `.groups` argument)
#> # A tibble: 3 x 2
#>   `school grade` `number of students`
#>   <fct>                         <int>
#> 1 preschool                        90
#> 2 1st class                        30
#> 3 2nd class                        68

hagelloch_raw %>%
   ggplot(aes(y = tI, z = tR, color = CL)) +
       geom_aggregate() + 
       coord_tern() +
       labs(x = "S", y = "I", z = "R",
            color = "School Grade",
            title = "Town analysis by grade") +
       theme_sir()
#> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

### Simulate SIR data

``` r
n_sims <- 100
n_time_steps <- 100
beta <- .1
gamma <- .03
init_SIR <- c(950, 50, 0)
    
out <- simulate_SIR_agents(n_sims = n_sims,
                           n_time_steps = n_time_steps,
                           beta = beta, gamma = gamma,
                           init_SIR = init_SIR)
                                 
df_groups <- out %>% dplyr::group_by(sim) %>%
    agents_to_aggregate(states = c(tI, tR)) %>%
    rename(S = X0, I = X1, R = X2)
  
df_groups %>% ggplot() +
    geom_prediction_band(aes(x = S, y = I, z = R, sim_group = as.numeric(sim)),
                         alpha = .2, fill = "blue", color = "blue") +
    geom_line(aes(x = S, y = I, z = R, group = sim), alpha = .1) +
    coord_tern() +
    theme_sir()
#> Warning: Ignoring unknown aesthetics: z
#> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
#> Due to dist_params$dist_approach = "equa_dist", this may take a little while - see `filament_compression` examples for a work-around if you're making this plot multiple times
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## Package Creation Notes:

\*\*We’re transferring to \~github actions\~ and away from Travis CI.
Thanks Travis CI for the long run (During undergrad - probably around
2015, I got introduced to Travis CI and it has been a really great tool
and CIs in general are great tools). Sadly, open source packages (like
ours) no longer gets infinite free resources on Travis. [Dean
Attali](https://deanattali.com/blog/migrating-travis-to-github/) and
[ROpenSci](https://ropensci.org/technotes/2020/11/19/moving-away-travis/)
have blog posts on the situation. As such, if you’re looking to learn
from our mistakes from Travis, then the comments below stop making sense
after December 15th, 2020.

1.  For writing code that works with `tidyverse` 1.0 vs `tidyverse` \<=
    0.8.3. We followed ideas found in [tidyr:
    in-packages](https://tidyr.tidyverse.org/articles/in-packages.html),
    for the code, and - when working with Travis CI (using a matrix for
    multiple builds) - we leverage ideas in [tidyverse travis on
    github](https://github.com/tidyverse/design/issues/95) and
    [tidyverse
    principles](https://principles.tidyverse.org/changes-multivers.html).
    \*\*This is no longer done (was removed 22 December 2020), as the
    rest of `tidyrverse` has moved on and now requires `tidyr >1.0.0`.
2.  For writing your own `geom`s and `stat`s that works with `ggtern`
    (which are generally restricted), the following 2 stack-exchange
    articles helped use do so with ease:
    1.  [stack-exchange: personal code
        idea](https://stackoverflow.com/questions/57175114/ternary-plot-scaling-opacity-across-groups)
    
    2.  [stack-exchange: being able to access ggtern’s element right
        away](https://stackoverflow.com/questions/31589479/how-to-fix-no-item-called-packagepkg-on-the-search-list-without-using-libra#comment51172891_31589479)
    
    3.  Finally, we’ve also leveraged ideas from [R-devel: on avoiding
        problems with
        `:::`](https://stat.ethz.ch/pipermail/r-devel/2013-August/067210.html)
        in `R/aaa.R` to overcome messages from CRAN relative to this
        hack (using `:::`). For some reason - when documenting for
        `pkgdown` website, we need to do `library(ggtern);
        EpiCompare:::update_approved_layers()`
3.  `geom_prediction_band` required not just `compute_group` but
    `compute_layer` - there is very little documentation on how to
    approach this correctly. Basically - there are problems when the
    `compute_group` wants to make multiple `pieces`/`groups` - and it is
    similar to the problem that if you do something like `aes(color =
    var1, group = var2)` you may actually want to do `aes(color = var1,
    group = paste(var1, var2))`, if there are the same `var2` values
    across different `var1` value, but should not necessarily be grouped
    together.
4.  Now that `R` has come out with version \>= 4.0.0, we now need to
    call `.S3method("method", "class")` to define the connection for
    `S3` methods (e.g. `method.class` function), which we have for the
    `contained` function.
5.  **Do two wrongs make a right?** As of 9/23 `ggtern` had an
    [issue](https://bitbucket.org/nicholasehamilton/ggtern/issues/13/ggtern-breaks-ggplots-plot)
    that it messed with `ggplot2`’s legends when loaded (it over-wrote
    the `print.ggplot` and other functions). We’ve over-written
    `ggtern`’s `print.ggplot` to correct this problem when not
    producing ternary plots (code in `aaa.R`).
6.  Useful Rstudio shortcuts for `Roxygen2`: (a) create `Roxygen2`
    comments template with `option` + `command` + `shift` + `R` (b) In
    `Roxygen2` comments do `control` + `shift` + `/` to format relative
    to 80 char limit.
7.  [stack
    overflow](https://stackoverflow.com/questions/11285496/r-cmd-check-warning-functions-methods-with-usage-in-documentation-object-bu)
    post on how to pass `check` for `rownames<-.tidy_dist_mat`.
8.  Transferring from Travis CI to github actions. We only use a single
    workflow file (although we use code ideas found in:
    [check-standard/`usethis::use_github_action_check_standard()`](https://github.com/r-lib/actions/blob/master/examples/check-standard.yaml)),
    [test-coverage/`usethis::use_github_action("test-coverage")`](https://github.com/r-lib/actions/blob/master/examples/test-coverage.yaml).
    In our `.github/workflows/R-CMD-check-coverage.yaml` you’ll find
    **a** (potentially not optional) approach to preform our complex
    checking approach (which tries to copy the ideas in our [old travis
    file](https://github.com/skgallagher/EpiCompare/blob/e1298f799d9248bb330885ba5e7b0fa17ea2c83e/.travis.yml)).
    The github action easily creates a *larger* matrix that allows us to
    run on all standard OS, R versions and our “tidyr-current” vs
    “tidyr-old” split. Our approach with 1 github action `yaml` file
    makes us only look at 8 matrix options - it’s possible that, in the
    future, we’ll go back to 3 files. We currently compile our `pkgdown`
    site on our own computers and then push - and will continue to do so
    due to complications with `reticulate` and `python` packages on the
    virtual machines in github actions.
9.  We’ve speed up some of the interval functions with `Rcpp`. The setup
    to use this functions was slightly complex, we drafted code using
    `cppFunction('...')` and then needed to pipe them over to the
    package setting. Approaches required (1) `usethis::use_rcpp()`, then
    copying C++ code into `src/code.cpp` (needed to make sure we had
    [`// [[Rcpp::export]]`](http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar)
    above each function. For make sure the functions would be compiled
    and available we needed to update the package documentation with [2
    roxygen tags](https://r-pkgs.org/src.html#cpp):

<!-- end list -->

    #' @useDynLib your-package-name
    #' @importFrom Rcpp sourceCpp

## Contributors

  - Shannon Gallagher ([`skgallagher`](https://github.com/skgallagher))
  - Benjamin LeRoy ([`benjaminleroy`](https://github.com/benjaminleroy))
