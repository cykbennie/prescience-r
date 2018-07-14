# prescience-r

R package of the Approximate Best Subset Maximum Binary Prediction Rule (PRESCIENCE) proposed by Chen and Lee (2018).

Description of this prediction method and its computation details can be found in the paper:
[Chen, Le-Yu and Lee, Sokbae (2018), "Best Subset Binary Prediction"](https://arxiv.org/abs/1610.02738).

The Matlab implementation of PRESCIENCE can be found at https://github.com/LeyuChen/Best-Subset-Binary-Prediction.

## Requirements

The `prescience` package requires the following three R packages:

* `gurobi`
* `slam`
* `stats`

**Important**: You must download the [Gurobi Optimizer](http://www.gurobi.com/) software before being able to install the `gurobi` package. The Gurobi Optimizer is free under academic liscence. See their [official documentation](http://www.gurobi.com/documentation/) for the software installation guide.

``` r
install.packages("slam")
install.packages("stats")
# For the installation of the gurobi package, see http://www.gurobi.com/documentation/.
```

## Installation

The package can be installed from github:

``` r
devtools::install_github("cykbennie/prescience-r")
```

## Reference Manual

Refer to the [prescience.pdf](prescience.pdf) file for details.

## Example

``` r
# Load the package
library(prescience)

# Create the PRESCIENCE object
results <- select(auto ~ dcost + cars + dovtt + divtt, data = transportation, nfoc = 1, q = 1, bound = 10, beta0 = 1, warmstart = TRUE, tau = 1.5, mio = 1, tlim = 86400)

# Summary of PRESCIENCE
summary(results)

# Estimated coefficients of PRESCIENCE
coef(results)
```

## Author

* Yankang (Bennie) Chen (<yankang.chen@columbia.edu>)

## License

This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Le-Yu Chen (<lychen@econ.sinica.edu.tw>)
* Sokbae Lee (<sl3841@columbia.edu>)

The codes in this R package are developed based on the Matlab implementation by Le-Yu Chen (https://github.com/LeyuChen/Best-Subset-Binary-Prediction).
