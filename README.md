# prescience-r

R package of the Approximate Best Subset Maximum Binary Prediction Rule (PRESCIENCE) proposed by Chen and Lee (2018).

Description of this prediction method and its computation details can be found in the paper:
Chen, Le-Yu and Lee, Sokbae (2018), "Best Subset Binary Prediction".

The above paper and the Matlab implementation of PRESCIENCE can be found at https://github.com/LeyuChen/Best-Subset-Binary-Prediction.

## Requirements

This R package requires the following R packages:

* [gurobi](http://www.gurobi.com/documentation/) - The Gurobi Optimizer is free under academic liscence
* [slam]
* [stats]

```
install.packages("slam")
install.packages("stats")
# See http://www.gurobi.com/documentation/ for the installation of the gurobi package
```

## Installation

The package can be installed from github:

```
devtools::install_github("cykbennie/prescience-r")
```

## Example

```
library(prescience)
results <- select(auto ~ dcost + cars + dovtt + divtt, data = transportation, nfoc = 1, q = 1, bound = 10)
summary(results)
coef(results)
```

## Author

* Yankang (Bennie) Chen <yankang.chen@columbia.edu>

## License

This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Le-Yu Chen <lychen@econ.sinica.edu.tw>
* Sokbae Lee <sl3841@columbia.edu>

The codes in this R package are developed based on the Matlab implementation by Le-Yu Chen (https://github.com/LeyuChen/Best-Subset-Binary-Prediction).
