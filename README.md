# Scalable Bayesian $p$-generalized probit and logistic regression

## Abstract
This R package implements a novel Bayesian binary regression framework using the $p$-generalized Gaussian distribution ($p$-GGD) to enhance model flexibility, particularly in tail behavior, thus offering improved fit over traditional logistic and probit models. The method employs Markov Chain Monte Carlo (MCMC) sampling for the efficient estimation of model parameters $\beta$ and the link function parameter $p$. We demonstrate the utility of this approach through both simulated and real-world datasets, and address large data challenges by integrating coreset-based data reduction techniques, ensuring efficient computational performance with minimal compromise on accuracy.

## Getting Started

### Prerequisites
Ensure you have R installed on your system (version 4.0.0 or newer recommended). The following R packages are required:
- `dplyr`: For data manipulation.
- `ggplot2`: For data visualization.

### Installation
You can install the package directly from GitHub using `devtools`:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("zeyudsai/BayesPprobit")
```

## Usage
Load the library and follow the vignettes for detailed examples on how to apply the Bayesian binary regression model to your data:
```r
library(BayesPprobit)

# For a detailed guide
vignette("bayesian_binary_regression_guide")
```

## Core Features and Parameters
- **$p$-Generalized Gaussian Distribution ($p$-GGD)**: Offers a flexible approach to modeling distribution tails through the parameter $p$.
- **Bayesian Estimation with MCMC Sampling**: Implements Bayesian posterior estimation using MCMC, detailed for both $\beta$ and $p$.
- **Coresets for Large Data**: Utilizes coreset-based techniques to facilitate efficient computations on large datasets with minimal accuracy loss.

## Contributing
We welcome contributions to improve this package. To contribute:
1. Fork the repo.
2. Create your feature branch (`git checkout -b feature/AmazingFeature`).
3. Commit your changes (`git commit -am 'Add some AmazingFeature'`).
4. Push to the branch (`git push origin feature/AmazingFeature`).
5. Create a new Pull Request.

## Citation
If this package aids in your research, kindly cite it as follows:
```
@article{xx,
  title={Scalable Bayesian $p$-generalized probit and logistic regression},
  author={Zeyu Ding and Simon Omlor and Katja Ickstadt and Alexander Munteanu},
  journal={Advances },
  year={2024},
  volume={xx},
  number={xx},
  pages={xx-xx},
  publisher={Publisher}
}
```

## License
This package is released under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
Zeyu Ding - zeyu.ding@tu-dortmund.de


