### CBModelTools.jl: A Julia package for working with COBRA model files
The [COBRA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3319681/) model format is a common
[MATLAB](https://www.mathworks.com/products/matlab.html)-based binary storage format used by constraint based metabolic modeling tools, including the [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) in [MATLAB](https://www.mathworks.com/products/matlab.html).
This [Julia](http://julialang.org) package contains a number of utility methods to manipulate/generate
[COBRA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3319681/) model files from sources outside
of [MATLAB](https://www.mathworks.com/products/matlab.html).

### Installation and requirements
``CBModelTools.jl`` can be installed in the ``package mode`` of Julia.

Start of the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/index.html) and enter the ``package mode`` using the ``]`` key (to get back press the ``backspace`` or ``^C`` keys). Then, at the prompt enter:

    (v1.1) pkg> add https://github.com/varnerlab/CBModelTools.git

This will install the `CBModelTools.jl` package and other all required packages.
``CBModelTools.jl`` requires Julia 1.x and above.

### Funding
The work described was supported by the [Center on the Physics of Cancer Metabolism at Cornell University](https://psoc.engineering.cornell.edu) through Award Number 1U54CA210184-01 from the [National Cancer Institute](https://www.cancer.gov). The content is solely the responsibility of the authors and does not necessarily
represent the official views of the [National Cancer Institute](https://www.cancer.gov) or the [National Institutes of Health](https://www.nih.gov).
