## gCNR_solver
<summary>Universal gCNR solver, developed with the recommendations from the work by Schlunk and Byram [1] in mind. gCNR is a robust CNR-like method, developed by Rodriguez-Molares et al. [2], [3]. It produces a value between 0 and 1, that represents lesion detectability.</summary>

By default, gCNR_solver(target,reference) will produce a gCNR estimate using histograms with k=ceil(2*N^(2/5)) variable-width bins, a technique recommended by others [4] and generally in literature [5],[6], that we found to be robust in our testing [1].

For other techniques, or more control over the estimation parameters, we list here example calls as well as additional inputs and outputs.

| Example function call | Description |
|----------------------:|-------------|
|gCNR = gCNR_solver(target,reference) | returns the gCNR estimated using histograms with k=ceil(2*N^(2/5)) variable-width bins|
|gCNR = gCNR_solver(target,reference,method) | specify the specific method used for estimating gCNR. For example, 'uniform', 'variable', 'ecdf'|
|gCNR = gCNR_solver(target,reference,binning) | specify the specific binning method used if applicable. For example, 'sqrt' or 'cuberoot'|
|[gCNR,gCNR_lower,gCNR_upper] = gCNR_solver(target,reference,'ecdf') | when specifying 'ecdf' for the estimation method, the lower and upper confidence bounds of the gCNR estimate can be returned|
|gCNR = gCNR_Solver(target,reference,Name,Value) | specify certain name-value pair arguments. For example, 'k',100 would specifies to use a fixed number of bins (k=100)|

### INPUTS 
target, reference - vectors or matrices of the target and reference regions being compared

varargin          - additional arguments listed below. These are divided between string flags and name-value pairs

|Flags         | just including the string will flag the code for the specific configuration|
|-------------:|---------------------------------------------------------------------------------------------------|
|'uniform'     | use histograms with uniform bin widths (will default to 'sqrt' binning if not otherwise specified)|
|'variable'    | (DEFAULT) use histograms with variable bin widths (will default to 'equiprobable' binning if not otherwise specified)|
|'ecdf'        | use eCDFs instead of histograms (does not require any additional input parameters)|
|'rank'        | specifies whether to rank-order the input data (this will estimate gCNR on the lists of the rankings rather than the raw data)|
|'sqrt'        | choose ceil(sqrt(N)) bins based on the total amount of data length(target). Note that ideally the size of the target and reference regions should be approximately equal|
|'cuberoot'    | choose ceil(N^(1/3)) bins|
|'equiprobable'| choose ceil(2*N^(2/5)) bins (recommended for variable bin widths [5],[6])|

| name | value |
|-----:|-------|
| 'k' | whole number, specifies an exact number of bins to use|

Note that the 'ecdf' method can produce different results compared to the histogram-based methods due to the difference in the estimation method.
The amount of variance can scale inversely based on the amount of data (less data -> more variance). The histogram-based estimate will still
fall within the expected error of the 'ecdf' result.

| Citations | |
|----------:|-|
|[1]| Schlunk, S., & Byram, B. (2023). Methods for enhancing the robustness of the generalized contrast-to-noise ratio. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 1–1. https://doi.org/10.1109/TUFFC.2023.3289157|
|[2]| Rodriguez-Molares, A., Rindal, O. M. H., D’Hooge, J., Måsøy, S.-E., Austeng, A., & Torp, H. (2018). The Generalized Contrast-to-Noise Ratio. IEEE International Ultrasonics Symposium (IUS), 1–4.|
|[3]| Rodriguez-Molares, A., Rindal, O. M. H., D’hooge, J., Masoy, S.-E., Austeng, A., Lediju Bell, M. A., & Torp, H. (2020). The Generalized Contrast-to-Noise Ratio: A Formal Definition for Lesion Detectability. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 67(4), 745–759. https://doi.org/10.1109/TUFFC.2019.2956855|
|[4]| Hyun, D., Kim, G. B., Bottenus, N., & Dahl, J. J. (2022). Ultrasound Lesion Detectability as a Distance Between Probability Measures. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 69(2), 732–743. https://doi.org/10.1109/TUFFC.2021.3138058|
|[5]| https://en.wikipedia.org/wiki/Histogram#Variable_bin_widths|
|[6]| J. Prins, D. McCormack, D. Michelson, and K. Horrell, "Chi-square goodness-of-fit test,” p. 7.2.1.1, 2003. [Online]. Available: https://itl.nist.gov/div898/handbook/prc/section2/prc211.htm|
