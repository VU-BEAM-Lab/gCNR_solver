function [gCNR, varargout] = gCNR_solver(target,reference,varargin)
% Universal gCNR solver, developed with the recommendations from the work by Schlunk and Byram [1] in mind.
% gCNR is a robust CNR-like method, developed by Rodriguez-Molares et al. [2], [3].
% It produces a value between 0 and 1, that represents lesion detectability
%
% By default, gCNR_solver(target,reference) will produce a gCNR estimate using histograms with k=ceil(2*N^(2/5)) variable-width bins,
% a technique recommended by others [4] and generally in literature [5],[6], that we found to be robust in our testing [1].
%
% For other techniques, or more control over the estimation parameters, we list here example calls as well as additional inputs and outputs.
%
% gCNR = gCNR_solver(target,reference) returns the gCNR estimated using histograms with k=ceil(2*N^(2/5)) variable-width bins
% gCNR = gCNR_solver(target,reference,method) specify the specific method used for estimating gCNR. For example, 'uniform', 'variable', 'ecdf'
% gCNR = gCNR_solver(target,reference,binning) specify the specific binning method used if applicable. For example, 'sqrt' or 'cuberoot'
% [gCNR, gCNR_lower, gCNR_upper] = gCNR_solver(target,reference,'ecdf') when specifying 'ecdf' for the estimation method, the lower and upper confidence bounds of the gCNR estimate can be returned
% gCNR = gCNR_Solver(target,reference,Name,Value) specify certain name-value pair arguments. For example, 'k',100 would specifies to use a fixed number of bins (k=100)
%
% INPUTS 
% target, reference - vectors or matrices of the target and reference regions being compared
% varargin          - additional arguments listed below. These are divided between string flags and name-value pairs
%
% Flags - i.e., just including the string will flag the code for the specific configuration
% 'uniform'     - use histograms with uniform bin widths (will default to 'sqrt' binning if not otherwise specified)
% 'variable'    - (DEFAULT) use histograms with variable bin widths (will default to 'equiprobable' binning if not otherwise specified)
% 'ecdf'        - use eCDFs instead of histograms (does not require any additional input parameters)
% 'rank'        - specifies whether to rank-order the input data (this will estimate gCNR on the lists of the rankings rather than the raw data)
% 'sqrt'        - choose ceil(sqrt(N)) bins based on the total amount of data length(target). Note that ideally the size of the target and reference regions should be approximately equal
% 'cuberoot'    - choose ceil(N^(1/3)) bins
% 'equiprobable'- choose ceil(2*N^(2/5)) bins (recommended for variable bin widths [5],[6])
%
% name-value argument pairs
% 'k'           - whole number, specifies an exact number of bins to use
%
% Note that the 'ecdf' method can produce different results compared to the histogram-based methods due to the difference in the estimation method.
% The amount of variance can scale inversely based on the amount of data (less data -> more variance). The histogram-based estimate will still
% fall within the expected error of the 'ecdf' result.
%
% Citations
% [1] Schlunk, S., & Byram, B. (2023). Methods for enhancing the robustness of the generalized contrast-to-noise ratio. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 1–1. https://doi.org/10.1109/TUFFC.2023.3289157
% [2] Rodriguez-Molares, A., Rindal, O. M. H., D’Hooge, J., Måsøy, S.-E., Austeng, A., & Torp, H. (2018). The Generalized Contrast-to-Noise Ratio. IEEE International Ultrasonics Symposium (IUS), 1–4.
% [3] Rodriguez-Molares, A., Rindal, O. M. H., D’hooge, J., Masoy, S.-E., Austeng, A., Lediju Bell, M. A., & Torp, H. (2020). The Generalized Contrast-to-Noise Ratio: A Formal Definition for Lesion Detectability. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 67(4), 745–759. https://doi.org/10.1109/TUFFC.2019.2956855
% [4] Hyun, D., Kim, G. B., Bottenus, N., & Dahl, J. J. (2022). Ultrasound Lesion Detectability as a Distance Between Probability Measures. IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 69(2), 732–743. https://doi.org/10.1109/TUFFC.2021.3138058
% [5] https://en.wikipedia.org/wiki/Histogram#Variable_bin_widths
% [6] J. Prins, D. McCormack, D. Michelson, and K. Horrell, "Chi-square goodness-of-fit test,” p. 7.2.1.1, 2003. [Online]. Available: https://itl.nist.gov/div898/handbook/prc/section2/prc211.htm
%
% Change Log
% 2023 10 04    Siegfried Schlunk
%               Original Version

%% Defaults
method = 'variable';
binning = 'equiprobable';
k = ceil(2*length(target(:))^(2/5));
used_rank_ordered_data = false;
num_additional_outputs = 0;

in_idx = 1;
while in_idx <= nargin - 2
    switch varargin{in_idx}
        case 'uniform'
            method = 'uniform';
            binning = 'sqrt';
            k = ceil(sqrt(length(target(:))));
        case 'variable'
            method = 'variable';
            binning = 'equiprobable';
            k = ceil(2*length(target(:))^(2/5));
        case 'ecdf'
            method = 'ecdf';
            binning = 'N/a';
            k = [];
        case 'rank'
            used_rank_ordered_data = true;
            [~,rank_idx] = sort([target(:); reference(:)]);
            target = find(rank_idx<=length(target(:)));
            reference = find(rank_idx>length(target(:)));
        case 'sqrt'
            binning = 'sqrt';
            k = ceil(sqrt(length(target(:))));
        case 'cuberoot'
            binning = 'cuberoot';
            k = ceil(length(target(:))^(1/3));
        case 'equiprobable'
            binning = 'equiprobable';
            k = ceil(2*length(target(:))^(2/5));
        case 'k'
            binning = 'manual';
            in_idx = in_idx + 1;
            k = varargin{in_idx};
    end
    in_idx = in_idx + 1;
end
params.method = method;
params.binning = binning;
params.k = k;
params.used_rank_ordered_data = used_rank_ordered_data;

switch method
    case 'uniform'
        [gCNR] = gCNR_histogram_solver(target,reference,k,method);
        params.gCNR = gCNR;
    case 'variable'
        [gCNR] = gCNR_histogram_solver(target,reference,k,method);
        params.gCNR = gCNR;
    case 'ecdf'
        [gCNR,gCNR_lower,gCNR_upper] = gCNR_eCDF_solver(target,reference);
        varargout{num_additional_outputs+1} = gCNR_lower;
        varargout{num_additional_outputs+2} = gCNR_upper;
        num_additional_outputs = num_additional_outputs + 2;
        params.gCNR = gCNR;
        params.gCNR_lower = gCNR_lower;
        params.gCNR_upper = gCNR_upper;
end

varargout{num_additional_outputs+1} = params;

end

function [gCNR] = gCNR_histogram_solver(target,reference,k,varargin)
lower_bound = min(min(target(:)),min(reference(:)));
upper_bound = max(max(target(:)),max(reference(:)));
Xbin = linspace(lower_bound,upper_bound,k+1);
if nargin > 3
    if strcmp(varargin{1},'variable')
        avg_counts = floor((length(target(:))+length(reference(:)))/k);
        total = sort([target(:); reference(:)]);
        Xbin = total((1:k-1)*avg_counts);
        Xbin = [lower_bound; Xbin; upper_bound];
    end
end
counts1 = histcounts(target,Xbin);
counts2 = histcounts(reference,Xbin);
p1 = counts1/sum(counts1);
p2 = counts2/sum(counts2);
minCount = min(p1,p2);
OVL = sum(minCount);
gCNR = 1-OVL;
end

function [gCNR,lower,upper] = gCNR_eCDF_solver(target,reference)
% calculate eCDF for target an reference data
% significance level 'Alpha' is default set to 0.05, this will control the
% confidence interval (f_lower,f_upper)
[ft,xt,ft_lower,ft_upper] = ecdf(target,'Alpha',0.05);
[fr,xr,fr_lower,fr_upper] = ecdf(reference,'Alpha',0.05);
% weird ecdf implementation behavior, fixing up bounds
xt(1) = xt(2)-(xt(3)-xt(2));
ft_lower(1) = 0;
ft_lower(end) = 1;
ft_upper(1) = 0;
ft_upper(end) = 1;
xr(1) = xr(2)-(xr(3)-xr(2));
fr_lower(1) = 0;
fr_lower(end) = 1;
fr_upper(1) = 0;
fr_upper(end) = 1;
% interpolate for consistent sampling between ft, fr. 'nearest' is
% specified to ensure extrapolated points are set to 0 or 1
x = unique(sort([xt(:); xr(:)]));
ft = interp1(xt,ft,x,'nearest','extrap');
ft_lower = interp1(xt,ft_lower,x,'nearest','extrap');
ft_upper = interp1(xt,ft_upper,x,'nearest','extrap');
fr = interp1(xr,fr,x,'nearest','extrap');
fr_lower = interp1(xr,fr_lower,x,'nearest','extrap');
fr_upper = interp1(xr,fr_upper,x,'nearest','extrap');
% estimate gCNR, lower, and upper bounds
% calculate max error between the upper and lower bounds, which we use to
% automatically choose minimum prominence for local extrema
max_error = max(max(ft_upper-ft_lower),max(fr_upper-fr_lower))/2;

% ecdf automatically rounds to ceil(log10(length(target(:))) decimals
max_decimals_to_round = ceil(log10(max(length(target(:)),length(reference(:)))));

% find difference of the target and reference eCDFs, round to value above
diff_fun = round(ft-fr,max_decimals_to_round);

% MATLAB peak finding functions can have bugs when we have two identical
% peak values in close proximity (e.g. 1-2-2-1 or 1-2-1-2-1) where multiple
% peaks of the same value should be returned rather than just the one. This
% can result in gCNR that is double the expected value if not corrected.

% We correct the peak issue by adding a small jitter: a monotonically
% increasing function to diff_fun that guarantees all the values of diff_fun
% are unique. Specifically, we add a linear function that is on the order of
% exp(-(max_decimals_to_round+1)) so that that added values are too small
% to impact the actual peak locations.
diff_fun = diff_fun + linspace(1,9,length(diff_fun))'*10^-(max_decimals_to_round+1);

% 'MinPominence' is based on the max_error calculated before. max_error
% generally scales with the size of the dataset (smaller datasets have more
% error in the eCDF estimate). This helps to prevent some issues that extra
% small datasets can cause.
maxlocal = islocalmax(diff_fun,'MinProminence',max_error);
minlocal = islocalmin(diff_fun,'MinProminence',max_error);

% if 'MinProminence' is too large due to a small dataset, will manually
% return only the global min and max of the data.
if ~sum(maxlocal)
    [~,idx] = max(diff_fun);
    maxlocal(idx) = 1;
end
if ~sum(minlocal)
    [~,idx] = min(diff_fun);
    minlocal(idx) = 1;
end

% re-calculate diff_fun without rounding or jitter
diff_fun = ft-fr;

% by our proof, gCNR is the difference of the sum of the max and min values
gCNR = sum(diff_fun(maxlocal)) - sum(diff_fun(minlocal));

% find the lower and upper estimate of gCNR based on the confidence
% interval initially estimated
diff_fun = ft_lower-fr_upper;
temp1 = sum(diff_fun(maxlocal)) - sum(diff_fun(minlocal));
diff_fun = ft_upper-fr_lower;
temp2 = sum(diff_fun(maxlocal)) - sum(diff_fun(minlocal));
lower = min([temp1 temp2]);
upper = max([temp1 temp2]);
end