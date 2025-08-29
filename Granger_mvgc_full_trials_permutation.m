%% Start ToolBox
startup_mvgc1 % path /Users/flavio/Documents/MATLAB/MVGC1-master
clc

clear('have_genvar_mex','mvgc_root','mvgc_version')

%%

acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
mvgc.parameters.regmode   = 'LWR';                                                            % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
mvgc.parameters.acmaxlags = [];                                                               % maximum autocovariance lags (empty for automatic calculation)
mvgc.parameters.acdectol  = [];

%    fP = permtest_tsdata_to_spwcgc(U,p,fres,bsize,nsamps,regmode,acmaxlags,acdectol)
%
% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     p          model order (number of lags)
%     fres       frequency resolution (default: automatic)
%     bsize      permutation block size (default: use model order)
%     nsamps     number of permutations
%     regmode    regression mode (default as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags  (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     fP         permutation test spectral Granger causalities (null distribution)
%


% Baseline
Granger_mvgc_permutation{1,ms} = permtest_tsdata_to_spwcgc(mvgc.data{1,1},mvgc.parameters.morder,mvgc.parameters.fres,mvgc.parameters.morder,1000,mvgc.parameters.regmode,mvgc.parameters.acmaxlags,mvgc.parameters.acdectol);

% CS_TONE
Granger_mvgc_permutation{2,ms} = permtest_tsdata_to_spwcgc(mvgc.data{2,1}(:,:,CSIT{ms)},mvgc.parameters.morder,mvgc.parameters.fres,mvgc.parameters.morder,1000,mvgc.parameters.regmode,mvgc.parameters.acmaxlags,mvgc.parameters.acdectol);

% ITI
Granger_mvgc_permutation{3,ms} = permtest_tsdata_to_spwcgc(mvgc.data{3,1}(:,:,CSIT_1{ms)},mvgc.parameters.morder,mvgc.parameters.fres,mvgc.parameters.morder,1000,mvgc.parameters.regmode,mvgc.parameters.acmaxlags,mvgc.parameters.acdectol);



