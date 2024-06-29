%% This script runs all the tests.
addpath('external','include',"testmats");
warning off
format compact
rng(0);

%% SECTION - test 1 - Runtime comparison of <Taylor of exp> in arbitrary precision arithmetic (AP)

fprintf('Running test of runtime comparison of <Taylor of exp> in arbitrary precision arithmetic...\n');
global_init = tic();
test_exp_taylor_ap_time
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);


%% SECTION - test 2 - Code profiling of <Taylor of exp> in AP

fprintf('Running test on code profiling of <Taylor of exp>...\n');
global_init = tic();
test_profile_table
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION - test 3 - The Cauchy matrix in AP

fprintf('Running test on the Cauchy matrix...\n');
global_init = tic();
test_cauchy_table
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION - test 4 - Accuracy and complexity reduction of <Taylor of exp> in AP

fprintf('Running test on accuracy and complexity reduction of <Taylor of exp> in AP...\n');
global_init = tic();
test_exp_taylor_ap
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);


%% SECTION - test 5 - Accuracy and complexity reduction of <Pade of exp> in double-single-half precision (DSH)

fprintf('Running test on accuracy and complexity reduction of <Pade of exp> in DSH...\n');
global_init = tic();
test_exp_pade_hsd
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION - test 6 - Accuracy and complexity reduction of <Taylor of cos> in AP

fprintf('Running test on accuracy and complexity reduction of <Taylor of cos> in AP...\n');
global_init = tic();
test_cos_taylor_ap
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION - test 7 - Accuracy and complexity reduction of <Pade of exp> in AP

fprintf('Running test on accuracy and complexity reduction of <Pade of exp> in AP...\n');
global_init = tic();
test_exp_pade_ap
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

% exit