%% This script runs all the tests.
warning off
format compact

%% SECTION 3.2.1 - test 1 - Vary m in variable precision arithmetic

fprintf('Running test of varing m in variable precision arithmetic...\n');
global_init = tic();
test_exp_m_mp
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);


%% SECTION 3.2.1 - test 2 - Test the Cauchy matrix

% fprintf('Running test on the Cauchy matrix...\n');
% global_init = tic();
% test_cauchy_mp
% fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);


%% SECTION 3.2.1 - test 3 - Test various matrices in variable precision arithmetic

fprintf('Running test on various matrices in 64 decimal digits...\n');
global_init = tic();
test_exp_mat_mp
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

% fprintf('Running test on various matrices in 256 decimal digits...\n');
% global_init = tic();
% test_exp_mat_mp_u256
% fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION 3.2.2 - test 1 - Test various matrices in double-single-half precision arithmetic

fprintf('Running test on various matrices in double-single-half precision...\n');
global_init = tic();
test_exp_mat_hsd
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION 4.1 - Test various matrices in variable precision arithmetic <Pade Approximants for the matrix exponential>

fprintf('Running test on various matrices in 64 decimal digits <Pade Approximants for the matrix exponential - numerator>...\n');
global_init = tic();
test_exp_pade_numerat_mat_mp
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

fprintf('Running test on various matrices in 64 decimal digits <Pade Approximants for the matrix exponential - denominator>...\n');
global_init = tic();
test_exp_pade_denomin_mat_mp
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

%% SECTION 4.2 - Test various matrices in variable precision arithmetic <Taylor Approximants for the matrix cosine>

fprintf('Running test on various matrices in 64 decimal digits <Taylor Approximants for the matrix cosine>...\n');
global_init = tic();
test_cos_taylor_mat_mp
fprintf('done\t[%5.2f min]\n', toc(global_init) / 60);

% exit