%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'rbcEST2';
M_.dynare_version = '5.5';
oo_.dynare_version = '5.5';
options_.dynare_version = '5.5';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e'};
M_.exo_names_tex(1) = {'e'};
M_.exo_names_long(1) = {'e'};
M_.endo_names = cell(8,1);
M_.endo_names_tex = cell(8,1);
M_.endo_names_long = cell(8,1);
M_.endo_names(1) = {'ly'};
M_.endo_names_tex(1) = {'ly'};
M_.endo_names_long(1) = {'ly'};
M_.endo_names(2) = {'lc'};
M_.endo_names_tex(2) = {'lc'};
M_.endo_names_long(2) = {'lc'};
M_.endo_names(3) = {'lk'};
M_.endo_names_tex(3) = {'lk'};
M_.endo_names_long(3) = {'lk'};
M_.endo_names(4) = {'li'};
M_.endo_names_tex(4) = {'li'};
M_.endo_names_long(4) = {'li'};
M_.endo_names(5) = {'lh'};
M_.endo_names_tex(5) = {'lh'};
M_.endo_names_long(5) = {'lh'};
M_.endo_names(6) = {'lw'};
M_.endo_names_tex(6) = {'lw'};
M_.endo_names_long(6) = {'lw'};
M_.endo_names(7) = {'Rk'};
M_.endo_names_tex(7) = {'Rk'};
M_.endo_names_long(7) = {'Rk'};
M_.endo_names(8) = {'z'};
M_.endo_names_tex(8) = {'z'};
M_.endo_names_long(8) = {'z'};
M_.endo_partitions = struct();
M_.param_names = cell(5,1);
M_.param_names_tex = cell(5,1);
M_.param_names_long = cell(5,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'delta'};
M_.param_names_tex(2) = {'delta'};
M_.param_names_long(2) = {'delta'};
M_.param_names(3) = {'chi'};
M_.param_names_tex(3) = {'chi'};
M_.param_names_long(3) = {'chi'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'rho'};
M_.param_names_tex(5) = {'rho'};
M_.param_names_long(5) = {'rho'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 8;
M_.param_nbr = 5;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
options_.varobs = cell(1, 1);
options_.varobs(1)  = {'ly'};
options_.varobs_id = [ 1  ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 3 0;
 0 4 11;
 1 5 0;
 0 6 0;
 0 7 12;
 0 8 0;
 0 9 0;
 2 10 13;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [7; 1; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , '3' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , '6' ;
  7 , 'name' , '7' ;
  8 , 'name' , 'z' ;
};
M_.mapping.ly.eqidx = [4 5 7 ];
M_.mapping.lc.eqidx = [1 2 4 ];
M_.mapping.lk.eqidx = [1 3 5 6 7 ];
M_.mapping.li.eqidx = [4 6 ];
M_.mapping.lh.eqidx = [1 2 3 5 ];
M_.mapping.lw.eqidx = [2 3 ];
M_.mapping.Rk.eqidx = [7 ];
M_.mapping.z.eqidx = [1 3 5 8 ];
M_.mapping.e.eqidx = [8 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 8 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(5, 1);
M_.endo_trends = struct('deflator', cell(8, 1), 'log_deflator', cell(8, 1), 'growth_factor', cell(8, 1), 'log_growth_factor', cell(8, 1));
M_.NNZDerivatives = [28; -1; -1; ];
M_.static_tmp_nbr = [6; 2; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
close all;
clc;
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 4, NaN, (-Inf), Inf, 1, 0.35, 0.02, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 1, NaN, (-Inf), Inf, 1, 0.99, 0.002, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, NaN, (-Inf), Inf, 1, 0.025, 0.003, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, NaN, (-Inf), Inf, 1, 0.9, 0.05, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, NaN, (-Inf), Inf, 4, 0.01, Inf, NaN, NaN, NaN ];
options_.mh_drop = 0.45;
options_.mh_jscale = 0.8;
options_.mh_nblck = 2;
options_.mh_replic = 2000;
options_.mode_compute = 6;
options_.order = 1;
options_.datafile = 'simuldataRBC';
options_.first_obs = 500;
options_.nobs = 400;
var_list_ = {'ly'};
oo_recursive_=dynare_estimation(var_list_);
options_.order = 1;
options_.periods = 1000;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcEST2_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
