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
M_.fname = 'rbcdata';
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
M_.endo_names = cell(9,1);
M_.endo_names_tex = cell(9,1);
M_.endo_names_long = cell(9,1);
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
M_.endo_names(6) = {'ly_l'};
M_.endo_names_tex(6) = {'ly\_l'};
M_.endo_names_long(6) = {'ly_l'};
M_.endo_names(7) = {'lw'};
M_.endo_names_tex(7) = {'lw'};
M_.endo_names_long(7) = {'lw'};
M_.endo_names(8) = {'Rk'};
M_.endo_names_tex(8) = {'Rk'};
M_.endo_names_long(8) = {'Rk'};
M_.endo_names(9) = {'z'};
M_.endo_names_tex(9) = {'z'};
M_.endo_names_long(9) = {'z'};
M_.endo_partitions = struct();
M_.param_names = cell(5,1);
M_.param_names_tex = cell(5,1);
M_.param_names_long = cell(5,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'chi'};
M_.param_names_tex(2) = {'chi'};
M_.param_names_long(2) = {'chi'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'rho'};
M_.param_names_tex(5) = {'rho'};
M_.param_names_long(5) = {'rho'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 9;
M_.param_nbr = 5;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
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
M_.orig_eq_nbr = 9;
M_.eq_nbr = 9;
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
 0 4 12;
 1 5 0;
 0 6 0;
 0 7 13;
 0 8 0;
 0 9 0;
 0 10 0;
 2 11 14;]';
M_.nstatic = 5;
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
  8 , 'name' , '8' ;
  9 , 'name' , 'z' ;
};
M_.mapping.ly.eqidx = [4 5 7 8 ];
M_.mapping.lc.eqidx = [1 2 4 ];
M_.mapping.lk.eqidx = [1 3 5 6 8 ];
M_.mapping.li.eqidx = [4 6 ];
M_.mapping.lh.eqidx = [1 2 3 5 7 ];
M_.mapping.ly_l.eqidx = [7 ];
M_.mapping.lw.eqidx = [2 3 ];
M_.mapping.Rk.eqidx = [8 ];
M_.mapping.z.eqidx = [1 3 5 9 ];
M_.mapping.e.eqidx = [9 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 9 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(5, 1);
M_.endo_trends = struct('deflator', cell(9, 1), 'log_deflator', cell(9, 1), 'growth_factor', cell(9, 1), 'log_growth_factor', cell(9, 1));
M_.NNZDerivatives = [31; -1; -1; ];
M_.static_tmp_nbr = [6; 2; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
close all;
[param,ss] = calibration;
M_.params(4) = param(1);
alpha = M_.params(4);
M_.params(1) = param(2);
beta = M_.params(1);
M_.params(3) = param(3);
delta = M_.params(3);
M_.params(2) = param(4);
chi = M_.params(2);
M_.params(5) = 0.95;
rho = M_.params(5);
sigma   = 0.01;
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(3) = log(ss(1));
oo_.steady_state(2) = log(ss(2));
oo_.steady_state(5) = log(ss(3));
oo_.steady_state(4) = log(ss(4));
oo_.steady_state(1) = log(ss(5));
oo_.steady_state(7) = log(ss(6));
oo_.steady_state(8) = log(ss(7));
oo_.steady_state(6) = oo_.steady_state(1)-oo_.steady_state(5);
oo_.steady_state(9) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = sigma^2;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.order = 1;
options_.periods = 1000;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
datatomfile('simuldataRBC',[]);
return;


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'rbcdata_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
