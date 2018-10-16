%% Example Script to validate Anduril
% Last Update:  8-June-2018
% Authors:  Georgios Leontaris and Oswaldo Morales-Nápoles
% email:    G.Leontaris@tudelft.nl & O.MoralesNapoles@tudelft.nl

% <ANDURIL: A toolbox for structured expert judgment>
%     Copyright (C) <2017>  <Georgios Leontaris & Oswaldo Morales-Nápoles>
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     The only official release of ANDURIL is by the authors. If you wish 
%     to contribute to the official release of ANDURIL please contact the 
%     authors. The authors will decide which contributions would enter the 
%     official release of ANDURIL. 
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>
%--------------------------------------------------------------------------

% clear all; 
close all; clc;

%% Load and formulate data
load('assessments_EC_Grow_Mex.mat')

for j = 1:19
    back_measure{j} = 'uni';
%     back_measure{j} = 'log_uni';
end
% back_measure{3} = 'uni'; % uncomment this line if log-uni background
                           % measure is used for this example
c = 0;
z = 0;
for i = 1:9 % expert
    c = c+1;
    for j = 1:3 % quantile
        for l = 1:19 % item
            z = (c-1)*19 + l;
            Cal_var(i,j,l) = Q(j,z);
        end
    end
end
% load realization
load('realizations_EC_Mex.mat')
Cal_var_new = Cal_var;
N_cal_var = 13;
f = 0;
for i = 1:9 % # of experts
    for j = 1:3 % # of quantiles
        for l = 1:19 % # of all items
            if isempty(realization{l})
                TQs(i,j,l-N_cal_var) = Cal_var(i,j,l);
                f = f + 1;
                excl(f) = l;
            end
        end
    end
end

excl_item = unique(excl);
Cal_var_new(:,:,excl_item) = [];

%% Calculate Decision makers with different weighting schemes
% Parameters
alpha = 0.05; % significance level
k = 0.1; % 
global cal_power
cal_power = 1; % this value should be between [0.1, 1]

% Calculation of DM using global weights
W = global_weights(Cal_var_new, TQs, realization, alpha, back_measure, k);
% a high value of significance level alpha could lead to zero weights for 
% every expert. In this case, a value lower than the highest calibration 
% score should be assigned
if isequal(W(:,4),zeros(size(W,1),1))
   error('Significance Level value should be smaller than the highest calibration score: %d',max(W(:,1))) 
end
[f_DM1, F_DM_out1, X_out1, DM1, W_incl_DM1] = calculate_DM_global(Cal_var_new,...
    TQs, realization, W(:,5)', k, back_measure, alpha);

% % Calculation of DM using item weights 
[unorm_w, W_itm, W_itm_tq] = item_weights(Cal_var_new, TQs, realization, alpha, back_measure, k);
% a high value of significance level alpha could lead to zero weights for 
% every expert for every item. In this case, a value lower than the highest
% calibration score should be assigned 
if isequal(W_itm(:,1),zeros(size(W_itm,1),1))
   error('Significance Level value should be smaller than the highest calibration score') 
end
[f_DM2, F_DM_out2, X_out2, DM2, W_incl_DM2] = calculate_DM_item(Cal_var_new,...
    TQs, realization, W_itm, W_itm_tq, k, back_measure, alpha);

% Calculation of DM using equal weights
for i = 1:size(Cal_var,1)
    eq_w(i) = 1/size(Cal_var,1);
end
alpha_eq = 0; % zero alpha value should be used to obtain W_inclDM 
[f_DM3, F_DM_out3, X_out3, DM3, W_incl_DM3] = calculate_DM_global(Cal_var_new, TQs, realization, eq_w, k, back_measure, alpha_eq);

% Calculate optimized DM (in terms of calibration score) using either
% global or item weights
weight_type = 'global';
% weight_type = 'item';
tic
[F_DM_out4, X_DM_out4, DM4_opt, W_opt, W_withDM, new_alpha] = DM_Optimization(Cal_var_new, TQs, realization, k, back_measure, weight_type);
toc

% Calculation of DM using user-defined weights
user_w = [0 0 0 0 0.4 0.6 0 0 0]; % example of user defined weights. Please note that these should add up to one.
alpha_ud = 0; % zero alpha value should be used to obtain W_inclDM 
[f_DM5, F_DM_out5, X_out5, DM5, W_incl_DM5] = calculate_DM_global(Cal_var_new, TQs, realization, user_w, k, back_measure, alpha_ud);

%% Check robustness (itemwise and expertwise)

% Check robustness excluding N_max items at most
N_max_it = 5;
alpha = 0.05;
weight_type = 'global';
optimization = 'no';
incl_cal_pwr = 'no'; %if this string is 'yes' then a different calibration power is considered in the calculation of the calibration score
tic
Robustness_table = Checking_Robustness_items(Cal_var_new,...
    TQs, realization, k ,alpha, back_measure, N_max_it, weight_type, optimization, incl_cal_pwr);
toc

% Check robustness excluding N_max_exp experts at most
tic
N_max_ex = 1;
alpha = 0.05;
weight_type = 'global';
optimization = 'no';
Robustness_table_ex = Checking_Robustness_experts(Cal_var_new,...
    TQs, realization, k ,alpha, back_measure, N_max_ex, weight_type, optimization);
toc
%% Plotting individual expert's and Decision maker's assessments itemwise

% example for global weigth DM only:
% DM_str = {'c-s'};
% ystr = {'','Exp. 1','Exp. 2','Exp. 3', 'Exp. 4','Exp. 5','Exp. 6','Exp. 7', 'Exp. 8','Exp. 9', 'DM1-global','Realization',''};

% example for global, item and equal weights DMs
DM_str2 = {'c-s','g-p','r-o'};
ystr2 = {'','Exp. 1','Exp. 2','Exp. 3', 'Exp. 4','Exp. 5','Exp. 6','Exp. 7', 'Exp. 8','Exp. 9', 'DM1-global','DM2-item','DM3-equal','Realization',''};
DM_set(:,:,1) = DM1;
DM_set(:,:,2) = DM2;
DM_set(:,:,3) = DM3;

plotting_itemwise(Cal_var_new, TQs, realization, DM_set, DM_str2, ystr2)

%% Plotting box plots for robustness itemwise 
robustness_plots(Cal_var_new, Robustness_table, W_incl_DM1, N_max_it)

%% Testing alternative calculations of DMs
% In order to check the alternative calculations of the DMs uncomment the
% following lines

% intrinsic range of every item by taking into account the realization and 
% the judgments of only those experts with non-zero weights:
W = global_weights(Cal_var_new, TQs, realization, alpha, back_measure, k);
[f_DM1_alt1, F_DM1_alt1, X_1_alt1, DM1_alt1, W_incl_DM1_alt1] = alter_calc_DM_global(Cal_var_new, TQs, realization, W(:,5)', k, back_measure, alpha, 'exp_realz');

% intrinsic range of every item by taking into account only the judgments 
% of the experts with non-zero weights:
W = global_weights(Cal_var_new, TQs, realization, alpha, back_measure, k);
[f_DM1_alt2, F_DM1_alt2, X_1_alt2, DM1_alt2, W_incl_DM1_alt2] = alter_calc_DM_global(Cal_var_new, TQs, realization, W(:,5)', k, back_measure, alpha, 'exp_only');
