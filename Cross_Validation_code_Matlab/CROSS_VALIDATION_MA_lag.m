%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 CROSS_VALIDATION_MA_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% Using one-fold cross-validation over the learning period in order to
% select the optimal moving average order for the dynanamic models. 
%
% WARNING: elapsed time is 52 seconds to run this script
%%
tic
clc;
clear all;
format long;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE DATASET                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA = xlsread('/Users/paulzaetta/Documents/MATLAB/Matlab_Thesis/Data_2016.xlsx');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We drop the missing values                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA(any(isnan(DATA),2),:) = [];

[T,N] = size(DATA);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We normalized the power average by Pn (nominal capacity)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pn = 13.56;

Normalised_power = DATA(:,3)*6/1000;
Normalised_power = Normalised_power/Pn;

DATA = [DATA, Normalised_power];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption on the parameter v and lambda                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_lambda = 2500;
v = 3.2;
lambda = 1 - (1/n_lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The generalised logit (GL) transformation                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Contraining the range of potential variations of the GL transformed     %
% variable                                                                %
%-------------------------------------------------------------------------%

threshold = 0.001;

% We have to create in DATA a seventh column in order to avoid the null values
DATA(:,5) = DATA(:,4);
for t=1:length(DATA(:,5))
    if DATA(t,5) <= threshold
        DATA(t,5) = threshold;
    elseif DATA(t,5) >= 1 - threshold
        DATA(t,5) = 1 - threshold;
    end
end

Y = DATA(:,5).^v;
Y = log(Y./(1-Y));
DATA = [DATA, Y];

%-------------------------------------------------------------------------%
% Clear previous variables to avoid errors                                %
%-------------------------------------------------------------------------%

clear Normalised_power_2 Normalised_power i Y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning and Testing datasets                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

one_third_size = round(T/3);

learning_set = DATA(1:one_third_size,5:6);
testing_set = DATA(one_third_size+1:end,5:6);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal number of lags for the MA part (without constant)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location parameter estimation via RLS with beta fixed                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 3;
beta_1 = 0.1;
BETA = ones(length(learning_set), 1)*beta_1;

for q = 0:6 %number of lags assumed
    % creating of a System object
    obj = recursiveARMA([p,q]);
    obj.ForgettingFactor = lambda;
    obj.InitialParameterCovariance = beta_1;
    
    EstimatedOutput = zeros(length(learning_set), 1);
    for i = 1:length(learning_set)
        [A,C,EstimatedOutput(i,1)] = step(obj,learning_set(i, 2));
    end
    learning_set(:,2+1+q) = EstimatedOutput; %the location parameter
    
    EPSI_2 = learning_set(:,2) - EstimatedOutput;
    EPSI_2 = (EPSI_2.^2);
    beta = ones(length(learning_set), 1)*beta_1;
    for t = 4:length(learning_set)-1
        beta(t+1) = lambda^3*beta(t-3)+lambda^2*(1-lambda)*EPSI_2(t-2)+lambda*(1-lambda)*EPSI_2(t-1)+(1-lambda)*EPSI_2(t); %the scale parameter
    end
BETA = [BETA, beta];
end

%-------------------------------------------------------------------------%
% Forecast densities for each ARMA(p) model with and without intercept    %
%-------------------------------------------------------------------------%

y_min = GL_transform(threshold, v);
y_max = GL_transform(1-threshold, v);

grid = linspace(y_min, y_max, 279);
grid0 = y_min-1:0.099:y_min-0.099;
y_0 = zeros(1, length(grid0));
grid1 = 0.099+y_max:0.099:y_max+1;
y_1 = ones(1, length(grid1));
grid_final = [grid0, grid, grid1];

T1 = length(learning_set);
T2 = length(grid_final);

% GL-Normal ARMA (with beta fixed) predictive density 
true_cdf_Y = zeros(T1,T2);
for j=7000:T1
    for i=1:T2
        if grid_final(1,i)<learning_set(j,2)
           true_cdf_Y(j,i)=0;
        else
           true_cdf_Y(j,i)=1;
        end
    end
end

CRPS_GLN_2 = zeros(1,q+1);
for z=1+2:q+1+2
    Y_GLN_2 = zeros(T1, T2);
for t = 7000:T1
        Y_GLN_2(t,:) = normcdf(grid_final, learning_set(t,z), BETA(t,z-1));
        Y_GLN_2(t,length(grid1)+1) = Y_GLN_2(t,length(grid1)+1) + normcdf(y_min, learning_set(t,z), BETA(t,z-1));
        Y_GLN_2(t,length(grid)+length(grid1)) = Y_GLN_2(t,length(grid)+length(grid1)) + 1 - normcdf(y_max, learning_set(t,z), BETA(t,z-1));
end
for j=7000:T1
    for i=1:T2
        if Y_GLN_2(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_2(1,z-2) = CRPS_GLN_2(1,z-2) + Y_GLN_2(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_2(1,z-2) = CRPS_GLN_2(1,z-2) - Y_GLN_2(j,i) + true_cdf_Y(j,i);
        end
    end
end
end

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

OUTCOME = array2table(CRPS_GLN_2, 'VariableNames', {'MA_lag_0', 'MA_lag_1', 'MA_lag_2', 'MA_lag_3', 'MA_lag_4', 'MA_lag_5', 'MA_lag_6'}, 'RowNames', {'CRPS'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear A beta BETA beta_1 C CRPS_GLN_2 DATA EPSI_2 EstimatedOutput grid grid0;
clear grid1 grid_final i j lambda learning_set N n_lambda obj one_third_size;
clear p Pn q t T T1 T2 testing_set threshold true_cdf_Y v y_0 y_1 Y_GLN_2;
clear y_max y_min z;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME is a row table with 7 columns. The first column gives the CRPS  % 
% of lag 0, the second column gives the CRPS of lag 1, the third column   %
% gives the CRPS of lag 2 etc... up to the seventh column which gives the % 
% CRPS of lag 6 (all without intercept).                                  %
%-------------------------------------------------------------------------%
% The optimal MA lag is two according to the cross-validation exercise    %
% which gives the smallest CRPS (column 3).                               %
%-------------------------------------------------------------------------%
toc