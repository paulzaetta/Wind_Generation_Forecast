%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 CROSS_VALIDATION_AR_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% Using one-fold cross-validation over the learning period in order to
% select the optimal autoregressive order for the dynanamic models.
%
% WARNING: elapsed time is 930 seconds to run this script
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

one_third_size = 17419; % approximately equal to round(T/3)

learning_set = DATA(1:one_third_size,5:6);
testing_set = DATA(one_third_size+1:end,5:6);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal number of lags for the AR part                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Location parameter estimation via RLS with beta fixed                   %
%-------------------------------------------------------------------------%

for lag = 1:4 %number of lags assumed
    % number of iterations before one-fold cross validation
    lim = 1999;
    % matrix containing the parameters estimated at each step
    phi_cap_1_int = zeros(lag+1, length(learning_set));
    phi_cap_1 = zeros(lag, length(learning_set));
    % covariance matrix initialisation
    cov_matrix_1_int = 0.1 * eye(lag+1);
    cov_matrix_1 = 0.1 * eye(lag);
    % parameters initialisation
    phi_1_int = zeros(lag+1,1);
    phi_1 = zeros(lag,1);
    % the scale parameter is supposed fixed
    beta_1 = 0.1;  
    % vector of lag
    y_1_int = [1, zeros(1,lag)];
    y_1_int_2 = [1, zeros(1,lag)];
    y_1 = zeros(1,lag);
    y_1_2 = zeros(1,lag);

%-------------------------------------------------------------------------%
% Covariance matrix updating for t = 2000                                 %
%-------------------------------------------------------------------------%

    for t = lag+1:lim
        for j = 2:lag+1
            y_1_int(1,j) = learning_set(t-j+1,2);
        end
            cov_matrix_1_int = lambda*cov_matrix_1_int + y_1_int'*y_1_int;
    end
    
    for t = lag+1:lim
        for j = 1:lag
            y_1(1,j) = learning_set(t-j,2);
        end
            cov_matrix_1 = lambda*cov_matrix_1 + y_1'*y_1;
    end
    
%-------------------------------------------------------------------------%
% Recursive Least Square (updating)                                       %
%-------------------------------------------------------------------------%

        for t = lim:length(learning_set)
            for j = 2:lag+1
                y_1_int(1,j) = learning_set(t-j+1,2);
            end
            epsi_int = learning_set(t,2) - phi_1_int'*y_1_int';
            cov_matrix_1_int = lambda*cov_matrix_1_int + y_1_int'*y_1_int;
            phi_cap_1_int(:, t-1) = phi_1_int;
            phi_1_int = phi_1_int + cov_matrix_1_int\y_1_int'*epsi_int;
            for j = 2:lag+1
                y_1_int_2(1,j) = learning_set(t-j+2,2);
            end
            mu_int = phi_1_int'*y_1_int_2';
            learning_set(t+1,2*lag+1) = mu_int;     
            for j = 1:lag
                y_1(1,j) = learning_set(t-j,2);
            end
            epsi = learning_set(t,2) - phi_1'*y_1';
            cov_matrix_1 = lambda*cov_matrix_1 + y_1'*y_1;
            phi_cap_1(:, t-1) = phi_1;
            phi_1 = phi_1 + cov_matrix_1\y_1'*epsi;
            for j = 1:lag
                y_1_2(1,j) = learning_set(t-j+1,2);
            end
            mu = phi_1'*y_1_2';
            learning_set(t+1,2*lag+2) = mu;         
        end
end

%-------------------------------------------------------------------------%
% Forecast densities for each AR(p) model with and without intercept      %
%-------------------------------------------------------------------------%

y_min = GL_transform(threshold, v);
y_max = GL_transform(1-threshold, v);

grid = linspace(y_min, y_max, 279);
grid0 = y_min-1:0.099:y_min-0.099;
y_0 = zeros(1, length(grid0));
grid1 = 0.099+y_max:0.099:y_max+1;
y_1 = ones(1, length(grid1));
grid_final = [grid0, grid, grid1];

T1 = 17419;
T2 = length(grid_final);
T = 7419;

% GL-Normal AR (with beta fixed) predictive density 
true_cdf_Y = zeros(T1,T2);
for j=T:T1
    for i=1:T2
        if grid_final(1,i)<learning_set(j,2)
           true_cdf_Y(j,i)=0;
        else
           true_cdf_Y(j,i)=1;
        end
    end
end

CRPS_GLN_1 = zeros(1,2*lag);
for z=1+2:lag*2+2
    Y_GLN_1 = zeros(T1, T2);
    YY2 = zeros(T1, 1);
for t = T:T1
        YY1 = normcdf(y_min, learning_set(t,z), beta_1);
        Y_GLN_1(t,1) = YY1 + (1-YY1)*normcdf(grid(1), learning_set(t,z), beta_1);
end
for t = T:T1
    for i = 2:length(grid)
        Y_GLN_1(t,i) = Y_GLN_1(t,1) + (1-Y_GLN_1(t,1))*normcdf(grid(i), learning_set(t,z), beta_1);
    end
end
for t = T:T1
    YY2(t,1) = 1 - normcdf(y_max, learning_set(t,z), beta_1);
end  
Y_GLN_1(:,end) = Y_GLN_1(:,end) + YY2;
Y_GLN_1 = [zeros(T1,length(grid0)), Y_GLN_1, ones(T1,length(grid1))];

for j=T:T1
    for i=1:T2
        if Y_GLN_1(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_1(1,z-2) = CRPS_GLN_1(1,z-2) + Y_GLN_1(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_1(1,z-2) = CRPS_GLN_1(1,z-2) - Y_GLN_1(j,i) + true_cdf_Y(j,i);
        end
    end
end
end

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

OUTCOME = array2table(CRPS_GLN_1, 'VariableNames', {'AR_lag_1_with_intercept', 'AR_lag_1_no_intercept', 'AR_lag_2_with_intercept', 'AR_lag_2_no_intercept', 'AR_lag_3_with_intercept', 'AR_lag_3_no_intercept', 'AR_lag_4_with_interncept', 'AR_lag_4_no_intercept'}, 'RowNames', {'CRPS'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear beta_1 cov_matrix_1 cov_matrix_1_int CRPS_GLN_1 DATA epsi epsi_int grid grid0 grid1;
clear grid_final i j lag lambda learning_set lim mu mu_int N n_lambda one_third_size;
clear phi_1 phi_1_int phi_cap_1 phi_cap_1_int Pn t T T1 T2 testing_set threshold;
clear true_cdf_Y v y_0 y_1 y_1_2 y_1_int y_1_int_2 Y_GLN_1 y_max y_min YY1 YY2 z;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME is a row table with 8 columns. The first column gives the       %
% CRPS of lag 1 with intercept, the second column gives the CRPS of lag 1 %
% without intercept, the third column gives the CRPS of lag 2 with        %
% intercept etc... up to the eighth column which gives the CRPS of lag 4  %
% without intercept.                                                      %
%-------------------------------------------------------------------------%
% The optimal AR lag is three without intercept (column 6, which gives    %
% the smallest CRPS) according to the cross-validation exercise.          %
%-------------------------------------------------------------------------%
toc