%% MASTER THESIS WIND POWER GENERATION ANALYTICS April-June 2018
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%

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

%Pn = 17.56;
Pn = 13.56;

Normalised_power = DATA(:,3)*6/1000;
Normalised_power = Normalised_power/Pn;

DATA = [DATA, Normalised_power];

%trick (with Pn equal to 13.56)
for t = 1:T
    if DATA(t, 4) >= 1
        DATA(t, 4) = 1;
    end
end

%-------------------------------------------------------------------------%
% Illustration of wind power measurements at Galicia (4 days)             %
%-------------------------------------------------------------------------%

figure(1)
plot(DATA((288*2:288*6),4), '-o');
title('\fontsize{20} Episode of 4 days with wind power measurements at Galicia');
legend('\fontsize{16} measurements');
xlabel('\fontsize{16} time steps [x10mins]');
ylabel('\fontsize{16} normalized power');
axis([0 288*4 0 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the histogram of the data                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
histogram(DATA(:,4));
title('\fontsize{20} Histogram of observed power measurements at Galicia');
xlabel('power');
ylabel('frequency');
axis tight;

%skew = skewness(DATA(:,4));
%kurt = kurtosis(DATA(:,4));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumption on the parameter v and lambda                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_lambda = 275;
v = 3.2;
lambda = 1 - (1/n_lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The generalised logit (GL) transformation                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Constraining the range of potential variations of the GL transformed    %
% variable                                                                %
%-------------------------------------------------------------------------%

threshold = 0.001;
%threshold = 0.01;

%I have to create in DATA a seventh column in order to avoid the null values
DATA(:,5) = DATA(:,4);
for t=1:T%t=1:length(DATA(:,5))
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

clear Normalised_power i Y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning and Testing datasets                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%one_third_size = round(T/3);
one_third_size = 17419;

learning_set = DATA(1:one_third_size,5:6);
testing_set = DATA(one_third_size+1:end,5:6);


%%
%-------------------------------------------------------------------------%
% Moving average and Persistence benchmarks                               %
%-------------------------------------------------------------------------%

% the lag moving average is according to the cross validation on the AR dynamic
Benchmark_MA = zeros(length(learning_set), 2);
Benchmark_Pe = zeros(length(learning_set), 2);

for t = 1:length(learning_set)-1
    Benchmark_Pe(t+1,:) = learning_set(t,1:2);
end
for t = 4:length(learning_set)
    Benchmark_MA(t,1) = mean((learning_set(t-3:t-1,1)'));
    Benchmark_MA(t,2) = mean((learning_set(t-3:t-1,2)'));
end

% matrix containing both predicted values (y and x) for each benchmark
BEN =[Benchmark_MA, Benchmark_Pe];

%%% TEST

BEN_2 = [BEN(4:end,2)-learning_set(4:end,2), BEN(4:end,4)-learning_set(4:end,2)];
BEN_2 = [zeros(3,2); BEN_2];
BEN_2 = BEN_2.^2;


% updating the scale parameter (MA & Pe benchmarks)
BEN(1:one_third_size, 5:6) = 0.9;
for t = 6:one_third_size-1
    BEN(t+1,5) = lambda^3*BEN(t-3,5)+lambda^2*(1-lambda)*BEN_2(t-2,1)+lambda*(1-lambda)*BEN_2(t-1,1)+(1-lambda)*BEN_2(t,1);
    BEN(t+1,6) = lambda^3*BEN(t-3,6)+lambda^2*(1-lambda)*BEN_2(t-2,2)+lambda*(1-lambda)*BEN_2(t-1,2)+(1-lambda)*BEN_2(t,2);
end


%-------------------------------------------------------------------------%
% Clear previous variables to avoid errors                                %
%-------------------------------------------------------------------------%

clear T1 N1 Benchmark_Pe Benchmark_MA t


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location parameter estimation via RLS with beta fixed                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Assumptions and iniatilisation                                          %
%-------------------------------------------------------------------------%

% number of lags assumed
lag = 3;
% number of iterations before one-fold cross validation
lim = 1999;
% matrix containing the parameters estimated at each step
phi_cap_1 = zeros(lag, length(learning_set));%phi_cap_1 = zeros(lag+1, length(learning_set));
% covariance matrix initialisation
cov_matrix_1 = 0.1 * eye(lag);%cov_matrix_1 = 0.1 * eye(lag+1);
% parameters initialisation
phi_1 = zeros(lag,1);%phi_1 = zeros(lag+1,1);
% the scale parameter is supposed fixed
beta_1 = 0.8;  
% vector of lag
y_1 = zeros(1,lag);%y_1 = [1, zeros(1,lag)];
%-------------------------------------------------------------------------%
% Covariance matrix updating for t = 2000                                 %
%-------------------------------------------------------------------------%

for t = lag+1:lim-1
    for j = 1:lag%j = 2:lag+1
    y_1(1,j) = learning_set(t-j,2);%y_1(1,j) = learning_set(t-j+1,2);
    end
    cov_matrix_1 = lambda*cov_matrix_1 + y_1'*y_1;
end

test_cov = cov_matrix_1; % test for the implemented matlab code
%-------------------------------------------------------------------------%
% Recursive Least Square Autoregressive (updating)                        %
%-------------------------------------------------------------------------%

for t = lim:length(learning_set)-1
    % updating the location parameter (AR-model)
    y_1 = [learning_set(t-1,2), learning_set(t-2,2), learning_set(t-3,2)];%y_1 = [1, learning_set(t-1,2), learning_set(t-2,2), learning_set(t-3,2)]; 
    epsi = learning_set(t,2) - phi_1'*y_1';
    %EPSI(t,1) = IGL_transform(epsi,v);
    %EPSI(t,1) = IGL_transform(learning_set(t,2),v) - IGL_transform(phi_1'*y_1',v);
    EPSI(t,1) = learning_set(t,2) - phi_1'*y_1';
    cov_matrix_1 = lambda*cov_matrix_1 + y_1'*y_1;
    phi_cap_1(:, t-1) = phi_1;
    phi_1 = phi_1 + inv(cov_matrix_1)*y_1'*epsi;
    mu = phi_1'*[learning_set(t,2), learning_set(t-1,2), learning_set(t-2,2)]';%mu = phi_1'*[1, learning_set(t,2), learning_set(t-1,2), learning_set(t-2,2)]';
    learning_set(t+1,3) = mu;  
end

% updating the scale parameter (AR-model)
EPSI = EPSI.^2;
learning_set(1:lim+3, 4) = beta_1;
for t = lim+3:length(learning_set)-1    
    learning_set(t+1,4) = lambda^3*learning_set(t-3,4)+lambda^2*(1-lambda)*EPSI(t-2)+lambda*(1-lambda)*EPSI(t-1)+(1-lambda)*EPSI(t);
end

%-------------------------------------------------------------------------%
% Recursive Least Square Moving Average Autoregressive (updating)         %
%-------------------------------------------------------------------------%

p = 3;
q = 1;
obj = recursiveARMA([p,q]);
obj.ForgettingFactor = lambda;
obj.InitialParameterCovariance = beta_1;

% updating the location parameter (ARMA-model)
EstimatedOutput = zeros(length(learning_set), 1);
for i = 1:length(learning_set)
    [A,C,EstimatedOutput(i)] = step(obj,learning_set(i, 2));
    %EPSI_2(i,1) = IGL_transform(learning_set(i,2),v) - IGL_transform(EstimatedOutput(i,1),v);
    EPSI_2(i,1) = learning_set(i,2) - EstimatedOutput(i,1);
end

% updating the scale parameter (ARMA-model)
EPSI_2 = EPSI_2.^2;
EstimatedOutput(1:lim+3, 2) = beta_1;
for t = lim+1:length(EstimatedOutput)-1
    EstimatedOutput(t+1,2) = lambda^3*EstimatedOutput(t-3,2)+lambda^2*(1-lambda)*EPSI_2(t-2)+lambda*(1-lambda)*EPSI_2(t-1)+(1-lambda)*EPSI_2(t);
end

%-------------------------------------------------------------------------%
% Forecast densities for each benchmark and model                         %
%-------------------------------------------------------------------------%

y_min = GL_transform(threshold, v);
y_max = GL_transform(1-threshold, v);

grid = y_min:0.1:y_max;

T1 = length(learning_set);
T = 7419;

% Moving Average GLN predictive density
% (elapsed time is 35 seconds to compute this loop)
% clear Y_MA;
tic
for t = T:T1
        Y_MA(t,:) = normcdf(grid, BEN(t,2), beta_1);
        Y_MA(t,1) = Y_MA(t,1) + normcdf(y_min, BEN(t,2), beta_1);
        Y_MA(t,end) = Y_MA(t,end) + 1 - normcdf(y_max, BEN(t,2), beta_1);
    end
toc

% Persistence GLN predictive density 
% (elapsed time is 35 seconds to compute this loop)
% clear Y_Pe;
tic
for t = T:T1
        Y_Pe(t,:) = normcdf(grid, BEN(t,4), beta_1);
        Y_Pe(t,1) = Y_Pe(t,1) + normcdf(y_min, BEN(t,4), beta_1);
        Y_Pe(t,end) = Y_Pe(t,end) + 1 - normcdf(y_max, BEN(t,4), beta_1);
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GL-Normal densities with beta fixed                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GL-Normal AR (with beta fixed) predictive density 
% (elapsed time is 35 seconds to compute this loop)
% clear Y_GLN_1;
tic
for t = T:T1
        Y_GLN_1(t,:) = normcdf(grid, learning_set(t,3), beta_1);
        Y_GLN_1(t,1) = Y_GLN_1(t,1) + normcdf(y_min, learning_set(t,3), beta_1);
        Y_GLN_1(t,end) = Y_GLN_1(t,end) + 1 - normcdf(y_max, learning_set(t,3), beta_1);
end
toc

% GL-Normal ARMA (with beta fixed) predictive density 
% (elapsed time is 35 seconds to compute this loop)
% clear Y_GLN_2;
tic
for t = T:T1
        Y_GLN_2(t,:) = normcdf(grid, EstimatedOutput(t,1), beta_1);
        Y_GLN_2(t,1) = Y_GLN_2(t,1) + normcdf(y_min, EstimatedOutput(t,1), beta_1);
        Y_GLN_2(t,end) = Y_GLN_2(t,end) + 1 - normcdf(y_max, EstimatedOutput(t,1), beta_1);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GL-Normal densities with beta not fixed                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GL-Normal AR (with beta varying) predictive density 
% (elapsed time is 38 seconds to compute this loop)
% clear Y_GLN_1;
tic
for t = T:T1
        Y_GLN_1_beta(t,:) = normcdf(grid, learning_set(t,3), learning_set(t,4));
        Y_GLN_1_beta(t,1) = Y_GLN_1_beta(t,1) + normcdf(y_min, learning_set(t,3), learning_set(t,4));
        Y_GLN_1_beta(t,end) = Y_GLN_1_beta(t,end) + 1 - normcdf(y_max, learning_set(t,3), learning_set(t,4));
end
toc

% GL-Normal ARMA (with beta varying) predictive density 
% (elapsed time is 38 seconds to compute this loop)
% clear Y_GLN_1;
tic
for t = T:T1
        Y_GLN_2_beta(t,:) = normcdf(grid, EstimatedOutput(t,1), EstimatedOutput(t,2));
        Y_GLN_2_beta(t,1) = Y_GLN_2_beta(t,1) + normcdf(y_min, EstimatedOutput(t,1), EstimatedOutput(t,2));
        Y_GLN_2_beta(t,end) = Y_GLN_2_beta(t,end) + 1 - normcdf(y_max, EstimatedOutput(t,1), EstimatedOutput(t,2));
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating of a new grid (the final grid)                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatilisation of the final grid
grid = y_min:0.1:y_max;
grid0 = y_min-1:0.1:y_min;
y_0 = zeros(1, length(grid0));
grid1 = y_max:0.1:y_max+1;
y_1 = ones(1, length(grid1));
grid_final = [grid0, grid, grid1];


% Adapting the final grid
tic %Moving-Average following a GL-Normal distribution
for t = T:T1
    Y_MA_final(t,:) = [y_0, Y_MA(t,:), y_1];
end
toc
tic %Persistence following a GL-Normal distribution
for t = T:T1
    Y_Pe_final(t,:) = [y_0, Y_Pe(t,:), y_1];
end
toc
tic %AR (with beta fixed) following a GL-Normal distribution
for t = T:T1
    Y_GLN_1_final(t,:) = [y_0, Y_GLN_1(t,:), y_1];
end
toc
tic %ARMA (with beta fixed) following a GL-Normal distribution
for t = T:T1
    Y_GLN_2_final(t,:) = [y_0, Y_GLN_2(t,:), y_1];
end
toc
tic %AR (with beta not fixed) following a GL-Normal distribution
for t = T:T1
    Y_GLN_1_beta_final(t,:) = [y_0, Y_GLN_1_beta(t,:), y_1];
end
toc
tic %ARMA (with beta not fixed) following a GL-Normal distribution
for t = T:T1
    Y_GLN_2_beta_final(t,:) = [y_0, Y_GLN_2_beta(t,:), y_1];
end
toc

[T1, T2] = size(Y_GLN_1_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncated normal pdf & cdf with beta fixed                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moving-Average (elapsed time is 90 seconds to compute this loop)
tic
for t = T:T1
        pd = makedist('norm');
        pd.sigma = beta_1;
        pd.mu = BEN(t,2);
        trunc = truncate(pd,y_min, y_max); 
        Y_MA_tr_pdf(t,:) = pdf(trunc,grid_final);
        Y_MA_tr_cdf(t,:) = cdf(trunc,grid_final);
end
toc

% Persistence (elapsed time is 90 seconds to compute this loop)
tic
for t = T:T1
        pd = makedist('norm');
        pd.sigma = beta_1;
        pd.mu = BEN(t,4);
        trunc = truncate(pd,y_min, y_max); 
        Y_Pe_tr_pdf(t,:) = pdf(trunc,grid_final); %= pdf(trunc,grid);
        Y_Pe_tr_cdf(t,:) = cdf(trunc,grid_final); %= cdf(trunc,grid);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncated normal pdf & cdf with beta not fixed                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moving-Average (elapsed time is 90 seconds to compute this loop)
tic
for t = T:T1
        pd = makedist('norm');
        pd.sigma = BEN(t,5);
        pd.mu = BEN(t,2);
        trunc = truncate(pd,y_min, y_max); 
        Y_MA_tr_pdf_2(t,:) = pdf(trunc,grid_final);
        Y_MA_tr_cdf_2(t,:) = cdf(trunc,grid_final);
end
toc

% Persistence (elapsed time is 90 seconds to compute this loop)
tic
for t = T:T1
        pd = makedist('norm');
        pd.sigma = BEN(t,6);
        pd.mu = BEN(t,4);
        trunc = truncate(pd,y_min, y_max); 
        Y_Pe_tr_pdf_2(t,:) = pdf(trunc,grid_final); 
        Y_Pe_tr_cdf_2(t,:) = cdf(trunc,grid_final);
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Countinuous Ranked Probability Score                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

true_cdf_Y = zeros(T1,T2);
for j=T:T1
    for i=1:T2
        if grid_final(1,i)<learning_set(j,2)%grid(1,i)<learning_set(j,2)
           true_cdf_Y(j,i)=0;
        else
           true_cdf_Y(j,i)=1;
        end
    end
end

%Continuous Ranked Probability Score for MA benchmark
CRPS_MA = 0;
for j=T:T1
    for i=1:T2
        if Y_MA_final(j,i)>=true_cdf_Y(j,i)
           CRPS_MA = CRPS_MA + Y_MA_final(j,i) - true_cdf_Y(j,i);
        else
           CRPS_MA = CRPS_MA - Y_MA_final(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for Persistence benchmark
CRPS_Pe = 0;
for j=T:T1
    for i=1:T2
        if Y_Pe_final(j,i)>=true_cdf_Y(j,i)
           CRPS_Pe = CRPS_Pe + Y_Pe_final(j,i) - true_cdf_Y(j,i);
        else
           CRPS_Pe = CRPS_Pe - Y_Pe_final(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for MA benchmark with (censored) Normal density with beta fixed 
CRPS_MA_tr = 0;
for j=T:T1
    for i=1:T2
        if Y_MA_tr_cdf(j,i)>=true_cdf_Y(j,i)
           CRPS_MA_tr = CRPS_MA_tr + Y_MA_tr_cdf(j,i) - true_cdf_Y(j,i);
        else
           CRPS_MA_tr = CRPS_MA_tr - Y_MA_tr_cdf(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for Persistence benchmark with (censored) Normal density with beta fixed 
CRPS_Pe_tr = 0;
for j=T:T1
    for i=1:T2
        if Y_Pe_tr_cdf(j,i)>=true_cdf_Y(j,i)
           CRPS_Pe_tr = CRPS_Pe_tr + Y_Pe_tr_cdf(j,i) - true_cdf_Y(j,i);
        else
           CRPS_Pe_tr = CRPS_Pe_tr - Y_Pe_tr_cdf(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for MA benchmark with (censored) Normal density with beta not fixed 
CRPS_MA_tr = 0;
for j=T:T1
    for i=1:T2
        if Y_MA_tr_cdf_2(j,i)>=true_cdf_Y(j,i)
           CRPS_MA_tr = CRPS_MA_tr + Y_MA_tr_cdf_2(j,i) - true_cdf_Y(j,i);
        else
           CRPS_MA_tr = CRPS_MA_tr - Y_MA_tr_cdf_2(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for Persistence benchmark with (censored) Normal density with beta not fixed 
CRPS_Pe_tr = 0;
for j=T:T1
    for i=1:T2
        if Y_Pe_tr_cdf_2(j,i)>=true_cdf_Y(j,i)
           CRPS_Pe_tr = CRPS_Pe_tr + Y_Pe_tr_cdf_2(j,i) - true_cdf_Y(j,i);
        else
           CRPS_Pe_tr = CRPS_Pe_tr - Y_Pe_tr_cdf_2(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for AR model with beta fixed
CRPS_GLN_1 = 0;
for j=T:T1
    for i=1:T2
        if Y_GLN_1_final(j,i)>=true_cdf_Y(j,i)%Y_GLN_1(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_1 = CRPS_GLN_1 + Y_GLN_1_final(j,i) - true_cdf_Y(j,i);%CRPS_GLN_1 = CRPS_GLN_1 + Y_GLN_1(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_1 = CRPS_GLN_1 - Y_GLN_1_final(j,i) + true_cdf_Y(j,i);%CRPS_GLN_1 = CRPS_GLN_1 - Y_GLN_1(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for ARMA model with beta fixed
CRPS_GLN_2 = 0;
for j=T:T1
    for i=1:T2
        if Y_GLN_2_final(j,i)>=true_cdf_Y(j,i)%Y_GLN_1(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_2 = CRPS_GLN_2 + Y_GLN_2_final(j,i) - true_cdf_Y(j,i);%CRPS_GLN_1 = CRPS_GLN_1 + Y_GLN_1(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_2 = CRPS_GLN_2 - Y_GLN_2_final(j,i) + true_cdf_Y(j,i);%CRPS_GLN_1 = CRPS_GLN_1 - Y_GLN_1(j,i) + true_cdf_Y(j,i);
        end
    end
end 

%Continuous Ranked Probability Score for AR model with beta not fixed
CRPS_GLN_1_beta = 0;
for j=T:T1
    for i=1:T2
        if Y_GLN_1_beta_final(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_1_beta = CRPS_GLN_1_beta + Y_GLN_1_beta_final(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_1_beta = CRPS_GLN_1_beta - Y_GLN_1_beta_final(j,i) + true_cdf_Y(j,i);
        end
    end
end

%Continuous Ranked Probability Score for ARMA model with beta not fixed
CRPS_GLN_2_beta = 0;
for j=T:T1
    for i=1:T2
        if Y_GLN_2_beta_final(j,i)>=true_cdf_Y(j,i)
           CRPS_GLN_2_beta = CRPS_GLN_2_beta + Y_GLN_2_beta_final(j,i) - true_cdf_Y(j,i);
        else
           CRPS_GLN_2_beta = CRPS_GLN_2_beta - Y_GLN_2_beta_final(j,i) + true_cdf_Y(j,i);
        end
    end
end



%-------------------------------------------------------------------------%
% GL-Normal density                                                       %
%-------------------------------------------------------------------------%

X = learning_set(52:end, 5); %normalized power generation estimated
Y = learning_set(51:end-1, 4); %scale parameter updated
grid2 = 0.001:0.001:0.999;

zz = (((((exp((-0.5./(Y.^2)).*((grid2-X).^2))))*v/sqrt(2*pi))./Y)./grid2)./(1-grid2.^v);

%-------------------------------------------------------------------------%
% Quantile forecasts (5% and 95%) and median values                       %
%-------------------------------------------------------------------------%

Q_0_05 = zeros(length(zz),1);
Q_0_50 = zeros(length(zz),1);
Q_0_95 = zeros(length(zz),1);

for i = 1:length(zz)
    Q_0_05(i) = quantile(zz(i,:), .05);
    Q_0_50(i) = quantile(zz(i,:), .50);
    Q_0_95(i) = quantile(zz(i,:), .95);
end

results_quantiles = [Q_0_05, Q_0_50, Q_0_95];


%-------------------------------------------------------------------------%
% Chart of wind power generation forecasted                               %
%-------------------------------------------------------------------------%

figure(5)
plot(DATA((288*2:288*6),4), '-o'); hold on;
plot(results_quantiles((288*2-51:288*6-51),2), '-o'); hold on;
title('Episode of 4 days with wind power measurements at Galicia');
legend('true values', 'forecasted values', '5% quantile forecasts', '95% quantile forecasts');
xlabel('time step [x10mins]');
ylabel('normalized power');
axis([0 288*4 0 1]);




%%
%-------------------------------------------------------------------------%
% NMAE and NRMSE (evaluation of point forecasts)                          %
%-------------------------------------------------------------------------%
%
% We can try in first with some benchmark such as 'Persistence', which
% takes the value observed at time t for t+1 (%Pn). (We can take another
% benchmark, which is 'Moving Average'. 
%-------------------------------------------------------------------------%
%
% TEST NMAE and NMRSE with Persistence and MA(3) on learning set instead on 
% the true evaluation set. (And also I don't take into account predictive 
% densities. And I evaluate on the GL transform variable instead of the 
% original variable.
%
% Evaluation of these scores overall and on a monthly basis.
%-------------------------------------------------------------------------%

