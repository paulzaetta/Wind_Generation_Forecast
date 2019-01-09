%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 CRPS_MA_Pe_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script computes the monthly and overall assessment of density forecasts
% obtained from the MA and Persistence models using a CRPS criterion. 
%
% WARNING: elapsed time is 1300 seconds to run this script
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

%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning and Testing datasets                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Benchmark_MA_T = zeros(length(testing_set), 2);
Benchmark_Pe_T = zeros(length(testing_set), 2);

for t = 1:length(learning_set)-1
    Benchmark_Pe(t+1,:) = learning_set(t,1:2);
end
for t = 1:length(testing_set)-1
    Benchmark_Pe_T(t+1,:) = testing_set(t,1:2);
end
for t = 4:length(learning_set)
    Benchmark_MA(t,1) = mean((learning_set(t-3:t-1,1)'));
    Benchmark_MA(t,2) = mean((learning_set(t-3:t-1,2)'));
end
for t = 4:length(testing_set)
    Benchmark_MA_T(t,1) = mean((testing_set(t-3:t-1,1)'));
    Benchmark_MA_T(t,2) = mean((testing_set(t-3:t-1,2)'));
end

% matrix containing both predicted values (y and x) for each benchmark
BEN =[Benchmark_MA, Benchmark_Pe];
BEN_T =[Benchmark_MA_T, Benchmark_Pe_T];

BEN_2 = [learning_set(4:end,2)-BEN(4:end,2), learning_set(4:end,2)-BEN(4:end,4)];
BEN_2 = [zeros(3,2); BEN_2];
BEN_2 = BEN_2.^2;
BEN_2_T = [testing_set(4:end,2)-BEN_T(4:end,2), testing_set(4:end,2)-BEN_T(4:end,4)];
BEN_2_T = [zeros(3,2); BEN_2_T];
BEN_2_T = BEN_2_T.^2;

lim = 1999;
% updating the scale parameter (MA & Pe benchmarks)
BEN(1:one_third_size, 5:6) = 0.1; %initialisation of the scale parameter
for t = lim+3:one_third_size
    BEN(t,5) = lambda*BEN(t-1,5)+(1-lambda)*BEN_2(t,1);
    BEN(t,6) = lambda*BEN(t-1,6)+(1-lambda)*BEN_2(t,2);
end
BEN_T(1,5) = lambda*BEN(end,5)+(1-lambda)*BEN_2_T(1,1); %initialisation of the scale parameter
BEN_T(1,6) = lambda*BEN(end,6)+(1-lambda)*BEN_2_T(1,2); %initialisation of the scale parameter
for t = 2:length(testing_set)
    BEN_T(t,5) = lambda*BEN_T(t-1,5)+(1-lambda)*BEN_2_T(t,1);
    BEN_T(t,6) = lambda*BEN_T(t-1,6)+(1-lambda)*BEN_2_T(t,2);
end

%-------------------------------------------------------------------------%
% Clear previous variables to avoid errors                                %
%-------------------------------------------------------------------------%

clear T1 N1 Benchmark_Pe Benchmark_MA t

%-------------------------------------------------------------------------%

y_min = GL_transform(threshold, v);
y_max = GL_transform(1-threshold, v);
T1 = length(learning_set);
T = 7419;

% Iniatilisation of the final grid
grid = y_min:0.1:y_max;
grid0 = y_min-1:0.1:y_min;
y_0 = zeros(1, length(grid0));
grid1 = y_max:0.1:y_max+1;
y_1 = ones(1, length(grid1));
grid_final = [grid0, grid, grid1];
T2 = length(grid_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncated normal pdf & cdf with beta fixed                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1T = length(testing_set);

% Moving-Average Testing set (elapsed time is 90 seconds to compute this loop)
for t = 2:T1T
        pd = makedist('norm');
        pd.sigma = BEN_T(t-1,5);
        pd.mu = BEN_T(t,3);
        trunc = truncate(pd,y_min, y_max); 
        Y_MA_tr_pdf_T(t,:) = pdf(trunc,grid_final);
        Y_MA_tr_cdf_T(t,:) = cdf(trunc,grid_final);
        t
end

% Persistence Testing set (elapsed time is 90 seconds to compute this loop)
for t = 2:T1T
        pd = makedist('norm');
        pd.sigma = BEN_T(t-1,6);
        pd.mu = BEN_T(t,4);
        trunc = truncate(pd,y_min, y_max); 
        Y_Pe_tr_pdf_T(t,:) = pdf(trunc,grid_final);
        Y_Pe_tr_cdf_T(t,:) = cdf(trunc,grid_final);
        t
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Countinuous Ranked Probability Score                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

true_cdf_Y = zeros(T1T,T2);
for j=2:T1T
    for i=1:T2
        if grid_final(1,i)<testing_set(j,2)
           true_cdf_Y(j,i)=0;
        else
           true_cdf_Y(j,i)=1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING SET CRPS                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CRPS_matrix = zeros(9,2); %initialisation

%-------------------------------------------------------------------------%
% CRPS from the first period                                              %
%-------------------------------------------------------------------------%

CRPS_MA_tr_T = 0;
for j=1:4338*1
    for i=1:T2
        CRPS_MA_tr_T = CRPS_MA_tr_T + abs(Y_MA_tr_cdf_T(j,i)-true_cdf_Y(j,i));
    end
end
CRPS_matrix(1,1) = CRPS_MA_tr_T;

CRPS_Pe_tr_T = 0;
for j=1:4338*1
    for i=1:T2
        CRPS_Pe_tr_T = CRPS_Pe_tr_T + abs(Y_Pe_tr_cdf_T(j,i)-true_cdf_Y(j,i));
    end
end
CRPS_matrix(1,2) = CRPS_Pe_tr_T;

%-------------------------------------------------------------------------%
% CRPS from the second period to the last period                          %
%-------------------------------------------------------------------------%

for ii = 2:8
% Continuous Ranked Probability Score for MA benchmark with (censored) Normal density with beta not fixed 
CRPS_MA_tr_T = 0;
for j=4338*(ii-1):4338*ii
    for i=1:T2
        CRPS_MA_tr_T = CRPS_MA_tr_T + abs(Y_MA_tr_cdf_T(j,i)-true_cdf_Y(j,i));
    end
end
CRPS_matrix(ii, 1) = CRPS_MA_tr_T;
% Continuous Ranked Probability Score for Persistence benchmark with (censored) Normal density with beta not fixed 
CRPS_Pe_tr_T = 0;
for j=4338*(ii-1):4338*ii
    for i=1:T2
        CRPS_Pe_tr_T = CRPS_Pe_tr_T + abs(Y_Pe_tr_cdf_T(j,i)-true_cdf_Y(j,i));
    end
end

CRPS_matrix(ii, 2) = CRPS_Pe_tr_T; 
end

%-------------------------------------------------------------------------%
% CRPS for the whole period considered                                    %
%-------------------------------------------------------------------------%

CRPS_matrix(end, 1) = sum(CRPS_matrix(:,1));
CRPS_matrix(end, 2) = sum(CRPS_matrix(:,2));

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

OUTCOME = array2table(CRPS_matrix, 'VariableNames', {'MA_benchmark', 'Pe_benchmark'}, 'RowNames', {'CRPS_MAY','CRPS_JUNE','CRPS_JULY','CRPS_AUGUST','CRPS_SEPTEMBER','CRPS_OCTOBER', 'CRPS_NOVEMBER', 'CRPS_DECEMBER', 'CRPS_ALL'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear BEN BEN_2 BEN_2_T BEN_T Benchmark_MA_T Benchmark_Pe_T DATA grid grid0 grid1;
clear grid_final i ii j lambda learning_set lim N n_lambda one_third_size pd Pn t;
clear T T1 T1T T2 testing_set threshold true_cdf_Y trunc v y_0 y_1 Y_MA_tr_cdf_T;
clear Y_MA_tr_pdf_T y_max y_min Y_Pe_tr_cdf_T Y_Pe_tr_pdf_T;
clear CRPS_MA_tr CRPS_MA_tr_T CRPS_Pe_tr_T CRPS_matrix;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME is a table with 9 rows and 2 columns. The first column gives    %
% the CRPS of the MA benchmark and the second column gives the CRPS of    %  
% the Persistence benchmark. Each row gives the CRPS for each period.     %
%-------------------------------------------------------------------------%
% We note that the CRPS of the Pe benchmark is smaller (and therefore     %
% better) for the period considered than the MA benchmark.                %
%-------------------------------------------------------------------------%
toc