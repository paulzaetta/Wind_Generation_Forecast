%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 PART MAE_RMSE_Benchmarks_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script computes the monthly and overall assessment of point forecasts
% obtained from the MA and Persistence models using the MAE and RMSE critera 
% for expectation and median respectively (with the test sample). 
%
% WARNING: elapsed time is 800 seconds to run this script
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

% Moving-Average Testing set (elapsed time is 600 seconds to run this loop)
for t = 2:T1T
        pd = makedist('norm');
        pd.sigma = BEN_T(t-1,5);
        pd.mu = BEN_T(t,2);
        trunc = truncate(pd,y_min, y_max); 
        Y_MA_tr_cdf_T(t,:) = cdf(trunc,grid_final);
end

% Persistence Testing set (elapsed time is 500 seconds to run this loop)
for t = 2:T1T
        pd = makedist('norm');
        pd.sigma = BEN_T(t-1,6);
        pd.mu = BEN_T(t,4);
        trunc = truncate(pd,y_min, y_max); 
        Y_Pe_tr_cdf_T(t,:) = cdf(trunc,grid_final);
end

%figure(1)
%plot(grid_final, Y_MA_tr_cdf_T(10000,:));
%title('\fontsize{16} Forecast CDF of truncated Normal with MA');
%axis([min(grid_final) max(grid_final) 0 1]);

%figure(2)
%plot(grid_final, Y_Pe_tr_cdf_T(10000,:));
%title('\fontsize{16} Forecast CDF of truncated Noraml with Persistence');
%axis([min(grid_final) max(grid_final) 0 1]);

%-------------------------------------------------------------------------%
% Quantiles for MA benchmark                                              %
%-------------------------------------------------------------------------%

Q_0_50_final_MA = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_MA_tr_cdf_T(i,j) < 0.5
        j = j+1;
    end
    Q_0_50_final_MA(i,1) = grid_final(1,j);
    j=1;
end

Q_0_05_final_MA = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_MA_tr_cdf_T(i,j) < 0.05
        j = j+1;
    end
    Q_0_05_final_MA(i,1) = grid_final(1,j);
    j=1;
end

Q_0_95_final_MA = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_MA_tr_cdf_T(i,j) < 0.95
        j = j+1;
    end
    Q_0_95_final_MA(i,1) = grid_final(1,j);
    j=1;
end
Quantiles_final_MA = [Q_0_05_final_MA,Q_0_50_final_MA,Q_0_95_final_MA];

%-------------------------------------------------------------------------%
% Quantiles for Persistence benchmark                                     %
%-------------------------------------------------------------------------%

Q_0_50_final_Pe = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_Pe_tr_cdf_T(i,j) < 0.5
        j = j+1;
    end
    Q_0_50_final_Pe(i,1) = grid_final(1,j);
    j=1;
end
Q_0_05_final_Pe = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_Pe_tr_cdf_T(i,j) < 0.05
        j = j+1;
    end
    Q_0_05_final_Pe(i,1) = grid_final(1,j);
    j=1;
end
Q_0_95_final_Pe = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_Pe_tr_cdf_T(i,j) < 0.95
        j = j+1;
    end
    Q_0_95_final_Pe(i,1) = grid_final(1,j);
    j=1;
end
Quantiles_final_Pe = [Q_0_05_final_Pe,Q_0_50_final_Pe,Q_0_95_final_Pe];

%-------------------------------------------------------------------------%
% Root Mean Square Error (using the median)                               %
%-------------------------------------------------------------------------%

DATA_T = DATA(one_third_size+1:end,6);

RMSE_all_tr_MA = sqrt(sum((DATA_T-Q_0_50_final_MA).^2)/T1T); %all the test period
RMSE_1_tr_MA = sqrt(sum((DATA_T(1:4338)-Q_0_50_final_MA(1:4338)).^2)/4338); %May
RMSE_2_tr_MA = sqrt(sum((DATA_T(4338:4338*2)-Q_0_50_final_MA(4338:4338*2)).^2)/4338); %June
RMSE_3_tr_MA = sqrt(sum((DATA_T(4338*2:4338*3)-Q_0_50_final_MA(4338*2:4338*3)).^2)/4338); %July
RMSE_4_tr_MA = sqrt(sum((DATA_T(4338*3:4338*4)-Q_0_50_final_MA(4338*3:4338*4)).^2)/4338); %August
RMSE_5_tr_MA = sqrt(sum((DATA_T(4338*4:4338*5)-Q_0_50_final_MA(4338*4:4338*5)).^2)/4338); %September
RMSE_6_tr_MA = sqrt(sum((DATA_T(4338*5:4338*6)-Q_0_50_final_MA(4338*5:4338*6)).^2)/4338); %October
RMSE_7_tr_MA = sqrt(sum((DATA_T(4338*6:4338*7)-Q_0_50_final_MA(4338*6:4338*7)).^2)/4338); %November
RMSE_8_tr_MA = sqrt(sum((DATA_T(4338*7:4338*8)-Q_0_50_final_MA(4338*7:4338*8)).^2)/4338); %December
RMSE_final_tr_MA = [RMSE_1_tr_MA,RMSE_2_tr_MA,RMSE_3_tr_MA,RMSE_4_tr_MA,RMSE_5_tr_MA,RMSE_6_tr_MA,RMSE_7_tr_MA,RMSE_8_tr_MA,RMSE_all_tr_MA];

RMSE_all_tr_Pe = sqrt(sum((DATA_T-Q_0_50_final_Pe).^2)/T1T); %all the test period
RMSE_1_tr_Pe = sqrt(sum((DATA_T(1:4338)-Q_0_50_final_Pe(1:4338)).^2)/4338); %May
RMSE_2_tr_Pe = sqrt(sum((DATA_T(4338:4338*2)-Q_0_50_final_Pe(4338:4338*2)).^2)/4338); %June
RMSE_3_tr_Pe = sqrt(sum((DATA_T(4338*2:4338*3)-Q_0_50_final_Pe(4338*2:4338*3)).^2)/4338); %July
RMSE_4_tr_Pe = sqrt(sum((DATA_T(4338*3:4338*4)-Q_0_50_final_Pe(4338*3:4338*4)).^2)/4338); %August
RMSE_5_tr_Pe = sqrt(sum((DATA_T(4338*4:4338*5)-Q_0_50_final_Pe(4338*4:4338*5)).^2)/4338); %September
RMSE_6_tr_Pe = sqrt(sum((DATA_T(4338*5:4338*6)-Q_0_50_final_Pe(4338*5:4338*6)).^2)/4338); %October
RMSE_7_tr_Pe = sqrt(sum((DATA_T(4338*6:4338*7)-Q_0_50_final_Pe(4338*6:4338*7)).^2)/4338); %November
RMSE_8_tr_Pe = sqrt(sum((DATA_T(4338*7:4338*8)-Q_0_50_final_Pe(4338*7:4338*8)).^2)/4338); %December
RMSE_final_tr_Pe = [RMSE_1_tr_Pe,RMSE_2_tr_Pe,RMSE_3_tr_Pe,RMSE_4_tr_Pe,RMSE_5_tr_Pe,RMSE_6_tr_Pe,RMSE_7_tr_Pe,RMSE_8_tr_Pe,RMSE_all_tr_Pe];

%-------------------------------------------------------------------------%
% Mean Absolute Error (using the median)                                  %
%-------------------------------------------------------------------------%

MAE_all_tr_MA = sum(abs(DATA_T-Q_0_50_final_MA))/T1T; %all the test period
MAE_1_tr_MA = sum(abs(DATA_T(1:4338)-Q_0_50_final_MA(1:4338)))/4338; %May
MAE_2_tr_MA = sum(abs(DATA_T(4338:4338*2)-Q_0_50_final_MA(4338:4338*2)))/4338; %June
MAE_3_tr_MA = sum(abs(DATA_T(4338*2:4338*3)-Q_0_50_final_MA(4338*2:4338*3)))/4338; %July
MAE_4_tr_MA = sum(abs(DATA_T(4338*3:4338*4)-Q_0_50_final_MA(4338*3:4338*4)))/4338; %August
MAE_5_tr_MA = sum(abs(DATA_T(4338*4:4338*5)-Q_0_50_final_MA(4338*4:4338*5)))/4338; %September
MAE_6_tr_MA = sum(abs(DATA_T(4338*5:4338*6)-Q_0_50_final_MA(4338*5:4338*6)))/4338; %October
MAE_7_tr_MA = sum(abs(DATA_T(4338*6:4338*7)-Q_0_50_final_MA(4338*6:4338*7)))/4338; %November
MAE_8_tr_MA = sum(abs(DATA_T(4338*7:4338*8)-Q_0_50_final_MA(4338*7:4338*8)))/4338; %December
MAE_final_tr_MA = [MAE_1_tr_MA, MAE_2_tr_MA, MAE_3_tr_MA, MAE_4_tr_MA, MAE_5_tr_MA, MAE_6_tr_MA, MAE_7_tr_MA, MAE_8_tr_MA, MAE_all_tr_MA];


MAE_all_tr_Pe = sum(abs(DATA_T-Q_0_50_final_Pe))/T1T; %all the test period
MAE_1_tr_Pe = sum(abs(DATA_T(1:4338)-Q_0_50_final_Pe(1:4338)))/4338; %May
MAE_2_tr_Pe = sum(abs(DATA_T(4338:4338*2)-Q_0_50_final_Pe(4338:4338*2)))/4338; %June
MAE_3_tr_Pe = sum(abs(DATA_T(4338*2:4338*3)-Q_0_50_final_Pe(4338*2:4338*3)))/4338; %July
MAE_4_tr_Pe = sum(abs(DATA_T(4338*3:4338*4)-Q_0_50_final_Pe(4338*3:4338*4)))/4338; %August
MAE_5_tr_Pe = sum(abs(DATA_T(4338*4:4338*5)-Q_0_50_final_Pe(4338*4:4338*5)))/4338; %September
MAE_6_tr_Pe = sum(abs(DATA_T(4338*5:4338*6)-Q_0_50_final_Pe(4338*5:4338*6)))/4338; %October
MAE_7_tr_Pe = sum(abs(DATA_T(4338*6:4338*7)-Q_0_50_final_Pe(4338*6:4338*7)))/4338; %November
MAE_8_tr_Pe = sum(abs(DATA_T(4338*7:4338*8)-Q_0_50_final_Pe(4338*7:4338*8)))/4338; %December
MAE_final_tr_Pe = [MAE_1_tr_Pe, MAE_2_tr_Pe, MAE_3_tr_Pe, MAE_4_tr_Pe, MAE_5_tr_Pe, MAE_6_tr_Pe, MAE_7_tr_Pe, MAE_8_tr_Pe, MAE_all_tr_Pe];

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

outcome = [RMSE_final_tr_MA; RMSE_final_tr_Pe; MAE_final_tr_MA; MAE_final_tr_Pe];
OUTCOME = array2table(outcome, 'VariableNames', {'MAY', 'JUNE', 'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER', 'ALL'}, 'RowNames', {'RMSE_MA', 'RMSE_Pe', 'MAE_MA', 'MAE_Pe'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear BEN BEN_2 BEN_2_T BEN_T Benchmark_MA_T Benchmark_Pe_T DATA DATA_T;
clear grid grid0 grid1 grid_final i j lambda learning_set lim;
clear N n_lambda one_third_size pd Pn t T T1 T1T T2 testing_set threshold;
clear trunc v y_0 y_1 Y_MA_tr_cdf_T y_max y_min Y_Pe_tr_cdf_T outcome;
clear Q_0_05_final_MA Q_0_05_final_Pe Q_0_50_final_MA Q_0_50_final_Pe; 
clear Q_0_95_final_MA Q_0_95_final_Pe Quantiles_final_MA Quantiles_final_Pe;
clear RMSE_all_tr_MA RMSE_1_tr_MA RMSE_2_tr_MA RMSE_3_tr_MA RMSE_4_tr_MA RMSE_5_tr_MA;
clear RMSE_6_tr_MA RMSE_7_tr_MA RMSE_8_tr_MA RMSE_final_tr_MA;
clear RMSE_all_tr_Pe RMSE_1_tr_Pe RMSE_2_tr_Pe RMSE_3_tr_Pe RMSE_4_tr_Pe RMSE_5_tr_Pe;
clear RMSE_6_tr_Pe RMSE_7_tr_Pe RMSE_8_tr_Pe RMSE_final_tr_Pe;
clear MAE_all_tr_MA MAE_1_tr_MA MAE_2_tr_MA MAE_3_tr_MA MAE_4_tr_MA MAE_5_tr_MA;
clear MAE_6_tr_MA MAE_7_tr_MA MAE_8_tr_MA MAE_final_tr_MA;
clear MAE_all_tr_Pe MAE_1_tr_Pe MAE_2_tr_Pe MAE_3_tr_Pe MAE_4_tr_Pe MAE_5_tr_Pe;
clear MAE_6_tr_Pe MAE_7_tr_Pe MAE_8_tr_Pe MAE_final_tr_Pe;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME contains both the RMSE and the MAE for each period considered   %
% using both the MA and the Persistence (noted Pe) benchmarks.            %
% The first and second rows are respectively the MA's RMSE and the Pe's   %
% RMSE.                                                                   %
% The third and fourth rows are respectively the MA's MAE and the Pe's    %
% MAE.                                                                    % 
%-------------------------------------------------------------------------%
toc