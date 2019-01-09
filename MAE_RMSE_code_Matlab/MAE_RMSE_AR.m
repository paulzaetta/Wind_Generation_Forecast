%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 MAE_RMSE_AR_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script computes the monthly and overall assessment of point forecasts
% obtained from the AR model using a the MAE and RMSE critera for expectation
% and median respectively (with the test sample). 
%
% WARNING: elapsed time is 25 seconds to run this script
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location parameter estimation via RLS with beta fixed                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_1 = 0.1;
BETA = ones(length(learning_set), 1)*beta_1;

p = 3; %number of lags assumed
    % creating of a System object
    obj = recursiveAR(p);
    obj.ForgettingFactor = lambda;
    obj.InitialParameterCovariance = beta_1;
    
    EstimatedOutput = zeros(length(learning_set), 1);
    for i = 1:length(learning_set)
        [A,EstimatedOutput(i,1)] = step(obj,learning_set(i, 2));
    end
    learning_set(:,3) = EstimatedOutput; %the location parameter
    
    EPSI_2 = learning_set(:,2) - EstimatedOutput;
    EPSI_2 = EPSI_2.^2;
    beta = ones(length(learning_set), 1)*beta_1;
    for t = 4:length(learning_set)
        beta(t) = lambda*beta(t-1)+(1-lambda)*EPSI_2(t);%the scale parameter
    end
BETA = [BETA, beta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT-OF-SAMPLE analysis                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BETAT = ones(length(testing_set),1);
p = 3; %number of lags assumed
    obj = recursiveAR(p); %creating of a System object
    obj.ForgettingFactor = lambda;
    obj.InitialParameterCovariance = BETA(end,2);
    
    EstimatedOutput2 = zeros(length(testing_set), 1);
    for i = 1:length(testing_set)
        [A,EstimatedOutput2(i,1)] = step(obj,testing_set(i, 2));
    end
    testing_set(:,3) = EstimatedOutput2; %the location parameter
    
    EPSI_2T = testing_set(:,2) - EstimatedOutput2;
    EPSI_2T = EPSI_2T.^2;
    betaT(1) = lambda*BETA(end,2)+(1-lambda)*EPSI_2T(1);
    for t = 2:length(testing_set)
        betaT(t) = lambda*betaT(t-1)+(1-lambda)*EPSI_2T(t);%the scale parameter
    end
BETAT = [BETAT, betaT'];

T1T = length(testing_set);

%-------------------------------------------------------------------------%
% Forecast predictive density                                             %
%-------------------------------------------------------------------------%

y_min = GL_transform(threshold, v);
y_max = GL_transform(1-threshold, v);
grid = linspace(y_min, y_max, 279);
grid_final = [y_min-1:0.1:y_min-0.1, grid, y_max+0.1:0.1:y_max+1];
T2 = length(grid_final);

Y_GLN_1 = zeros(T1T, length(grid));
for t=2:T1T
    Y_GLN_1_0(t) = normcdf(y_min, testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1_E(t) = 1 - normcdf(y_max, testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1(t,:) = normcdf(grid', testing_set(t,3), BETAT(t-1,2));
end
Y_GLN_1_0 = Y_GLN_1_0';
Y_GLN_1_E = Y_GLN_1_E';
Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;
Y_GLN_1 = [zeros(T1T,10), Y_GLN_1, ones(T1T,10)];

%figure(1)
%plot(grid_final, Y_GLN_1(10000,:));
%title('\fontsize{20} Forecast CDF of GL-Normal distribution');
%axis([min(grid_final) max(grid_final) 0 1]);

Q_0_50_final = zeros(T1T,1);
j = 1;
for i = 1:T1T
    while Y_GLN_1(i,j) < 0.5
        j = j+1;
    end
    Q_0_50_final(i,1) = grid_final(1,j);
    j=1;
end
Q_0_05_final = zeros(T1T,1);
j = 1;
for i = 1:T1T
    while Y_GLN_1(i,j) < 0.05
        j = j+1;
    end
    Q_0_05_final(i,1) = grid_final(1,j);
    j=1;
end
Q_0_95_final = zeros(T1T,1);
j = 1;
for i = 1:T1T
    while Y_GLN_1(i,j) < 0.95
        j = j+1;
    end
    Q_0_95_final(i,1) = grid_final(1,j);
    j=1;
end
Quantiles_final = [Q_0_05_final,Q_0_50_final,Q_0_95_final];

%-------------------------------------------------------------------------%
%  Root Mean Square Error (using the median)                              %
%-------------------------------------------------------------------------%

RMSE_AR_all = sqrt(sum((testing_set(:,2)-Q_0_50_final).^2)/T1T); %all the test period
RMSE_AR_1 = sqrt(sum((testing_set(1:4338,2)-Q_0_50_final(1:4338)).^2)/4338); %May
RMSE_AR_2 = sqrt(sum((testing_set(4338:4338*2,2)-Q_0_50_final(4338:4338*2)).^2)/4338); %June
RMSE_AR_3 = sqrt(sum((testing_set(4338*2:4338*3,2)-Q_0_50_final(4338*2:4338*3)).^2)/4338); %July
RMSE_AR_4 = sqrt(sum((testing_set(4338*3:4338*4,2)-Q_0_50_final(4338*3:4338*4)).^2)/4338); %August
RMSE_AR_5 = sqrt(sum((testing_set(4338*4:4338*5,2)-Q_0_50_final(4338*4:4338*5)).^2)/4338); %September
RMSE_AR_6 = sqrt(sum((testing_set(4338*5:4338*6,2)-Q_0_50_final(4338*5:4338*6)).^2)/4338); %October
RMSE_AR_7 = sqrt(sum((testing_set(4338*6:4338*7,2)-Q_0_50_final(4338*6:4338*7)).^2)/4338); %November
RMSE_AR_8 = sqrt(sum((testing_set(4338*7:4338*8,2)-Q_0_50_final(4338*7:4338*8)).^2)/4338); %December
RMSE_AR_final = [RMSE_AR_1,RMSE_AR_2,RMSE_AR_3,RMSE_AR_4,RMSE_AR_5,RMSE_AR_6,RMSE_AR_7,RMSE_AR_8,RMSE_AR_all];

%-------------------------------------------------------------------------%
% Mean Absolute Error (using the expectation)                             %
%-------------------------------------------------------------------------%

MAE_AR_all = sum(abs(testing_set(:,2)-testing_set(:,3)))/T1T; %all the test period
MAE_AR_1 = sum(abs(testing_set(1:4338,2)-testing_set(1:4338,3)))/4338; %May
MAE_AR_2 = sum(abs(testing_set(4338:4338*2,2)-testing_set(4338:4338*2,3)))/4338; %June
MAE_AR_3 = sum(abs(testing_set(4338*2:4338*3,2)-testing_set(4338*2:4338*3,3)))/4338; %July
MAE_AR_4 = sum(abs(testing_set(4338*3:4338*4,2)-testing_set(4338*3:4338*4,3)))/4338; %August
MAE_AR_5 = sum(abs(testing_set(4338*4:4338*5,2)-testing_set(4338*4:4338*5,3)))/4338; %September
MAE_AR_6 = sum(abs(testing_set(4338*5:4338*6,2)-testing_set(4338*5:4338*6,3)))/4338; %October
MAE_AR_7 = sum(abs(testing_set(4338*6:4338*7,2)-testing_set(4338*6:4338*7,3)))/4338; %November
MAE_AR_8 = sum(abs(testing_set(4338*7:4338*8,2)-testing_set(4338*7:4338*8,3)))/4338; %December
MAE_AR_final = [MAE_AR_1,MAE_AR_2,MAE_AR_3,MAE_AR_4,MAE_AR_5,MAE_AR_6,MAE_AR_7,MAE_AR_8,MAE_AR_all];

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

outcome = [RMSE_AR_final; MAE_AR_final];
OUTCOME = array2table(outcome, 'VariableNames', {'MAY', 'JUNE', 'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER', 'ALL'}, 'RowNames', {'RMSE_AR', 'MAE_AR'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear grid_final i j lambda learnin_set N n_lambda one_third_size;
clear A beta BETA beta_1 betaT BETAT DATA EPSI_2 EPSI_2T grid learning_set;
clear EstimatedOutput EstimatedOutput2 gir P Pn t T T1T T2 obj p;
clear testing_set threshold v Y_GLN_1 Y_GLN_1_0 Y_GLN_1_E y_max y_min;
clear Q_0_05_final Q_0_50_final Q_0_95_final Quantiles_final outcome;
clear RMSE_AR_all RMSE_AR_1 RMSE_AR_2 RMSE_AR_3 RMSE_AR_4 RMSE_AR_5 RMSE_AR_6 RMSE_AR_7 RMSE_AR_8 RMSE_AR_final;
clear MAE_AR_all MAE_AR_1 MAE_AR_2 MAE_AR_3 MAE_AR_4 MAE_AR_5 MAE_AR_6 MAE_AR_7 MAE_AR_8 MAE_AR_final;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME contains both the RMSE and the MAE for each period considered   %
% using the AR model (with lag 3, which is the optimal).                  %
%-------------------------------------------------------------------------%
toc