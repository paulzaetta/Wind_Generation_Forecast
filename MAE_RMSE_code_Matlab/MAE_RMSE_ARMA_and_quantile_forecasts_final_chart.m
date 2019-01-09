%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 MAE_RMSE_ARMA_PART AND FINAL_FIGURE
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script computes the monthly and overall assessment of point forecasts
% obtained from the ARMA model using the MAE and RMSE criteria for expectation
% and median respectively (with the test sample). It computes also the
% final figure (prediction forecast with ARMA dynamics and GL-Normal distribution). 
%
% WARNING: elapsed time is 50 seconds to run this script
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

%We have to create in DATA a seventh column in order to avoid the null values
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
% Location parameter estimation via RLS with beta not fixed               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 3;
beta_1 = 0.1;
BETA = ones(length(learning_set), 1)*beta_1;

q=2; %number of lags assumed
obj = recursiveARMA([p,q]); %creating of a System object
obj.ForgettingFactor = lambda;
obj.InitialParameterCovariance = beta_1;
    
EstimatedOutput = zeros(length(learning_set), 1);
for i = 1:length(learning_set)
    [A,C,EstimatedOutput(i,1)] = step(obj,learning_set(i, 2));
end
learning_set(:,3) = EstimatedOutput; %the location parameter
    
EPSI_2 = learning_set(:,2) - EstimatedOutput;
EPSI_2 = (EPSI_2.^2);
beta = ones(length(learning_set), 1)*beta_1;
for t = 4:length(learning_set)
        beta(t) = lambda*beta(t-1)+(1-lambda)*EPSI_2(t); %the scale parameter
end
BETA = [BETA, beta];

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
T = 7419;

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

Y_GLN_1 = zeros(T1, length(grid));
for t=T:T1
    Y_GLN_1_0(t) = normcdf(y_min, learning_set(t,3), BETA(t-1,2));
    Y_GLN_1_E(t) = 1 - normcdf(y_max, learning_set(t,3), BETA(t-1,2));
    Y_GLN_1(t,:) = normcdf(grid', learning_set(t,3), BETA(t-1,2));
end
Y_GLN_1_0 = Y_GLN_1_0';
Y_GLN_1_E = Y_GLN_1_E';
Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;
Y_GLN_1 = [zeros(T1,10), Y_GLN_1, ones(T1,10)];
grid_final = [y_min-1:0.1:y_min-0.1, grid, y_max+0.1:0.1:y_max+1];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT-OF-SAMPLE analysis                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BETAT = ones(length(testing_set), 1);

% creating of a System object
obj = recursiveARMA([p,q]);
obj.ForgettingFactor = lambda;
obj.InitialParameterCovariance = BETA(end,2);
    
EstimatedOutput2 = zeros(length(testing_set), 1);
for i = 1:length(testing_set)
    [A,C,EstimatedOutput2(i,1)] = step(obj,testing_set(i, 2));
end
testing_set(:,3) = EstimatedOutput2; %the location parameter
    
EPSI_2T = testing_set(:,2) - EstimatedOutput2;
EPSI_2T = (EPSI_2T.^2);
betaT(1) = lambda*BETA(end,2)+(1-lambda)*EPSI_2T(1);
for t = 2:length(testing_set) 
    betaT(t) = lambda*betaT(t-1)+(1-lambda)*EPSI_2T(t);%the scale parameter
end
BETAT = [BETAT, betaT'];

T1T = length(testing_set);

Y_GLN_1 = zeros(T1T, length(grid));
for t=2:T1T
    Y_GLN_1_0(t) = normcdf(y_min, testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1_E(t) = 1 - normcdf(y_max, testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1(t,:) = normcdf(grid', testing_set(t,3), BETAT(t-1,2));
end
Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;
Y_GLN_1 = [zeros(T1T,10), Y_GLN_1, ones(T1T,10)];

%figure(1)
%plot(grid_final, Y_GLN_1(10043,:),'b','LineWidth',2);
%legend('\fontsize{20} Forecast CDF');
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
% Root Mean Square Error (using the median)                               %
%-------------------------------------------------------------------------%

RMSE_ARMA_all = sqrt(sum((testing_set(:,2)-Q_0_50_final).^2)/T1T); %all the test period
RMSE_ARMA_1 = sqrt(sum((testing_set(1:4338,2)-Q_0_50_final(1:4338)).^2)/4338); %May
RMSE_ARMA_2 = sqrt(sum((testing_set(4338:4338*2,2)-Q_0_50_final(4338:4338*2)).^2)/4338); %June
RMSE_ARMA_3 = sqrt(sum((testing_set(4338*2:4338*3,2)-Q_0_50_final(4338*2:4338*3)).^2)/4338); %July
RMSE_ARMA_4 = sqrt(sum((testing_set(4338*3:4338*4,2)-Q_0_50_final(4338*3:4338*4)).^2)/4338); %August
RMSE_ARMA_5 = sqrt(sum((testing_set(4338*4:4338*5,2)-Q_0_50_final(4338*4:4338*5)).^2)/4338); %September
RMSE_ARMA_6 = sqrt(sum((testing_set(4338*5:4338*6,2)-Q_0_50_final(4338*5:4338*6)).^2)/4338); %October
RMSE_ARMA_7 = sqrt(sum((testing_set(4338*6:4338*7,2)-Q_0_50_final(4338*6:4338*7)).^2)/4338); %November
RMSE_ARMA_8 = sqrt(sum((testing_set(4338*7:4338*8,2)-Q_0_50_final(4338*7:4338*8)).^2)/4338); %December
RMSE_ARMA_final = [RMSE_ARMA_1,RMSE_ARMA_2,RMSE_ARMA_3,RMSE_ARMA_4,RMSE_ARMA_5,RMSE_ARMA_6,RMSE_ARMA_7,RMSE_ARMA_8,RMSE_ARMA_all];

%-------------------------------------------------------------------------%
% Mean Absolute Error (using the expectation)                             %
%-------------------------------------------------------------------------%

MAE_ARMA_all = sum(abs(testing_set(:,2)-testing_set(:,3)))/T1T; %all the test period
MAE_ARMA_1 = sum(abs(testing_set(1:4338,2)-testing_set(1:4338,3)))/4338; %May
MAE_ARMA_2 = sum(abs(testing_set(4338:4338*2,2)-testing_set(4338:4338*2,3)))/4338; %June
MAE_ARMA_3 = sum(abs(testing_set(4338*2:4338*3,2)-testing_set(4338*2:4338*3,3)))/4338; %July
MAE_ARMA_4 = sum(abs(testing_set(4338*3:4338*4,2)-testing_set(4338*3:4338*4,3)))/4338; %August
MAE_ARMA_5 = sum(abs(testing_set(4338*4:4338*5,2)-testing_set(4338*4:4338*5,3)))/4338; %September
MAE_ARMA_6 = sum(abs(testing_set(4338*5:4338*6,2)-testing_set(4338*5:4338*6,3)))/4338; %October
MAE_ARMA_7 = sum(abs(testing_set(4338*6:4338*7,2)-testing_set(4338*6:4338*7,3)))/4338; %November
MAE_ARMA_8 = sum(abs(testing_set(4338*7:4338*8,2)-testing_set(4338*7:4338*8,3)))/4338; %December
MAE_ARMA_final = [MAE_ARMA_1,MAE_ARMA_2,MAE_ARMA_3,MAE_ARMA_4,MAE_ARMA_5,MAE_ARMA_6,MAE_ARMA_7,MAE_ARMA_8,MAE_ARMA_all];

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

outcome = [RMSE_ARMA_final; MAE_ARMA_final];
OUTCOME = array2table(outcome, 'VariableNames', {'MAY', 'JUNE', 'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER', 'ALL'}, 'RowNames', {'RMSE_ARMA', 'MAE_ARMA'});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of prediction forecast with ARMA dynamics and GL-Normal dist     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BETAT = ones(length(testing_set), 1);

% creating of a System object
obj = recursiveARMA([p,q]);
obj.ForgettingFactor = lambda;
obj.InitialParameterCovariance = beta_1;
    
EstimatedOutput2 = zeros(length(testing_set), 1);
for i = 1:length(testing_set)
    [A,C,EstimatedOutput2(i,1)] = step(obj,testing_set(i, 2));
end
testing_set(:,3) = EstimatedOutput2; %the location parameter
    
EPSI_2T = testing_set(:,1) - IGL_transform(EstimatedOutput2,v);
EPSI_2T = (EPSI_2T.^2);
betaT(1) = lambda*beta_1+(1-lambda)*EPSI_2T(1);
for t = 2:length(testing_set) 
    betaT(t) = lambda*betaT(t-1)+(1-lambda)*EPSI_2T(t);%the scale parameter
end
BETAT = [BETAT, betaT'];

T1T = length(testing_set);

grid = linspace(threshold, (1-threshold), 279);
T1 = length(learning_set);
T2 = length(grid);
T = 7419;
testing_set(:,3) = IGL_transform(testing_set(:,3),v);
Y_GLN_1 = zeros(T1T, length(grid));
for t=2:T1T
    Y_GLN_1_0(t) = normcdf(threshold, testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1_E(t) = 1 - normcdf((1-threshold), testing_set(t,3), BETAT(t-1,2));
    Y_GLN_1(t,:) = normcdf(grid', testing_set(t,3), BETAT(t-1,2));
end

Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;

%figure(1)
%plot(grid, Y_GLN_1(10000,:));
%title('\fontsize{20} Forecast CDF of GL-Normal distribution');
%axis([min(grid) max(grid) 0 1]);

Q_0_50_final = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_GLN_1(i,j) < 0.5
        j = j+1;
    end
    Q_0_50_final(i,1) = grid(1,j);
    j=1;
end
Q_0_05_final = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_GLN_1(i,j) < 0.05
        j = j+1;
    end
    Q_0_05_final(i,1) = grid(1,j);
    j=1;
end
Q_0_95_final = zeros(T1T,1);
j = 1;
for i = 2:T1T
    while Y_GLN_1(i,j) < 0.95
        j = j+1;
    end
    Q_0_95_final(i,1) = grid(1,j);
    j=1;
end

Quantiles_final = [Q_0_05_final,Q_0_50_final,Q_0_95_final];

figure(2)
plot(testing_set(325:755,1),'o-');hold on;
plot(Q_0_05_final(325:755),'m--'); hold on;
plot(Q_0_95_final(325:755),'m--'); hold on;
title('\fontsize{20} Episode of 32 hours with quantile forecasts of wind power measurements at Galicia');
legend('\fontsize{16} measurements','\fontsize{16} 5% and 95% quantile forecasts');
xlabel('\fontsize{16} time steps [x10mins]');
ylabel('\fontsize{16} normalized power');

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear A beta BETA beta_1 betaT BETAT C DATA EPSI_2 EPSI_2T EstimatedOutput EstimatedOutput2; 
clear grid grid0 grid1 grid_final i j lambda learning_set N n_lambda obj one_third_size p Pn q;
clear Q_0_05_final Q_0_50_final Q_0_95_final Quantiles_final t T T1 T1T T2 testing_set threshold;
clear true_cdf_Y v y_0 y_1 Y_GLN_1 Y_GLN_1_0 Y_GLN_1_E y_max y_min outcome;
clear RMSE_ARMA_all RMSE_ARMA_1 RMSE_ARMA_2 RMSE_ARMA_3 RMSE_ARMA_4 RMSE_ARMA_5 RMSE_ARMA_6 RMSE_ARMA_7 RMSE_ARMA_8 RMSE_ARMA_final;
clear MAE_ARMA_all MAE_ARMA_1 MAE_ARMA_2 MAE_ARMA_3 MAE_ARMA_4 MAE_ARMA_5 MAE_ARMA_6 MAE_ARMA_7 MAE_ARMA_8 MAE_ARMA_final;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME contains both the RMSE and the MAE for each period considered   %
% using the ARMA model (with AR lag equal to 3 and MA lag equal to 2,     %
% which are the optimal lags).                                            %
%-------------------------------------------------------------------------%
toc