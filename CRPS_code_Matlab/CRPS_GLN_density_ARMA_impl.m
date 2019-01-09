%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 CRPS_ARMA_PART
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script computes the monthly and overall assessment of density forecasts
% obtained from the Autoregressive Moving Average models using a CRPS criterion. 
%
% WARNING: elapsed time is 80 seconds to run this script
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
% Location parameter estimation via RLS with beta not fixed               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 3;
beta_1 = 0.1;
BETA = ones(length(learning_set), 1)*beta_1;

for q = 1:5 %number of lags assumed
    % creating of a System object
    obj = recursiveARMA([p,q]);
    obj.ForgettingFactor = lambda;
    obj.InitialParameterCovariance = beta_1;
    
    EstimatedOutput = zeros(length(learning_set), 1);
    for i = 1:length(learning_set)
        [A,C,EstimatedOutput(i,1)] = step(obj,learning_set(i, 2));
    end
    learning_set(:,2+q) = EstimatedOutput; %the location parameter
    
    EPSI_2 = learning_set(:,2) - EstimatedOutput;
    EPSI_2 = (EPSI_2.^2);
    beta = ones(length(learning_set), 1)*beta_1;
    for t = 4:length(learning_set)
        beta(t) = lambda*beta(t-1)+(1-lambda)*EPSI_2(t);%the scale parameter
    end
BETA = [BETA, beta];
q
end

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


CRPS_GLN_1_final = zeros(1,q);
for z=1+2:q+2
    Y_GLN_1 = zeros(T1, T2);
    YY2 = zeros(T1, 1);
    clear Y_GLN_1_0 Y_GLN_1_E Y_GLN_1;
    Y_GLN_1 = zeros(T1, length(grid));
    for t=T:T1
    Y_GLN_1_0(t) = normcdf(y_min, learning_set(t,z), BETA(t-1,z-1));
    Y_GLN_1_E(t) = 1 - normcdf(y_max, learning_set(t,z), BETA(t-1,z-1));
    Y_GLN_1(t,:) = normcdf(grid', learning_set(t,z), BETA(t-1,z-1));
    end
    Y_GLN_1_0 = Y_GLN_1_0';
    Y_GLN_1_E = Y_GLN_1_E';
    Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
    Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
    Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;
    Y_GLN_1 = [zeros(T1,10), Y_GLN_1, ones(T1,10)];
    grid_final = [y_min-1:0.1:y_min-0.1, grid, y_max+0.1:0.1:y_max+1];

    CRPS_GLN_1 = 0;
    for j=T:T1
        for i=1:length(grid_final)
            if Y_GLN_1(j,i)>=true_cdf_Y(j,i)
               CRPS_GLN_1 = CRPS_GLN_1 + Y_GLN_1(j,i) - true_cdf_Y(j,i);
            else
               CRPS_GLN_1 = CRPS_GLN_1 - Y_GLN_1(j,i) + true_cdf_Y(j,i);
            end
        end
    end
CRPS_GLN_1_final(1,z-2) = CRPS_GLN_1;
z-2
end
    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT-OF-SAMPLE analysis                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 3;
BETAT = ones(length(testing_set), 1);

for q = 1:2 %number of lags assumed
    % creating of a System object
    obj = recursiveARMA([p,q]);
    obj.ForgettingFactor = lambda;
    obj.InitialParameterCovariance = BETA(end,q+1);
    
    EstimatedOutput2 = zeros(length(testing_set), 1);
    for i = 1:length(testing_set)
        [A,C,EstimatedOutput2(i,1)] = step(obj,testing_set(i, 2));
    end
    testing_set(:,2+q) = EstimatedOutput2;%the location parameter
    
    EPSI_2T = testing_set(:,2) - EstimatedOutput2;
    EPSI_2T = (EPSI_2T.^2);
    betaT(1) = lambda*BETA(end,q+1)+(1-lambda)*EPSI_2T(1);
    for t = 2:length(testing_set) 
        betaT(t) = lambda*betaT(t-1)+(1-lambda)*EPSI_2T(t);%the scale parameter
    end
BETAT = [BETAT, betaT'];
q
end

T1T = length(testing_set);

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

CRPS_GLN_1_final_testing = zeros(1,q);
for z=1+2:q+2
    clear Y_GLN_1_0 Y_GLN_1_E Y_GLN_1;
    Y_GLN_1 = zeros(T1T, length(grid));
    for t=2:T1T
    Y_GLN_1_0(t) = normcdf(y_min, testing_set(t,z), BETAT(t-1,z-1));
    Y_GLN_1_E(t) = 1 - normcdf(y_max, testing_set(t,z), BETAT(t-1,z-1));
    Y_GLN_1(t,:) = normcdf(grid', testing_set(t,z), BETAT(t-1,z-1));
    end
    Y_GLN_1_0 = Y_GLN_1_0';
    Y_GLN_1_E = Y_GLN_1_E';
    Y_GLN_1 = Y_GLN_1.*(1-Y_GLN_1_0(:));
    Y_GLN_1 = Y_GLN_1 + Y_GLN_1_0;
    Y_GLN_1(:,end) = Y_GLN_1(:,end) + Y_GLN_1_E;
    Y_GLN_1 = [zeros(T1T,10), Y_GLN_1, ones(T1T,10)];
    grid_final = [y_min-1:0.1:y_min-0.1, grid, y_max+0.1:0.1:y_max+1];

    CRPS_GLN_1_T = 0;
    for j=2:T1T
        for i=1:length(grid_final)
            CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
        end
    end
CRPS_GLN_1_final_testing(1,z-2) = CRPS_GLN_1_T;
z-2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monthly results                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CRPS_GLN_1_T = 0;
for j=1:4338
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_MAY = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*1:4338*2
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_JUNE = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*2:4338*3
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_JULY = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*2:4338*3
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_JULY = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*2:4338*3
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_JULY = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*3:4338*4
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_AUGUST = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*4:4338*5
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_SEPTEMBER = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*5:4338*6
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_OCTOBER = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*6:4338*7
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_NOVEMBER = CRPS_GLN_1_T;

CRPS_GLN_1_T = 0;
for j=4338*7:4338*8
     for i=1:length(grid_final)
         CRPS_GLN_1_T = CRPS_GLN_1_T + abs(Y_GLN_1(j,i)-true_cdf_Y(j,i));
     end
 end
CRPS_ARMA_DECEMBER = CRPS_GLN_1_T;

CRPS_ARMA_ALL = CRPS_ARMA_MAY+CRPS_ARMA_JUNE+CRPS_ARMA_JULY+CRPS_ARMA_AUGUST+CRPS_ARMA_SEPTEMBER+CRPS_ARMA_OCTOBER+CRPS_ARMA_NOVEMBER+CRPS_ARMA_DECEMBER;

CRPS_GLN_1 = [CRPS_ARMA_MAY; CRPS_ARMA_JUNE; CRPS_ARMA_JULY; CRPS_ARMA_AUGUST; CRPS_ARMA_SEPTEMBER; CRPS_ARMA_OCTOBER; CRPS_ARMA_NOVEMBER; CRPS_ARMA_DECEMBER; CRPS_ARMA_ALL];

%-------------------------------------------------------------------------%
% Final Result                                                            %
%-------------------------------------------------------------------------%

OUTCOME = array2table(CRPS_GLN_1, 'VariableNames', {'CRPS_ARMA'}, 'RowNames', {'MAY', 'JUNE', 'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER', 'ALL'});

%-------------------------------------------------------------------------%
% Clear previous useless variables                                        %
%-------------------------------------------------------------------------%

clear A ans beta BETA beta_1 betaT BETAT C CRPS_ARMA_ALL CRPS_ARMA_AUGUST CRPS_ARMA_DECEMBER CRPS_ARMA_JULY;
clear CRPS_ARMA_JUNE CRPS_ARMA_MAY CRPS_ARMA_NOVEMBER CRPS_ARMA_OCTOBER CRPS_ARMA_SEPTEMBER CRPS_GLN_1;
clear CRPS_GLN_1_final CRPS_GLN_1_final_testing CRPS_GLN_1_T DATA EPSI_2 EPSI_2T EstimatedOutput;
clear EstimatedOutput2 grid grid0 grid1 grid_final i j lambda learning_set N n_lambda obj one_third_size;
clear p Pn q t T T1 T1T T2 testing_set threshold true_cdf_Y v Y_GLN_1 Y_GLN_1_0 Y_GLN_0 Y_GLN_1_E;
clear y_0 y_1 y_max y_min YY2 z;

%-------------------------------------------------------------------------%
% RESULT:                                                                 %
%-------------------------------------------------------------------------%
% OUTCOME is a column table with 9 rows. The first row gives the          %
% ARMA(3,2)'s CRPS for the month of May, the second row gives the         %
% ARMA(3;2)'s CRPS for the month of June etc... up to the row 8, which    %
% gives the ARMA(3,2)'s CRPS for the month of December and finaly the     %
% last row gives the ARMA(3,2)'s CRPS for the whole period considered.    %
%-------------------------------------------------------------------------%
toc