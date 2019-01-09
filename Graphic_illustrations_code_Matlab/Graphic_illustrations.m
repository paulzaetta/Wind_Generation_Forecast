%% MASTER THESIS WIND POWER GENERATION ANALYSIS April-June 2018 Graphic_Illustrations_Matlab
%%
%% ZAETTA Paul
%% Matriculation number: 872113
%%
%
% This script presents different graphic illustrations. The first illustration 
% is an episode of 4 days with wind power measurement at the Galicia firm.
% The second illustration is a histogram of the observed wind energy production at 
% Galicia, and the last one is a visual illustration of the CRPS criterion. 
%
% WARNING: elapsed time is 6 seconds to run this script
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

% trick (with Pn equal to 13.56)
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
% Visual illustration of CRPS                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% The probabilistic forecast is assumed to be a standard normal           %
% distribution.                                                           %
%-------------------------------------------------------------------------%

mu = 0;
obs = 1;
sigma_2 = 1;

grid = -4*sigma_2:0.01:4*sigma_2; 

res = normpdf(grid, mu, sigma_2);

true_cdf = zeros(length(grid),1);
for i=1:length(grid)
    if grid(i)<obs
        true_cdf(i)=0;
    else
        true_cdf(i)=1;
    end
end

true_cdf = true_cdf';

res2 = normcdf(grid, mu, sigma_2);

cdf_2 = [res2; true_cdf; grid];

y = ylim;

%-------------------------------------------------------------------------%
% Subplot for pdf and cdf illustrations                                   %
%-------------------------------------------------------------------------%

figure(3)
subplot(1,2,1)
plot(grid, res, 'b', 'linewidth', 3); hold on;
line([obs;obs],repmat(y(:),1,numel(obs)),'color','r','linestyle','--', 'linewidth', 2);
title('\fontsize{20} \bf Forecast PDF and measurement');
legend('\fontsize{10} Forecast PDF', '\fontsize{10} Measurement');
axis([min(grid) max(grid) 0 0.5]);
subplot(1,2,2)
plot(grid, res2, 'b', 'linewidth', 3); hold on; 
plot(grid, true_cdf, 'r', 'linewidth', 3);
x2 = [grid, fliplr(grid)];
inBetween = [res2, fliplr(true_cdf)];
fill(x2, inBetween, 'y');
title('\fontsize{20} \bf Forecast and measurement CDFs');
legend('\fontsize{10} Forecast CDF', '\fontsize{10} Measurement CDF', '\fontsize{10} Deviation below and above');
hold off;

%-------------------------------------------------------------------------%
toc