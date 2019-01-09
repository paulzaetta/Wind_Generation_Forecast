function res = IGL_transform(x,v)
    % IGL_transfrom: This function gives the invers generalised logit transform for an original time-series x
    %
    % INPUTS
    %----------------
    % x : the original time-series
    % v : a parameter, which aims to influence the evolution of variance
    %     and skweness of these distributions as a function of their mean
    %
    % OUTPUT
    %----------------
    % res : the IGL transformation
    %
    %
    % Copyright P.ZAETTA 2018

    res = (1+(1./exp(x))).^(-1/v);
end