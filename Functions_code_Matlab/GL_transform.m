function res = GL_transform(x,v)
    % GL_transfrom: This function gives the generalised logit transform for an original time-series x
    %
    % INPUTS
    %----------------
    % x : the original time-series
    % v : a parameter, which aims to influence the evolution of variance
    %     and skweness of these distributions as a function of their mean
    %
    % OUTPUT
    %----------------
    % res : the GL transformation
    %
    %
    % Copyright P.ZAETTA 2018
   
    x = x.^v;
    res = log(x./(1-x));
end