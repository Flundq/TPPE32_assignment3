function [ EWMA ] = EWMA_serie( lambda, returns )
%EWMA Generatesa EWMA Serie for the returns using lambda.
%   
    EWMA=[];
    EWMA(1) = returns(1)^2;
    for i=2:length(returns)-1
        EWMA(i)=EWMA(i-1)*lambda+returns(i)^2*(1-lambda);
    end
    
end

