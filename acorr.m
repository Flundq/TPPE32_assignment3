function acorr = acorr( serie, lag )
% acorr Determines autocorrelation
   
    acorr = corr(serie(1:length(serie)-lag), serie(1+lag:length(serie)));

end

