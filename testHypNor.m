function [ XT,m_ ] = testHypNor( sign_lvl,conf_lvl, noSides, startReturn,serie, VaR_serie)
%testHypNor Test hypothesis for normal distribution
%

if noSides==2
    XT=length(find(serie(startReturn:end)<-VaR_serie));
    m_=norminv(1-sign_lvl/2,length(VaR_serie)*(1-conf_lvl),sqrt(length(VaR_serie)*(1-conf_lvl)*conf_lvl));
    m_=m_-norminv(sign_lvl/2,length(VaR_serie)*(1-conf_lvl),sqrt(length(VaR_serie)*(1-conf_lvl)*conf_lvl));   
end

end

