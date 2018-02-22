function [ XT,m_,Z,N ] = testHypNor( sign_lvl,conf_lvl, noSides,serie, VaR_serie)
%testHypNor Test hypothesis for normal distribution
%

if noSides==2
    XT=length(find(serie < -VaR_serie));
    m_=norminv(1-sign_lvl/2, length(VaR_serie)*(1-conf_lvl), sqrt(length(VaR_serie)*(1-conf_lvl)*conf_lvl));
    m_=m_ - norminv(sign_lvl/2, length(VaR_serie)*(1-conf_lvl), sqrt(length(VaR_serie)*(1-conf_lvl)*conf_lvl));
    
    Z = abs((XT-length(VaR_serie)*(1-conf_lvl))/(sqrt(length(VaR_serie)*conf_lvl*(1-conf_lvl))));
    N = norminv(1-sign_lvl/2);
    
end

end

