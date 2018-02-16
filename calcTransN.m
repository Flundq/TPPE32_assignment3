function [ testStatistic, Fscore ] = calcTransN( returns, VaR, sign_lvl)
%calcTransN Calculate VaR in/out transistions for return serie
n11 = 0; n10 = n11; n01 = n10; n00 = n01;
if returns(1) < -VaR(1)
    startState=1;
    n01=n01+1;
else
    startState=0;
    n00=n00+1;
end

for i=2:length(VaR)
    if returns(i) < -VaR(i)
        if returns(i-1) < -VaR(i-1)
           n11=n11+1; 
        end
        n01=n01+1;
    else
        if returns(i-1) < -VaR(i-1)
           n10=n10+1; 
        end
        n00=n00+1;
    end
end

PI=(n01+11)/(n00+n01+n10+n11);
PI01=n01/(n00+n01);
PI11=n11/(n10+n11);

testStatistic=-2*log((1-PI)^(n00+n10)*PI^(n01+n11))+2*log((1-PI01)^(n00)*PI01^(n01)*(1-PI11)^(n10)*PI11^(n11));

Fscore = chi2inv(1-sign_lvl,1);

end

