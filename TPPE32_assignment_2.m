%% Clear all
close all
clear all
clc

%% Get Data

[Vp, TXT_Vp, RAW_Vp] = xlsread('timeSeries2018.xlsx', 'Problem 1 and 4');
dates = datenum(cell2mat(TXT_Vp(3:1442,1)));

%% Get Portfolio
Vp_C = 10;
Vp_R = zeros(1439,15);
for i=1:15
    Vp_R(1:1439,i) = Vp(2:end,i)./Vp(1:end-1,i)-1;
end

performance = [100/15 * ones(1,15)];

for i=1:length(Vp_R)
    
   performance = [performance ; performance(i)*(Vp_R(i,1:15)+1)]; 
    
end

perf=ones(1440,1);
for i=1:length(performance)
    perf=sum(performance,2);
end

semilogy(dates, perf);
title('Performance of Equal Weighted Portfolio');
datetick('x','YYYY-mm')

%% Task 1

% 1a
% Compute Return Volatility
w = 1/15;
w  = w*ones(15,1);

vol = ones(1,15);
for i=1:15
    vol(1,i) = sqrt(sum((Vp_R(1:end,i)-mean(Vp_R(1:end,i))).^2)/(length(Vp_R)-1));
end

rho = zeros(15,15);
for i=1:15
    for j=1:15
        rho(j,i)=corr(Vp_R(1:end,j),Vp_R(1:end,i));
    end
end

vol_C = ones(15,15);
for i=1:15
    vol_C(i,1:15)=vol;
end


C = vol_C'.*rho.*vol_C;

sigma_ret = sqrt(w'*C*w);

VaR_95 = sigma_ret * norminv(0.95);
VaR_975 = sigma_ret * norminv(0.975);

% 1b
% EWMA
avg_ret = mean(Vp_R,2);

x0=[0.94];
A=[1];
b=[1];
Aeq=[0];
beq=[0];
lb=0;
ub=1;

fun=@(x0)-sum((-log(EWMA_serie(x0, avg_ret)')-avg_ret(2:end).^2./EWMA_serie(x0, avg_ret)'));
[lambda, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);

EWMA = EWMA_serie(lambda,avg_ret)';

VaR_ewma95 = 10*norminv(0.95) * sqrt(EWMA(501:end));
VaR_ewma99 = 10*norminv(0.99) * sqrt(EWMA(501:end));

figure()
plot(dates(503:end), VaR_ewma95);
hold on
plot(dates(503:end), VaR_ewma99);
title('VaR-estimates using EWMA, optimal lambda')
legend('95-VaR','99-VaR')
datetick('x','YYYY-mm');

% 1c
% Historical Simulation | Expected Shortfall
figure()
subplot(2,1,1)
hist(avg_ret*Vp_C,100);
sorted_VpLoss = sort(10*avg_ret*(-1),'descend');
subplot(2,1,2)
bar(sorted_VpLoss)

% VaR
VaR95_his = zeros(length(avg_ret)-500,1);
VaR99_his = VaR95_his;
for i=1:length(VaR95_his)
    VaR95_his(i) = -quantile(avg_ret(i:i+499)*10,0.05);
    VaR99_his(i) = -quantile(avg_ret(i:i+499)*10,0.01);
end

figure()
subplot(2,1,1)
plot(dates(502:end), VaR95_his);
title('His. Sim. VaR 95%');
datetick('x','YYYY-mm')
hold on
subplot(2,1,2)
plot(dates(502:end), VaR99_his)
title('His. Sim. VaR 99%');
datetick('x','YYYY-mm')

% ES
E_fall95 = mean(sorted_VpLoss(1:floor((length(Vp_R)*0.05))));
E_fall99 = mean(sorted_VpLoss(1:floor((length(Vp_R)*0.01))));

% 1d
sigma_t = sqrt(sum((avg_ret(1:end)-mean(avg_ret(1:end))).^2)/(length(avg_ret)-1));
for i=1:length(EWMA)
    avg_retStd=avg_ret(2:end)*(sigma_t/sqrt(EWMA(i)));
end
figure();
hist(avg_retStd,100);

% VaR
VaR95_hisStd = zeros(length(avg_retStd)-500,1);
VaR99_hisStd = VaR95_hisStd;
for i=1:length(VaR95_hisStd)
    VaR95_hisStd(i) = -quantile(avg_retStd(i:i+499)*10,0.05);
    VaR99_hisStd(i) = -quantile(avg_retStd(i:i+499)*10,0.01);
end

figure()
subplot(2,1,1)
plot(dates(503:end), VaR95_hisStd);
title('His. Sim. VaR 95% STD return');
datetick('x','YYYY-mm')
hold on
subplot(2,1,2)
plot(dates(503:end), VaR99_hisStd)
title('His. Sim. VaR 99% STD return');
datetick('x','YYYY-mm')

% 1e
% FOR 1B
%XT=length(find(avg_ret(502:end)<-VaR_ewma95));
%m_=norminv(0.975,length(VaR_ewma95)*0.05,sqrt(length(VaR_ewma95)*0.05*0.95));
% EWMA
[ XT,m_ ] = testHypNor(0.05, 0.95, 2, 502, avg_ret, VaR_ewma95);
[ XT,m_ ] = testHypNor(0.01, 0.99, 2, 502, avg_ret, VaR_ewma99);
% Historical Simulation
[ XT,m_ ] = testHypNor(0.01, 0.95, 2, 501, avg_ret, VaR95_his);
[ XT,m_ ] = testHypNor(0.01, 0.99, 2, 501, avg_ret, VaR99_his);
% Historical Simulation
[ XT,m_ ] = testHypNor(0.01, 0.95, 2, 501, avg_retStd, VaR95_hisStd);
[ XT,m_ ] = testHypNor(0.01, 0.99, 2, 501, avg_retStd, VaR99_hisStd);

