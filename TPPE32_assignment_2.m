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

performance = [10/15 * ones(1,15)];

for i=1:length(Vp_R)
    
   performance = [performance ; performance(i,:).*(Vp_R(i,1:15)+1)]; 
    
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

% x0=[0.94];
% A=[1];
% b=[1];
% Aeq=[0];
% beq=[0];
% lb=0;
% ub=1;
% 
% fun=@(x0)-sum((-log(EWMA_serie(x0, avg_ret)')-avg_ret(2:end).^2./EWMA_serie(x0, avg_ret)'));
% [lambda, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);

lambda=0.94;

EWMA = EWMA_serie(lambda,avg_ret)';

VaR_ewma95 = norminv(0.95) * sqrt(EWMA(501:end));
VaR_ewma99 = norminv(0.99) * sqrt(EWMA(501:end));

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
    VaR95_his(i) = -quantile(avg_ret(i:i+499),0.05);
    VaR99_his(i) = -quantile(avg_ret(i:i+499),0.01);
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
    VaR95_hisStd(i) = -quantile(avg_retStd(i:i+499),0.05);
    VaR99_hisStd(i) = -quantile(avg_retStd(i:i+499),0.01);
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
% EWMA
[ XT1,m_1,Z_1,N_1 ] = testHypNor(0.05, 0.95, 2, avg_ret(502:end), VaR_ewma95);
[ XT2,m_2,Z_1,N_2 ] = testHypNor(0.05, 0.99, 2, avg_ret(502:end), VaR_ewma99);
% Historical Simulation
[ XT3,m_3,Z_3,N_3 ] = testHypNor(0.05, 0.95, 2, avg_ret(501:end), VaR95_his);
[ XT4,m_4,Z_4,N_4 ] = testHypNor(0.05, 0.99, 2, avg_ret(501:end), VaR99_his);
% Historical Simulation STD
[ XT5,m_5,Z_5,N_5 ] = testHypNor(0.05, 0.95, 2, avg_retStd(501:end), VaR95_hisStd);
[ XT6,m_6,Z_6,N_6 ] = testHypNor(0.05, 0.99, 2, avg_retStd(501:end), VaR99_hisStd);
output.struct.hTest=[m_1 m_2 m_3 m_4 m_5 m_6]-[XT1 XT2 XT3 XT4 XT5 XT6];

% 1f
% TEST EWMA
[ TS1, FS1 ] = calcTransN( avg_ret(500:end), VaR_ewma95, 0.05)
[ TS2 ] = calcTransN( avg_ret(500:end), VaR_ewma99, 0.05)
% Historical Simulation
[ TS3 ] = calcTransN( avg_ret(500:end), VaR95_his, 0.05)
[ TS4 ] = calcTransN( avg_ret(500:end), VaR99_his, 0.05)
% Historical Simulation STD
[ TS5] = calcTransN( avg_retStd(500:end), VaR95_hisStd, 0.05)
[ TS6] = calcTransN( avg_retStd(500:end), VaR99_hisStd, 0.05)
output.struct.chris=FS1*ones(1,6)-[TS1 TS2 TS3 TS4 TS5 TS6];

figure()
scatter(dates(1+501:end), avg_ret(501:end),'.')
hold on
plot(dates(2+501:end), -VaR_ewma95)


%% Task 2 in EXCEL
% clear
% close all
% clc
% [Bn, TXT_Bn, RAW_Bn] = xlsread('timeSeries2018.xlsx', 'Problem 2');
% cfDates=[datenum(cell2mat(RAW_Bn(4:8,32))), datenum(cell2mat(RAW_Bn(4:8,33)))] ;
% vertices=[1 3 6 12 24 36 48 60 84 108 120 180 240 360];  
% mTillExp = ceil((cfDates(end,1)-datenum(date()))/30);
% vertNeed=(length(find(mTillExp>vertices))+1);

%% Task 3
[Op, TXT_Op, RAW_Op] = xlsread('timeSeries2018.xlsx', 'Problem 3');
dates=datenum(cell2mat(RAW_Op(3:3687,1)));
SPX=cell2mat(RAW_Op(3:3687,2));
VIX=cell2mat(RAW_Op(3:3687,3))/100;
RF3m=cell2mat(RAW_Op(3:3687,4))/100;
RF3m = log(1+RF3m*0.25)/(0.25);

K_C16m = 2765;
K_C20a = 2765;
K_P16m = 2775;

iv_C16m = Op(1,7)/100;
iv_C20a = Op(3,7)/100;
iv_P16m = Op(2,7)/100;

expiry_C16m = datenum(TXT_Op(7,7));
expiry_P16m = datenum(TXT_Op(8,7)');
expiry_C20a = datenum(TXT_Op(9,7));

if ismac
   dates=dates+693960;
   'It is a MAC'
end

C16m_d1 = (log(SPX(1)/K_C16m)+(RF3m(1)+(iv_C16m^2)/2)*((expiry_C16m-dates(1))/365))/(iv_C16m*sqrt(((expiry_C16m-dates(1))/365)));
C16m_d2 = C16m_d1-iv_C16m*sqrt(((expiry_C16m-dates(1))/365));

C20a_d1 = (log(SPX(1)/K_C20a)+(RF3m(1)+(iv_C20a^2)/2)*((expiry_C20a-dates(1))/365))/(iv_C20a*sqrt(((expiry_C20a-dates(1))/365)));
C20a_d2 = C20a_d1-iv_C20a*sqrt(((expiry_C20a-dates(1))/365));

P16m_d1 = (log(SPX(1)/K_P16m)+(RF3m(1)+(iv_P16m^2)/2)*((expiry_P16m-dates(1))/365))/(iv_P16m*sqrt(((expiry_P16m-dates(1))/365)));
P16m_d2 = P16m_d1-iv_P16m*sqrt(((expiry_P16m-dates(1))/365));

delta_C16m = normcdf(C16m_d1);
delta_C20a = normcdf(C20a_d1);
delta_P16m = normcdf(P16m_d1)-1;

vega_C16m = SPX(1)*sqrt(((expiry_C16m-dates(1))/365))*normpdf(C16m_d1);
vega_C20a = SPX(1)*sqrt(((expiry_C20a-dates(1))/365))*normpdf(C20a_d1);
vega_P16m = SPX(1)*sqrt(((expiry_P16m-dates(1))/365))*normpdf(P16m_d1);

rho_C16m=K_C16m*((expiry_C16m-dates(1))/365)*(exp(-RF3m(1)*((expiry_C16m-dates(1))/365)))*normcdf(C16m_d2); 
rho_C20a=K_C20a*((expiry_C20a-dates(1))/365)*(exp(-RF3m(1)*((expiry_C20a-dates(1))/365)))*normcdf(C20a_d2); 
rho_P16m=-K_P16m*((expiry_P16m-dates(1))/365)*(exp(-RF3m(1)*((expiry_P16m-dates(1))/365)))*normcdf(-P16m_d2); 

G=[delta_C16m vega_C16m rho_C16m ; delta_C20a vega_C20a rho_C20a ; delta_P16m vega_P16m rho_P16m];

% for i=1:length(delta_C16m')-1
%     value_cng(i) = delta_C16m(i)*(SPX(i)-SPX(i+1))+vega_C16m(i)*(VIX(i)-VIX(i+1))+rho_C16m(i)*(RF3m(i)-RF3m(i+1));
% end

Price_C16m = SPX(1)*normcdf(C16m_d1)-K_C16m*exp(-RF3m(1)*((expiry_C16m-dates(1))/365))*normcdf(C16m_d2);
Price_C20a = SPX(1)*normcdf(C20a_d1)-K_C20a*exp(-RF3m(1)*((expiry_C20a-dates(1))/365))*normcdf(C20a_d2);
Price_P16m = K_P16m*exp(-RF3m(1)*((expiry_P16m-dates(1))/365))*normcdf(-P16m_d2)-SPX(1)*normcdf(-P16m_d1);

Price_P16m_v2 = Price_C16m + K_P16m*exp(-RF3m(1)*((expiry_C16m-dates(1))/365))-SPX(1);

h=[10^4; 2*10^4 ; 10^4];

SPX_dif = SPX(1:end-1)-SPX(2:end);
VIX_dif = VIX(1:end-1)-VIX(2:end);
RF3m_dif = RF3m(1:end-1)-RF3m(2:end);

A=[SPX_dif,VIX_dif,RF3m_dif];
C=cov(A);

V = [Price_C16m Price_C20a Price_P16m]*h; % Portfolio Value

sigma_sq = (1/V^2)*h'*G*C*G'*h;

VaR = V*norminv(0.99)*sqrt(sigma_sq);

VaR_procent=VaR/V;

% 3b

asset_contr = (norminv(0.99)*sqrt(1)*G*C*G'*h)/sqrt((V^2*sigma_sq));
fact_contr = (norminv(0.99)*C*G'*h)/sqrt((V^2*sigma_sq));

% FRÅGA PONTUS! 

% 4
% x0=[0.94];
% A=[1];
% b=[1];
% Aeq=[0];
% beq=[0];
% lb=0;
% ub=1;
% 
% fun=@(x0)-sum((-log(EWMA_serie(x0, avg_ret)')-avg_ret(2:end).^2./EWMA_serie(x0, avg_ret)'));
% [lambda, fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)

fun = @(beta, zeta, y)-sum(ln((1+zeta*y/beta)/beta)^((-1/zeta)-1)); 








