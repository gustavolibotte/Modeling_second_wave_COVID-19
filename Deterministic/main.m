
% Simulation of the SIDR model with parameters obtained from the 
% solution of the proposed inverse problem using epidemiological data 
% from Germany in the period between February 2 and November 25, 2020.

% Nomenclature

% s = susceptible
% i = infected
% d = dead
% r = recovered

% params(1)= transmission rate - beta_0(first wave)
% params(2)= recovery rate (first wave)
% params(3)= number used to normalize all populations
% params(4)= death rate (first wave)
% params(5)= transmission rate - beta_1(first wave)
% params(6)= transmission rate - beta_2(first wave)

% params(7)= transmission rate - beta_0(second wave)
% params(8)= recovery rate (second wave)
% params(9)= death rate(second wave)
% params(10)= transmission rate - beta_1(second wave)
% params(11)= transmission rate - beta_2(second wave)

% params(12)= parameter tau (interface between the end of the first 
% wave and the start of the second wave)

% A complete description of model and parameters can be found in 
% "Mathematical modelling of the second wave of COVID-19
% infections using deterministic and stochastic SIDR models" 
% by Fran S. Lobato*, Gustavo B. Libotte and Gustavo M. Platt
% *fslobato@ufu.br

function main

clc 
close all
clear all

% Input parameters (see nomenclature section)
params=[0.313492 0.002907 678899.137368 0.045564 0.021907 0.101406 ...
    0.095174 0.000581 0.008549 0.000132 0.000338 166.180157];

% Date
%--------------------------------------------------------------------------
% Infected
Data=[13 13 9 7 7 3 2 2 2 2 3 11 32 58 63 114 149 187 246 528 652 782 1022 1204 1545 1938 2714 3621 4544 5754 7188 9274 12194 15161 19600 22071 24513 28480 29542 33570 37998 43862 48781 52683 52740 54933 58350 61247 65309 68248 69839 72865 69566 64647 63221 65522 65181 64532 62578 60515 58349 56646 53931 53786 53100 50703 48167 45933 44254 42439 40836 39794 38132 36198 34672 32886 30441 29155 28198 26459 24914 23191 22138 21378 20475 19910 19298 18233 17537 16747 15998 15617 15202 14566 13934 13361 12712 12361 11720 11657 11161 10790 10562 10682 10325 9794 9689 9247 9017 8426 8387 8151 8027 7993 7822 7485 6966 6744 6788 6656 6601 6559 6372 6977 7080 7300 7555 7713 7850 8092 7951 7973 8273 8163 8135 8251 7680 7463 7353 6927 7037 6772 6765 6552 6950 6473 6458 6178 6216 6197 6122 6104 6279 5985 5910 5882 6514 6610 6688 6955 6359 6530 6938 6774 6744 7599 8432 9141 8251 8636 8388 9148 10159 9758 10861 10236 10621 11335 11362 11674 12188 11935 12288 12807 14490 15900 16486 17160 17893 15576 15557 17181 18627 15415 15711 15693 15978 16089 15521 14820 14898 14815 16280 15447 16115 17220 14947 15039 15388 17002 15819 16170 16235 16299 18316 18285 18780 19342 20007 19770 19785 20197 22326 24676 25993 26004 26673 26710 27340 28044 29267 29631 30069 31341 31884 33761 36347 38991 40262 41889 44473 46839 50071 54406 59356 61880 65215 69032 72643 79256 87730 96960 103597 110399 118576 126456 136362 148718 161497 175507 177824 186652 194748 204931 216426 227398 235046 241602 247815 252580 259294 269562 278619 283180 287852 289435 290384 294541 301721 308483 310932 313968 313265 311341 311850 313440 315053];
%--------------------------------------------------------------------------
% Dead
Data_m=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2 3 6 8 9 13 17 26 28 44 68 84 94 123 159 206 267 351 433 541 645 775 931 1107 1275 1444 1584 1810 2016 2349 2607 2736 2871 3022 3194 3495 3804 4052 4352 4538 4642 4862 5086 5315 5575 5760 5877 5976 6126 6314 6467 6623 6736 6812 6866 6993 6993 7275 7392 7510 7549 7569 7661 7738 7861 7928 8001 8027 8049 8123 8193 8270 8309 8352 8366 8371 8428 8498 8533 8570 8594 8600 8605 8618 8674 8699 8736 8763 8769 8776 8783 8831 8844 8851 8863 8867 8870 8885 8910 8927 8946 8960 8961 8962 8969 8986 9003 9012 9026 9026 9029 9041 9052 9061 9064 9073 9081 9086 9092 9103 9115 9125 9130 9134 9134 9139 9144 9148 9157 9160 9162 9163 9173 9180 9182 9187 9201 9202 9203 9205 9207 9212 9221 9224 9226 9226 9232 9232 9245 9252 9254 9260 9260 9265 9268 9276 9281 9289 9290 9290 9296 9305 9314 9324 9328 9331 9332 9336 9345 9352 9359 9360 9363 9364 9371 9381 9393 9399 9401 9401 9401 9405 9409 9410 9419 9423 9427 9428 9436 9445 9449 9457 9464 9466 9470 9481 9491 9508 9519 9530 9532 9534 9545 9556 9571 9586 9596 9597 9602 9616 9635 9652 9667 9687 9691 9702 9721 9740 9771 9810 9836 9853 9866 9899 9955 9999 10044 10090 10111 10138 10182 10263 10359 10435 10523 10583 10622 10734 10883 11028 11190 11364 11435 11505 11657 11860 12082 12276 12503 12619 12692 12891 13248 13492 13788 14076 14239 14343 14583 14965 15381 15767 16172];
%--------------------------------------------------------------------------

npontos=length(Data);
data = cumsum(Data(1:npontos))'; % Infected
data_m = Data_m(1:npontos)'; % Dead

tspan = 1 : length(data); 
tspan0 = 1 : length(Data);
x0fcn = @(params) [1-data(1)/params(3); data(1)/params(3); 0; 0]; % S I D R

% Rotine to integrate the SIDR model
[test,x] = ode45(@sirdODE,tspan0,x0fcn(params),[],params);

y1=x(:,1)*params(3); % s
y2=x(:,2)*params(3); % i
y3=x(:,3)*params(3); % d
y4=x(:,4)*params(3); % r

nn=length(test);
icont=0;
for i=1:nn
 if test(i)<=params(12)
  beta0 = params(1);
  gamma = params(2);
  alpha = params(4);
  beta1 = params(5);
  tau = params(6);
  beta(i)=beta0*exp(-(beta1*(test(i)-tau))^2);
  Ro(i)=(beta(i)/(gamma+alpha));
  Rt(i)=(beta(i)/(gamma+alpha))*x(i,1);
  icont=icont+1;
 else
  beta0 = params(7);
  gamma = params(8);
  alpha = params(9);
  beta1 = params(10);
  tau = params(11);   
  beta(i)=beta0*exp(-(beta1*(test(i)-tau))^2);
  Ro(i)=(beta(i)/(gamma+alpha));
  Rt(i)=(beta(i)/(gamma+alpha))*x(i,1);
end    
end

% Plot results

% Infecteds (cumulative)
figure(1) 
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(tspan,data,'ko','LineWidth',2);
plot(test,cumsum(y2),'r','LineWidth',2);
hold off
ylabel('I (Cumulative Number of Inf.)');
xlabel('Time (Days)');

% Susceptible
figure(2)
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(test,y1,'k','LineWidth',2);
hold off
ylabel('S (Number of Individuals)');
xlabel('Time (Days)');

% Death
figure(3)
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(tspan,data_m,'ko','LineWidth',2);
plot(test,y3,'r','LineWidth',2);
hold off
ylabel('D (Number of Individuals)');
xlabel('Time (Days)');

% Recovered
figure(4)
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(test,y4,'k','LineWidth',2);
hold off
ylabel('R (Number of Individuals)');
xlabel('Time (Days)');

% Infecteds
figure(5)  
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(tspan,Data,'ko','LineWidth',2);
plot(test,y2,'r-','LineWidth',2);
hold off
ylabel('I (Number of Individuals)');
xlabel('Time (Days)');

% Transmission rate model
figure(6)  
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(test,beta,'ko','LineWidth',2);
hold off
ylabel('\beta (1/Day)');
xlabel('Time (Days)');

% Effective Reproduction Number
figure(7)  
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
hold on
plot(test,Rt,'ko','LineWidth',2);
hold off
ylabel('R_t (Effective Reproduction Number)');
xlabel('Time (Days)');

