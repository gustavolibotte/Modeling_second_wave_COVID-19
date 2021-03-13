
% Simulation of the stochastic SIDR model with parameters obtained 
% from the solution of the proposed inverse problem using epidemiological 
% data from Germany in the period between February 2 and November 25, 2020.

% Nomenclature

% s = susceptible
% i = infecteds
% d = dead
% r = recovered
% sigma = proportionality constant

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

clc
close all
clear all

N=1000;
params=[0.313492 0.002907 678899.137368 0.045564 0.021907 0.101406 ...
    0.095174 0.000581 0.008549 0.000132 0.000338 166.180157];

% Initial parameters obtained from the first wave
lambda0 = 0.094914291760573;
gamma0 = params(8);
alpha0 = params(9);
tau = params(12);

load Population.txt % Simulated profiles for the first wave 
% this profile does not change during the stochastic analysis

[npop,mpop]=size(Population);
tempo=Population(:,1);
[~,icont_tau]=min(abs(tau-tempo)); 
CI=Population(icont_tau,2:end); % Initial Condition
tinicial=tempo(icont_tau); % Initial Time
tfinal=tempo(end); % Final Time
T=linspace(tinicial,tfinal,N); 
dt=T(2)-T(1); 

sigma_par=[0.1 0.2 0.3 0.4]; % Values considered for the 
% proportionality constant
cor_v=['r' 'b' 'm' 'g'];

for ij=1:length(sigma_par)
 sigma=sigma_par(ij);
 mu_lambda=lambda0;
 mu_gamma=gamma0;
 mu_alpha=alpha0;
 sigma_lambda=sigma;
 sigma_gamma=sigma;
 sigma_alpha=sigma;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iv=1:100 % number of runs

S=zeros(1,N);
I=zeros(1,N);
D=zeros(1,N);
R=zeros(1,N);

gamma(1)=gamma0;
alpha(1)=alpha0;
lambda(1)=lambda0;

Rt(1)=lambda(1)/(gamma(1)+alpha(1));

% INITIAL CONDITIONS
S(1)=CI(1);
I(1)=CI(2);
D(1)=CI(3);
R(1)=CI(4);
% EULER SCHEME FOR THE SDEs
for t=1:N-1
Rt(t+1)=S(t)*(lambda(t)/(gamma(t)+alpha(t)));
S(t+1)=S(t)-(lambda(t)*S(t)*I(t)).*dt;
I(t+1)=I(t)+((lambda(t)*S(t)*I(t))-(alpha(t)+gamma(t))*I(t)).*dt;
D(t+1)=D(t)+(gamma(t)*I(t)).*dt;
R(t+1)=R(t)+(alpha(t)*I(t)).*dt;
lambda(t+1)=lambda0*dt+dt*exp(random('Normal',log(mu_lambda-sigma_lambda^2/2),sigma_lambda));
gamma(t+1)=gamma0*dt+dt*random('Normal',mu_gamma,mu_gamma*sigma_lambda);%*dt+gamma0/5*randn*sqrt(dt);
alpha(t+1)=alpha0*dt+dt*random('Normal',mu_alpha,mu_alpha*sigma_lambda);%*dt+alpha0./5*randn*sqrt(dt);
end

figure(1)
subplot(4,2,1)
hold on
plot(T,S*params(3),cor_v(ij));
plot(Population(:,1),Population(:,2)*params(3),'k-');
xlabel('Time (Days)');
ylabel('S');

subplot(4,2,2)
hold on
plot(T,I*params(3),cor_v(ij));
plot(Population(:,1),Population(:,3)*params(3),'k-');
xlabel('Time (Days)');
ylabel('I');

subplot(4,2,3)
hold on
plot(T,D*params(3),cor_v(ij));
plot(Population(:,1),Population(:,4)*params(3),'k-');
xlabel('Time (Days)');
ylabel('D');

subplot(4,2,4)
hold on
plot(T,R*params(3),cor_v(ij));
plot(Population(:,1),Population(:,5)*params(3),'k-');
xlabel('Time (Days)');
ylabel('R');

if ij==1
 subplot(4,2,5)
 hold on
 plot(T,Rt,cor_v(ij));
 xlabel('Time (Days)');
 ylabel('R_t');
 legend('\sigma=0.1');
elseif ij==2
 subplot(4,2,6)
 hold on
 plot(T,Rt,cor_v(ij));
 xlabel('Time (Days)');
 ylabel('R_t');
 legend('\sigma=0.2');
elseif ij==3
 subplot(4,2,7)
 hold on
 plot(T,Rt,cor_v(ij));
 xlabel('Time (Days)');
 ylabel('R_t');
 legend('\sigma=0.3');
else
 subplot(4,2,8)
 hold on
 plot(T,Rt,cor_v(ij));
 xlabel('Time (Days)');
 ylabel('R_t');
 legend('\sigma=0.4');
end
end
end

