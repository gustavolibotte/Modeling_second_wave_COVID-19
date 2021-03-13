% Model equations for the scaled SIDR model

function dxdt = sirdODE(t,x,params)

%x(1) = susceptible s
%x(2) = infected i
%x(3) = dead d
%x(4) = recovered r

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

if t<=params(12)
beta0 = params(1);
gamma = params(2);
alpha = params(4);
beta1 = params(5);
tau = params(6);
beta=beta0*exp(-(beta1*(t-tau))^2);
else
beta0 = params(7);
gamma = params(8);
alpha = params(9);
beta1 = params(10);
tau = params(11);   
beta=beta0*exp(-(beta1*(t-tau))^2);
end    

% mathematical model
dxdt = zeros(4,1); %column vector for the state variables

dxdt(1) = -beta*x(1)*x(2);
dxdt(2) = beta*x(1)*x(2)-alpha*x(2)-gamma*x(2);
dxdt(3) = gamma*x(2);
dxdt(4) = alpha*x(2);