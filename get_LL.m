function LL=get_LL(mu,mu_EMS,TU)
T=2*pi*TU; % 
option_LL = optimoptions('fsolve','StepTolerance',1e-12,'FunctionTolerance',1e-12,'OptimalityTolerance',1e-12,'Display','none');
gamma1=fsolve(@(x)(x^5-(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2+2*mu*x-mu)-0,(mu)^(1/3),option_LL);
gamma2=fsolve(@(x)(x^5+(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2-2*mu*x-mu)-0,(mu)^(1/3),option_LL);
gamma3=fsolve(@(x)(x^5+(2+mu)*x^4+(1+2*mu)*x^3-(1-mu)*x^2-2*(1-mu)*x-1+mu)-0,1-7/12*mu,option_LL);
LL.L1=[1-mu-gamma1;0;0];
LL.L2=[1-mu+gamma2;0;0];
LL.L3=[-mu-gamma3;0;0];
LL.L4=[0.5-mu;sqrt(3)/2;0];
LL.L5=[0.5-mu;-sqrt(3)/2;0];
gamma1=fsolve(@(x)(x^5-(3-mu_EMS)*x^4+(3-2*mu_EMS)*x^3-mu_EMS*x^2+2*mu_EMS*x-mu_EMS)-0,(mu_EMS)^(1/3),option_LL);
gamma2=fsolve(@(x)(x^5+(3-mu_EMS)*x^4+(3-2*mu_EMS)*x^3-mu_EMS*x^2-2*mu_EMS*x-mu_EMS)-0,(mu_EMS)^(1/3),option_LL);
gamma3=fsolve(@(x)(x^5+(2+mu_EMS)*x^4+(1+2*mu_EMS)*x^3-(1-mu_EMS)*x^2-2*(1-mu_EMS)*x-1+mu_EMS)-0,1-7/12*mu_EMS,option_LL);
LL.SEL1=[1-mu_EMS-gamma1;0;0];
LL.SEL2=[1-mu_EMS+gamma2;0;0];
LL.SEL3=[-mu_EMS-gamma3;0;0];
LL.SEL4=[0.5-mu_EMS;sqrt(3)/2;0];
LL.SEL5=[0.5-mu_EMS;-sqrt(3)/2;0];


end