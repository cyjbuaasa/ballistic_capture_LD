function dzdt = LD_a2_rot(t,Z,diffstr,keta,~,~,~)
X = Z(1:6);
% dzdt = zeros(length(x)+1,1);
% X_RSM=trans_REM_RSE(X,t);
% X_RSM = X;X_RSM0 = X0;
% 
mu_EMS = 3.040423577705280e-06;mu=0.012150668300000;
if strcmp(diffstr,'diff_sbcm')
    yy = diff_sbcm(t,X);
elseif strcmp(diffstr,'EM')
    yy = orbitdym(t,X,mu);
else
    yy=zeros(6,1);
end

E2=E_secondary(X,mu);
x = X(1:3);
v = [X(4:6)-[X(2);-X(1);0]];

u2=mu;u1=1-mu;
%任意点位置与受力
r=[((x(1)+u2)^2+x(2)^2+x(3)^2)^0.5;((x(1)-u1)^2+x(2)^2+x(3)^2)^0.5];
a2=[-u2*(x(1)-u1)/(r(2)^3);-u2*x(2)/(r(2)^3);-u2*x(3)/(r(2)^3)];

f= norm(a2)^(1/1)*exp(-keta*E2);

dzdt = [yy;f];

end