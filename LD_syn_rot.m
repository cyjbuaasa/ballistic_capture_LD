function dzdt = LD_syn_rot(t,Z,diffstr,LD_str,eps,alp,keta,~)
X = Z(1:6);
% dzdt = zeros(length(x)+1,1);
% X_RSM=trans_REM_RSE(X,t);
% X_RSM = X;X_RSM0 = X0;
% 
mu_EMS = 3.040423577705280e-06;
DU = 384405;
yy = diff_sbcm(t,X);
mu=0.012150668300000;
if strcmp(diffstr,'diff_sbcm')
    yy = diff_sbcm(t,X);
elseif strcmp(diffstr,'EM')
    yy = orbitdym(t,X,mu);
else
    yy=zeros(6,1);
end

u2=mu;u1=1-mu;
x = X(1:3);
v = [X(4:6)-[X(2);-X(1);0]];
x2 = X(1:3)-[1-mu;0;0];
v2 = [X(4:6)-[X(2);-X(1)+1-mu;0]];
%任意点位置与受力
r=[((x(1)+u2)^2+x(2)^2+x(3)^2)^0.5;((x(1)-u1)^2+x(2)^2+x(3)^2)^0.5];
a2=[-u2*(x(1)-u1)/(r(2)^3);-u2*x(2)/(r(2)^3);-u2*x(3)/(r(2)^3)];
if strcmp(LD_str,'E2')
    E2=E_secondary(X,mu);
    rm_hill=66180/DU;
    r2n = norm(X(1:3)-[1-mu;0;0]);
    f1 = exp(-alp*E2);
    % f = exp(-1.*E2)*exp(-0.5*(r2n-1.*rm_hill)/1);
    % f = exp(-1.*E2)*rm_hill/r2n;
elseif strcmp(LD_str,'H2')
    H2 = cross(x2,v2);
    f1= exp(-alp*norm(H2));
elseif strcmp(LD_str,'r2')
    f1= exp(-alp*norm(x2));
elseif strcmp(LD_str,'v2')
    f1= norm(v2);
elseif strcmp(LD_str,'onlyv')
    f1= norm(X(4:6));
else
    f1= norm(X(4:6));
end
E2=E_secondary(X,mu);
f2= norm(a2)^(1/1)*exp(-keta*E2);

dzdt = [yy;f1;f2];

end