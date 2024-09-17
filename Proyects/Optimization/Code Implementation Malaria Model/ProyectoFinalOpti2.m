%% Cond Iniciales

Sh0 = 100000
Ih0 = 30000
Rh0 = 20000
Sv0 = 50000
Iv0 = 10000

Sh2_0 = 100000
Ih2_0 = 0
Rh2_0 = 0
Sv2_0 = 50000
Iv2_0 = 0

%% Parametros

lambdaH = 200
omega = 0.02
betaH = 0.8
muH = 0.002

xi1 = 0.8
theta1 = [0.6, 0.5]
q1 = [0.1, 0.11]

delta = 0.017
rho = 0
e=1

lambdaV = 1000
betaV = 0.7
muV = 0.12

xi2 = 0.8
theta2 = [0.3, 0.4]
q2 = [0.05,0.04]

lambda = [ 1, 0 ;
           0, 1 ]
% Control Optimo

d1=0.01
d2=0.01
d3=0.01
d4=0.01

c1 = 0.00001
c2 = 0.001
%% Modelo 1

dShdt = @(Sh,Ih,Rh,Sv,Iv,t) lambdaH + omega*Rh - Sh*betaH*e*(Iv/(Sh+Ih+Rh)) - muH*Sh
dIhdt = @(Sh,Ih,Rh,Sv,Iv,t) Sh*betaH*e*(Iv/(Sh+Ih+Rh)) - xi1*theta1(1)*(1-q1(1))*Ih - (delta+rho+muH)*Ih
dRhdt = @(Sh,Ih,Rh,Sv,Iv,t) xi1*theta1(1)*(1-q1(1))*Ih + delta*Ih - (omega + muH)*Rh
dSvdt = @(Sh,Ih,Rh,Sv,Iv,t) lambdaV - Sv*betaV*e*(Ih/(Sh+Ih+Rh)) - xi2*theta2(1)*(1-q2(1))*Sv - muV*Sv
dIvdt = @(Sh,Ih,Rh,Sv,Iv,t) Sv*betaV*e*(Ih/(Sh+Ih+Rh)) - xi2*theta2(1)*(1-q2(1))*Iv - muV*Iv

% Resultados
t0=0;
tf=100;
h=1;

[Sh,Ih,Rh,Sv,Iv,t] = rk4Siste5(dShdt,dIhdt,dRhdt,dSvdt,dIvdt,Sh0,Ih0,Rh0,Sv0,Iv0,t0,tf,h);

t(end+1)= tf+h;

% Visualización de los Datos
hold on 
plot(t, Sh, 'b');
plot(t, Ih, '--r');
plot(t, Rh, ':k');
plot(t, Sv, '-.m');
plot(t, Iv, 'c');

title("Modelo con una población");
xlabel("Time(Days)");
ylabel("Population");
legend('S_h(t)','I_h(t)','R_h(t)','S_v(t)','I_v(t)');
hold off
%% Modelo 2


dShdt1 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) lambdaH + omega*Rh1 - Sh1*(sum(lambda(1,:)))*(2*betaH)*(2*e)*((Iv1+Iv2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - muH*Sh1 %Red flag
dIhdt1 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) Sh1*(sum(lambda(1,:)))*(2*betaH)*(2*e)*((Iv1+Iv2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - xi1*theta1(1)*(1-q1(1))*Ih1 - (delta+rho+muH)*Ih1 %Red flag
dRhdt1 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) xi1*theta1(1)*(1-q1(1))*Ih1 + delta*Ih1 - (omega + muH)*Rh1
dSvdt1 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) lambdaV - Sv1*(sum(lambda(1,:)))*(2*betaV)*(2*e)*((Ih1+Ih2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2))- xi2*theta2(1)*(1-q2(1))*Sv1 - muV*Sv1 %Red flag
dIvdt1 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) Sv1*(sum(lambda(1,:)))*(2*betaV)*(2*e)*((Ih1+Ih2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - xi2*theta2(1)*(1-q2(1))*Iv1 - muV*Iv1

dShdt2 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) lambdaH + omega*Rh2 - Sh2*(sum(lambda(2,:)))*(2*betaH)*(2*e)*((Iv1+Iv2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - muH*Sh2 %Red flag
dIhdt2 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) Sh2*(sum(lambda(2,:)))*(2*betaH)*(2*e)*((Iv1+Iv2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - xi1*theta1(2)*(1-q1(2))*Ih2 - (delta+rho+muH)*Ih2 %Red flag
dRhdt2 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) xi1*theta1(2)*(1-q1(2))*Ih2 + delta*Ih2 - (omega + muH)*Rh2
dSvdt2 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) lambdaV - Sv2*(sum(lambda(2,:)))*(2*betaV)*(2*e)*((Ih1+Ih2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2))- xi2*theta2(2)*(1-q2(2))*Sv2 - muV*Sv2 %Red flag
dIvdt2 = @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t) Sv2*(sum(lambda(2,:)))*(2*betaV)*(2*e)*((Ih1+Ih2)/(Sh1+Ih1+Rh1+Sh2+Ih2+Rh2)) - xi2*theta2(2)*(1-q2(2))*Iv2 - muV*Iv2

% Resultados
t0=0;
tf=1000;
h=1;

[Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t] = rk4Siste10(dShdt1,dIhdt1,dRhdt1,dSvdt1,dIvdt1,dShdt2,dIhdt2,dRhdt2,dSvdt2,dIvdt2,Sh0,Ih0,Rh0,Sv0,Iv0,Sh2_0,Ih2_0,Rh2_0,Sv2_0,Iv2_0,t0,tf,h);

t(end+1)= tf+h;

% Visualización de los Resultados
% Población 1
figure
hold on 
plot(t, Sh1, 'b');
plot(t, Ih1, '--r');
plot(t, Rh1, ':k');
plot(t, Sv1, '-.m');
plot(t, Iv1, 'c');

title("Población #1");
xlabel("Time(Days)");
ylabel("Population");
legend('S_h1(t)','I_h1(t)','R_h1(t)','S_v1(t)','I_v1(t)');
hold off

% Población 2
hold on 
plot(t, Sh2, 'b');
plot(t, Ih2, '--r');
plot(t, Rh2, ':k');
plot(t, Sv2, '-.m');
plot(t, Iv2, 'c');

title("Población #2");
xlabel("Time(Days)");
ylabel("Population");
legend('S_h2(t)','I_h2(t)','R_h2(t)','S_v2(t)','I_v2(t)');
hold off

%% Control Optimo
% Metodo forward-backward sweep

test = -1;


mu = 0.001;

M = tf;

hh = 1/M;
h2 = hh/2;


x1= Sh1;
x2= Ih1;
x3= Rh1;
x4= Sv1;
x5= Iv1;
x6= Sh2;
x7= Ih2;
x8= Rh2;
x9= Sv2;
x10= Iv2;


u1= zeros(1,M+1);
u2= zeros(1,M+1);
u3= zeros(1,M+1);
u4= zeros(1,M+1);
u5= zeros(1,M+1);
u6= zeros(1,M+1);
u7= zeros(1,M+1);
u8= zeros(1,M+1);
u9= zeros(1,M+1);
u10= zeros(1,M+1);

% z1= zeros(1,M+1);
% z2= zeros(1,M+1);
% z3= zeros(1,M+1);
% z4= zeros(1,M+1);
% z5= zeros(1,M+1);
% z6= zeros(1,M+1);
% z7= zeros(1,M+1);
% z8= zeros(1,M+1);
% z9= zeros(1,M+1);
% z10= zeros(1,M+1);


while(test < 0)
    
%     oldu1 = u1;
%     oldu2 = u2;
%     oldu3 = u3;
%     oldu4 = u4;
%     oldu5 = u5;
%     oldu6 = u6;
%     oldu7 = u7;
%     oldu8 = u8;
%     oldu9 = u9;
%     oldu10 = u10;
%     
%     oldx1 = x1;
%     oldx2 = x2;
%     oldx3 = x3;
%     oldx4 = x4;
%     oldx5 = x5;
%     oldx6 = x6;
%     oldx7 = x7;
%     oldx8 = x8;
%     oldx9 = x9;
%     oldx10 = x10;
%     
%     oldlambda1 = lambda1;
%     oldlambda2 = lambda2;
%     oldlambda3 = lambda3;
%     oldlambda4 = lambda4;
%     oldlambda5 = lambda5;
%     oldlambda6 = lambda6;
%     oldlambda7 = lambda7;
%     oldlambda8 = lambda8;
%     oldlambda9 = lambda9;
%     oldlambda10 = lambda10;




    dz1= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) muH*z1 + lambda(1,2)* betaH*e*(Iv2/((Sh2+Ih2+Rh2)))*(z1-z2)  + lambda(1,1)*betaH*e*Iv1*(((Sh1+Ih1+Rh1)-Sh1)/(Sh1+Ih1+Rh1)^2)*(z1-z2) + lambda(1,1)* betaV*e*(Ih1*Sv1/(Sh1+Ih1+Rh1)^2)*(z5-z4) + lambda(2,1)*betaH*e*(Iv1*Sh2/(Sh1+Ih1+Rh1)^2)*(z7-z6) + lambda(1,2)*betaH*e*(Ih1*Sv2/(Sh1+Ih1+Rh1)^2)*(z10-z9); 
    dz2= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -c1 - (xi1*theta1(1)*(1-q1(1))+delta)*z3 + (delta + rho + muH - xi1*theta1(1)*(1-q1(1)))*z2 + lambda(1,1)*betaH*e*(Iv1*Sh1/(Sh1+Ih1+Rh1)^2)*(z2-z1) + lambda(1,1)*betaV*e*Sv1*(((Sh1+Ih1+Rh1)-Ih1)/(Sh1+Ih1+Rh1)^2)*(z4-z5) + lambda(2,1)*betaH*e*(Iv1*Sh2/(Sh1+Ih1+Rh1)^2)*(z7-z6) + lambda(1,2)*betaV*e*Sv2*(((Sh1+Ih1+Rh1)-Ih1)/(Sh1+Ih1+Rh1)^2)*(z9-z10);
    dz3= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -omega*z1 + (omega+muH)*z3 + lambda(1,1)*betaH*e*(Iv1*Sh1/(Sh1+Ih1+Rh1)^2)*(z2-z1) + lambda(1,1)*betaV*e*(Ih1*Sv1/(Sh1+Ih1+Rh1)^2)*(z5-z4) + lambda(2,1)*betaH*e*(Iv1*Sh2/(Sh1+Ih1+Rh1)^2)*(z7-z6) + lambda(1,2)*betaV*e*(Ih1*Sv2/(Sh1+Ih1+Rh1)^2)*(z10-z9);
    dz4= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) (xi2*theta2(1)*(1-q2(1)) + muV)*z4 + (lambda(1,1)*betaV*e*(Ih1/Sh1+Ih1+Rh1) + lambda(2,1)*betaV*e*(Ih2/Sh2+Ih2+Rh2))*(z4-z5);
    dz5= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -c1 + (xi2*theta2(1)*(1-q2(1))+muV)*z5 + (lambda(1,1)*betaH*e*(Sh1/Sh1+Ih1+Rh1)*(z1-z2)) + (lambda(2,1)*betaH*e*(Sh2/Sh1+Ih1+Rh1)*(z6-z7));
    dz6= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -muH*z6 + lambda(1,2)*betaH*e*(Iv2*Sh1/(Sh2+Ih2+Rh2)^2)*(z2-z1) + lambda(2,1)*betaV*e*(Ih2*Sv1/(Sh2+Ih2+Rh2)^2)*(z5-z4) + lambda(2,2)*betaH*e*Iv2*((Sh2+Ih2+Rh2)-Sh2/(Sh2+Ih2+Rh2)^2)*(z6-z7) + lambda(2,1)*betaH*e*(Iv1/(Sh1+Ih1+Rh1))*(z6-z7) + lambda(2,2)*betaV*e*(Ih2*Sv2/(Sh2+Ih2+Rh2)^2)*(z9-z10);
    dz7= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -c1 - (xi1*theta1(2)*(1-q1(2))+delta)*z8 + (xi1*theta1(2)*(1-q1(2))+delta+rho+muH)*z7 + lambda(1,2)*betaH*e*(Iv2*Sh1/(Sh2+Ih2+Rh2)^2)*(z2-z1) + lambda(2,1)*betaV*e*Sv1*((Sh2+Ih2+Rh2)-Ih2/(Sh2+Ih2+Rh2)^2)*(z5-z4) + lambda(2,2)*betaH*e*(Iv2*Sh2/(Sh2+Ih2+Rh2)^2)*(z7-z6) + lambda(2,2)*betaV*e*Sv2*((Sh2+Ih2+Rh2)-Ih2/(Sh2+Ih2+Rh2)^2)*(z9-z10);
    dz8= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -omega*z6 + (omega+muH)*z8 + lambda(1,2)*betaH*e*(Iv2*Sh1/(Sh2+Ih2+Rh2)^2)*(z2-z1) + lambda(2,1)*betaV*e*(Ih2*Sv1/(Sh2+Ih2+Rh2)^2)*(z5-z4) + lambda(2,2)*betaH*e*(Iv2*Sh2/(Sh2+Ih2+Rh2)^2)*(z7-z6) + lambda(2,2)*betaV*e*(Ih2*Sv2/(Sh2+Ih2+Rh2)^2)*(z10-z9);
    dz9= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) (xi2*theta2(2)*(1-q2(2)) + muV)*z9 + (lambda(2,2)*betaV*e*(Ih2/Sh2+Ih2+Rh2) + lambda(1,2)*betaV*e*(Ih2/Sh2+Ih2+Rh2))*(z9-z10);
    dz10= @(Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t) -c2 + (xi2*theta2(2)*(1-q2(2))+muV)*z10 + (lambda(2,2)*betaH*e*(Sh2/Sh2+Ih2+Rh2)*(z6-z7)) + (lambda(1,2)*betaH*e*(Sh1/Sh2+Ih2+Rh2)*(z1-z2));
    
    t0=0;
    tf=10;
    h=1;
    
    [z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,t] = Backward_rk4Sjste10(dz1,dz2,dz3,dz4,dz5,dz6,dz7,dz8,dz9,dz10,Sh1,Ih1,Rh1,Sv1,Iv1,Sh2,Ih2,Rh2,Sv2,Iv2,t0,tf,h);
    %%% ERROR DE LA IMPLEMENTCIÓN RESULTADOS DE LOS Zs NO SON LOS CORRECTOS
    
    
    u11 = min(max(0,(e.*(1-q1(1)).*Ih1.*(z2-z3))/d1),1);    
    
    u21 = min(max(0,(e*(1-q2(1))*(Sv1*z4 + Iv1*z5))/d2),1);
    
    u12 = min(max(0,(e*(1-q1(2))*Ih2*(z7-z8))/d3),1);    
    
    u22 = min(max(0,(e*(1-q2(2))*(Sv2*z9 + Iv2*z10))/d4),1);
    
    
   
    
   J=sum(c1*(Ih1 + Iv1) + c2*(Ih2 + Iv2) + 0.5*(d1*(u11)^2 + d2*(u21)^2 +d3*(u12)^2 +d4*(u22)^2 )); 
    
    
%     temp1 = mu*sum(abs(u1)) - sum(abs(oldu1 - u1));
%     temp2 = mu*sum(abs(u2)) - sum(abs(oldu2 - u2));
%     temp3 = mu*sum(abs(u3)) - sum(abs(oldu3 - u3));
%     temp4 = mu*sum(abs(x1)) - sum(abs(oldx1 - x1));
%     temp5 = mu*sum(abs(x2)) - sum(abs(oldx2 - x2));
%     temp6 = mu*sum(abs(x3)) - sum(abs(oldx3 - x3));
%     temp7 = mu*sum(abs(lambda1)) - sum(abs(oldlambda1 - lambda1));
%     temp8 = mu*sum(abs(lambda2)) - sum(abs(oldlambda2 - lambda2));
%     temp9 = mu*sum(abs(lambda3)) - sum(abs(oldlambda3 - lambda3));
%     
%     test = min(temp1, min(temp2, min(temp3, min(temp4, min(temp5, min(temp6, min(temp7, min(temp8, temp9))))))));
end

