function expr = Penduloo(t,in2,in3)
%Penduloo
%    EXPR = Penduloo(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    23-Nov-2022 18:41:55

param40 = in3(:,1);
param41 = in3(:,2);
param42 = in3(:,3);
theta = in2(3,:);
x = in2(1,:);
y = in2(2,:);
expr = [y;x-param41.*y+param40.*cos(theta)-x.^3;param42];
