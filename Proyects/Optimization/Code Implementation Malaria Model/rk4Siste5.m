function [x,y,z,a,b,t] =  rk4Siste5(f,g,j,p,k,x0,y0,z0,a0,b0,to,tf,h)
t = to:h:tf ;
n = length(t);
y = [y0] ;
x = [x0] ;
z = [z0] ;
a = [a0] ;
b = [b0] ;
for i = 1:n
    k1 = h*f(x(i),y(i),z(i),a(i),b(i),t(i));
    l1 = h*g(x(i),y(i),z(i),a(i),b(i),t(i));
    m1 = h*j(x(i),y(i),z(i),a(i),b(i),t(i));
    n1 = h*p(x(i),y(i),z(i),a(i),b(i),t(i));
    u1 = h*k(x(i),y(i),z(i),a(i),b(i),t(i));
    
    k2 = h*f(x(i)+0.5*k1,y(i)+0.5*l1, z(i)+0.5*m1 ,a(i)+0.5*n1,b(i)+0.5*u1,t(i) + h*0.5);
    l2 = h*g(x(i)+0.5*k1,y(i)+0.5*l1, z(i)+0.5*m1 ,a(i)+0.5*n1,b(i)+0.5*u1,t(i) + h*0.5);
    m2 = h*j(x(i)+0.5*k1,y(i)+0.5*l1, z(i)+0.5*m1 ,a(i)+0.5*n1,b(i)+0.5*u1,t(i) + h*0.5);
    n2 = h*p(x(i)+0.5*k1,y(i)+0.5*l1, z(i)+0.5*m1 ,a(i)+0.5*n1,b(i)+0.5*u1,t(i) + h*0.5);
    u2 = h*k(x(i)+0.5*k1,y(i)+0.5*l1, z(i)+0.5*m1 ,a(i)+0.5*n1,b(i)+0.5*u1,t(i) + h*0.5);
    
    k3 = h*f(x(i)+0.5*k2,y(i)+0.5*l2, z(i)+0.5*m2 ,a(i)+0.5*n2,b(i)+0.5*u2,t(i) + h*0.5);
    l3 = h*g(x(i)+0.5*k2,y(i)+0.5*l2, z(i)+0.5*m2 ,a(i)+0.5*n2,b(i)+0.5*u2,t(i) + h*0.5);
    m3 = h*j(x(i)+0.5*k2,y(i)+0.5*l2, z(i)+0.5*m2 ,a(i)+0.5*n2,b(i)+0.5*u2,t(i) + h*0.5);
    n3 = h*p(x(i)+0.5*k2,y(i)+0.5*l2, z(i)+0.5*m2 ,a(i)+0.5*n2,b(i)+0.5*u2,t(i) + h*0.5);
    u3 = h*k(x(i)+0.5*k2,y(i)+0.5*l2, z(i)+0.5*m2 ,a(i)+0.5*n2,b(i)+0.5*u2,t(i) + h*0.5);
    
    k4 = h*f(x(i)+k3,y(i)+ k3, z(i) + k3,a(i)+ k3,b(i)+k3, t(i)+h);
    l4 = h*g(x(i)+l3,y(i)+ l3,z(i) + l3,a(i) + l3,b(i) + l3, t(i)+h);
    m4 = h*j(x(i)+m3,y(i)+ m3,z(i) + m3,a(i) + m3,b(i) + m3, t(i)+h);
    n4 = h*p(x(i)+n3,y(i)+ n3,z(i) + n3,a(i) + n3,b(i) + n3, t(i)+h);
    u4 = h*k(x(i)+u3,y(i)+ u3,z(i) + u3,a(i) + u3,b(i) + u3, t(i)+h);
    
    x(i+1) = x(i) + (1/6)*(k1+2*k2+2*k3+k4);
    y(i+1) = y(i) + (1/6)*(l1+2*l2+2*l3+l4);
    z(i+1) = z(i) + (1/6)*(m1+2*m2+2*m3+m4);
    a(i+1) = a(i) + (1/6)*(n1+2*n2+2*n3+n4);
    b(i+1) = b(i) + (1/6)*(u1+2*u2+2*u3+u4);
end

end

