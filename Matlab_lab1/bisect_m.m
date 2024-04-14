function out=bisect_m(f, a, b, tol,n)

if f(a)*f(b)>0.0
error('function has same sign at end points')
end

a(1)=a;
b(1)=b;

i=1;

while i<=n
     c(i)=(a(i)+b(i))/2;
%      if f(c(i))==0 || (b(i)-a(i))/2<tol
     if  abs(f(c(i)))<tol      
         break;
     end
     if f(a(i))*f(c(i))>0
        a(i+1)=c(i);
        b(i+1)=b(i);
     else
        a(i+1)=a(i);
        b(i+1)=c(i);
     end
     i=i+1;
end

if i>n
    k=1:n;
else
    k=1:i;
end

out=[k', a(k)', b(k)', c(k)', f(c(k))'];

x=[a:0.1:b]
y1=x.^3-x.^2-x+1;
y2=0;
figure
plot(x,y1,'r-',x,y2,'b-*','MarkerSize',2)
xlabel ( 'Axonas x','FontName','Arial','FontSize',10); 
ylabel ( 'Axonas y','FontName','Arial','FontSize',10);
title( 'Grafikes parastaseis y1, y2','FontSize',10 );

