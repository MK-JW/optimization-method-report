function [x_next,f_next]=armijo_mk(f_test,g_test,x_current,d_current)
beta=0.5; sigma=0.2;
m=0; mmax=20;
xk=x_current;
dk=d_current;
mk=m;
if(norm(dk)>=1.e-6)
    while (m<=mmax)
      if(f_test(xk+beta^m*dk)<=f_test(xk)+sigma*beta^m*g_test(xk)'*dk)
        mk=m; break;
      end
    m=m+1;
   end
end
alpha=beta^mk;
x_next=xk+alpha*dk;
f_next=f_test(x_next);
end