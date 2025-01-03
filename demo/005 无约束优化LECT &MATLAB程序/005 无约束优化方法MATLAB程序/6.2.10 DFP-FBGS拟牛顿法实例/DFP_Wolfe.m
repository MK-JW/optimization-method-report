function [x_optimal,f_optimal,k]=DFP_Wolfe(f_test,g_test,x_initial,tolerance)
%==========================================================================
%函数调用格式：
%[x_optimal,f_optimal,k]=DFP_Wolfe(f_test,g_test,x_initial,tolerance)
%--------------------------------------------------------------------------
%输入参数说明
%--------------------------------------------------------------------------
%f_test：目标函数
%g_test：f对变量x的梯度
%x_initial：指定的初始点
%tolerance:指定误差
%--------------------------------------------------------------------------
%输出参数
%--------------------------------------------------------------------------
%x_optimal：最优点
%f_optimal：f_test对应x_optimal的函数值
%k：完成DFP法所需的迭代次数
%==========================================================================
%==========================================================================
%主程序及说明
%--------------------------------------------------------------------------
%在调用本程序的上级程序中定义f_test、g_test的表达式
%x_current：DFP法迭代过程中变量的当前点
%g_current：f_test在x_current处的梯度
%d_current：在x_current处的拟牛顿方向
%x_next:DFP法搜索到的下一点
%f_next:f_next对应x_next处的函数值
%收敛准则：变量变化量小于tol，（判断梯度值范数小于tol，停止）
%搜索步长采用Wolfe__search非精确方法，其中步长最小值可精确到1e-12
%--------
%%[x_next,f_next]=Wolfe__Search(@f_test,@g_test,x_current,d_current)
%----
%--------------------------------------------------------------------------
k=1;
%------------open a file to save Data-----------------------
fileID2=fopen('testdata.txt','w');
fprintf(fileID2,'%5s %20s %25s\r\n','k','x^(k)','f^(k)');
%----------------------------------------------------------
n=length(x_initial);
x_current=x_initial;
f_current=f_test(x_current);
fprintf(fileID2,'%4.0f%15.4f%15.4f%15.4f\r\n',k,x_current,f_current);
g_current=g_test(x_current);
D_current=eye(n,n);
d_current=-D_current*g_current;
if(norm(g_current)<=tolerance)
  x_optimal=x_current;
  f_optimal=f_current;
else
  [x_next,f_next]=Wolfe__Search(f_test,g_test,x_current,d_current);
  fprintf(fileID2,'%4.0f%15.4f%15.4f%15.4f\r\n',k+1,x_next,f_next);
  while(norm(x_next-x_current)>tolerance)
       k=k+1;
       x_previous=x_current;
       g_previous=g_current;
       D_previous=D_current;
       x_current=x_next;
       f_current=f_test(x_current);
       g_current=g_test(x_current);
       delta_current=x_current-x_previous;
       gama_current=g_current-g_previous;
     if(norm(g_current)>1e-3*tolerance)
         if((delta_current')*gama_current<=0)
            D_current=eye(n,n);
         else
             Dg=D_previous*gama_current;
           D_current=D_previous+(delta_current*(delta_current'))/((delta_current')*gama_current)-(Dg*(Dg'))/((gama_current')*Dg);
         end
         d_current=-D_current*g_current;
        [x_next,f_next]=Wolfe__Search(f_test,g_test,x_current,d_current);
        fprintf(fileID2,'%4.0f%15.4f%15.4f%15.4f\r\n',k+1,x_next,f_next);
     end
     x_optimal=x_next;
     f_optimal=f_next;
  end
end
end
%==========================================================================

