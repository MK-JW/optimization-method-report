function [x_optimal,f_optimal,k]=Conjugate_gradient_DY(f_test,g_test,x_initial,tolerance)
%==========================================================================
%函数调用格式：
%[x_optimal,f_optimal,k]=Conjugate_gradient_DY(f_test,g_test,x_initial,tolerance)
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
%k：完成共轭梯度法—DY所需的迭代次数
%==========================================================================
%==========================================================================
%主程序及说明
%--------------------------------------------------------------------------
%在调用本程序的上级程序中定义f_test、g_test的表达式
%x_current：共轭梯度法迭代过程中变量的当前点
%g_current：f_test在x_current处的梯度
%d_current：在x_current处的共轭梯度方向
%x_next:共轭梯度法搜索到的下一点
%f_next:f_next对应x_next处的函数值
%收敛准则：变量变化量小于tol，（判断梯度值范数小于tol，停止）
%搜索步长采用Wolfe_search非精确方法，其中步长最小值可精确到1e-12
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
d_current=-g_current;
if(norm(g_current)<=tolerance)
  x_optimal=x_current;
  f_optimal=f_current;
else
  [x_next,f_next]=Wolfe__Search(f_test,g_test,x_current,d_current);
  fprintf(fileID2,'%4.0f%15.4f%15.4f%15.4f\r\n',k+1,x_next,f_next);
  while(norm(x_next-x_current)>tolerance)
       k=k+1;
       g_previous=g_current;
       d_previous=d_current;
       x_current=x_next;
       f_current=f_test(x_current);
       g_current=g_test(x_current);
     if(norm(g_current)>tolerance)
         if(mod(k,n)==0)
            d_current=-g_current;
         else
           beta_current=((g_current')*g_current)/((d_previous')*(g_current-g_previous));
           d_current=-g_current+beta_current*d_previous;
         end
        [x_next,f_next]=Wolfe__Search(f_test,g_test,x_current,d_current);
        fprintf(fileID2,'%4.0f%15.4f%15.4f%15.4f\r\n',k+1,x_next,f_next);
     end
     x_optimal=x_next;
     f_optimal=f_next;
  end
end
end
%==========================================================================

