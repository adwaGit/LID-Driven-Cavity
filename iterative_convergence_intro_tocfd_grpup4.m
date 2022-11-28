% continuity residual 
for i = 2 : imax-1 
for j = 2 : jmax-1

res1(r,r) = (rho/(2*(del_x)))*(oldu(i+1,j)-oldu(i-1,j))  +  (rho/(2*(del_y)))*(oldv(i,j+1)-oldv(i,j-1)) - S(i,j) - f_mass;

% x momentum residual 
res2(r,r) = (rho*oldu(i,j))*((oldu(i+1,j) - oldu(i-1,j))/(2*del_x)) + (rho*oldv(i,j))*((oldu(i,j+1)-oldu(i,j-1))/(2*del_y)) + (oldp(i+1,j)-oldp(i-1,j)/(2*del_x)) -(mu* (oldu(i+1,j)-2*oldu(i,j) + oldu(i-1,j))/del_x^2) - (mu* (oldu(i,j+1)-2*oldu(i,j) + oldu(i,j-1))/del_y^2) - f_xmtm;

% y momentum residual
res3(r,r) = (rho*oldv(i,j))*((oldv(i,j+1) - oldv(i,j-1))/(2*del_x)) + (rho*oldu(i,j))*((oldv(i+1,j)-oldv(i-1,j))/(2*del_y)) + (oldp(i+1,j)-oldp(i-1,j)/(2*del_x)) -(mu* (oldv(i,j+1)-2*oldv(i,j) + oldv(i,j-1))/del_x^2) - (mu* (oldv(i+1,j)-2*oldv(i,j) + oldv(i-1,j))/del_y^2) - f_xmtm;
end 
end 


if (q == ninit)
    
res1_init = res1;    
res2_init = res2;
res3_init = res3;

S1_init = sum(res1_init,'all');
S2_init = sum(res2_init,'all');
S3_init = sum(res3_init,'all');

N = imax*jmax;           % no. of nodes 

L1norm_res1_init = S1_init/ N;
L1norm_res2_init = S2_init/ N;
L1norm_res3_init = S3_init/ N;

% checking convergence  
if rem(q,1000) == 0
 
 S1 = sum(res1,'all');
 S2 = sum(res2,'all');
 S3 = sum(res3,'all');

 N = imax*jmax;  % no. of nodes 

 L1norm_res1 = S1/ N;
 L1norm_res2 = S2/ N;
 L1norm_res3 = S3/ N;
% potting residuals v/s iterations 
residual_plot1(1,q/1000) = L1norm_res1/L1norm_res1_init; 
residual_plot2(1,q/1000) = L1norm_res2/L1norm_res2_init;
residual_plot3(1,q/1000) = L1norm_res3/L1norm_res3_init;


% norm of residual at current time step/ norm of residual at 1st time
 % step. <=10^-8 
 if ( (L1norm_res1/L1norm_res1_init<=10^-8) && (L1norm_res2/L1norm_res2_init<=10^-8) && (L1norm_res3/L1norm_res3_init<=10^-8))
       fprintf(' The solution has converged');
   else 
       fprintf(' The solution is not converged');

 end 
end                                                     
end
r = r +1; 