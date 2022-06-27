
clear all
clc

i=1;

for K=0.3:0.2:0.7
D10=10;
D20=10*K;
D12=10*(1-K);
% K=0.3:0.2:0.9;

P0=10^(20/10);
eff=0.67;%energy conversion efficiency ay receiver of Ui
band=10^6;%
a=2;
p=0.75;%a fixed portion of harvested energy
sigma = -160 + 10*log10(10^6); % noise power,


h10=10^(-3)*D10^(-a);
h20=10^(-3)*D20^(-a);
h12=10^(-3)*D12^(-a);

rho1_10=h10^2*eff*p*P0/10^((sigma)/10);
rho1_12=h10*h12*eff*p*P0/10^((sigma)/10);
rho2=h20^2*eff*p*P0/10^((sigma)/10);
rho=[rho1_10 rho1_12 rho2];



[t, R] = opt_co_P2_w(rho);
% 
% 
% w1=1;
% w2=1;
% w=[w1 w2];
% R_max=10;
% R_min=0;
% err=0.01;
% 
% w_new=w;
% 
% 
% 

% 
% 
% while R_max - R_min > err
%     
%     
%     R_max;
%     R_min;
%     
%     R_m = mean([R_min, R_max]);%%(1)
%     w=w_new;%%(2)
%     
%     [tou0 ,tou1 ,tou21, tou22, t21, t22 ,lam]=find_t_lam(rho,w);
%     
%     R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
%     R1_12=tou1*log2(1+rho(2)*tou0/tou1);
%     R1_20=tou21*log2(1+rho(3)*t21/tou21);
%     R1_co=min(R1_10+R1_20,R1_12);
%     R2_co=tou22*log2(1+rho(3)*t22/tou22);
%     R=[R1_co, R2_co];
%     sum(R);
%     g_w = - sum(w.* R) + sum(w) * R_m;
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if g_w > 0 %infeasible%%(5)
%         R_max = R_m;
%         
%         
%         
%     elseif g_w <=0 %feasible
%         
%         %%%%%%%%%%% min of f_lam is same max of -g_lam
%         w_new = w;
%         f_w_new = sum(w_new .* R) - sum(w_new) * R_m;
%         f_w = 1000;
%         a = 0;
%         g = R - R_m;
%         while abs(f_w_new - f_w) > 0.01%(5)update lambda
%             w = w_new;
%             a = a + 1;
%             step=0.1/(a+1);%%step size
%             
%             f_w = f_w_new;
%             g = R - R_m; %%subgradient of -g(lambda)
%             w_new = w -step * g; %x(k+1)=x(k)-a*g(k)
%             
%             [tou0 ,tou1 ,tou21, tou22, t21, t22 ,lam]=find_t_lam(rho,w_new);
%             t=[tou0 ,tou1 ,tou21, tou22, t21, t22];
%             R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
%             R1_12=tou1*log2(1+rho(2)*tou0/tou1);
%             R1_20=tou21*log2(1+rho(3)*t21/tou21);
%             R1_co=min(R1_10+R1_20,R1_12);
%             R2_co=tou22*log2(1+rho(3)*t22/tou22);
%             R=[R1_co, R2_co];
%             
%             
%             f_w_new = sum(w_new .* R) - sum(w_new) * R_m;%%3
%             
%         end
%         
%         
%         g_w = - sum(w_new .* R) + sum(w_new) * R_m;
%         
%         
%         
%         if g_w <= 0%% if g(lam*)<0 ,then feasible thus R_m<R_m*
%             R_min = R_m;%(6)
%             
%         end
%         
%         
%         
%     end
%     
%     
% end

   co_R_common(i)=min(R);

disp([' [t0, t1, t2] = ', num2str([t])])
disp([' [R1, R2] = ', num2str([R])])
disp(['min(R1, R2) = ', num2str(min(R))])

hold on
subplot(3,1,i)
bar(t)
axis([0,7,0,0.55])
grid on
xlabel(' tou0     tou1      tou21      tou22      t21      t22')
i=i+1;
hold off
end
ylabel(' K= 0.7      0.5      0.3')