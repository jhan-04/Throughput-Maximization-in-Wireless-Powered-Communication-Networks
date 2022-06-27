
clear all
clc

K=0.2;
D10=10;
D20=10*K;
D12=10*(1-K);
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
R_opt=0;




w1=1;
w2=1;
%w=[w1 w2];
w=[50/log(1+rho(1)) ,50/log(1+rho(3))]
R_max=10;
R_min=0;
err=0.001;

w_opt=w;


[R_opt,t_opt]=t_co_opt(w,rho);


while R_max - R_min > err
    
    R_max
    R_min
    R_m = mean([R_min, R_max]);%%(1)
    w=w_opt;%%(2)
    R=R_opt;
    
    g_w = - sum(w.* R) + sum(w) * R_m;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if g_w > 0 %infeasible%%(5)
        R_max = R_m;
        
        
        
    elseif g_w <=0 %feasible
        
        %%%%%%%%%%% min of f_lam is same max of -g_lam
        w_new = w;
        f_w_new = sum(w_new .* R) - sum(w_new) * R_m;
        f_w = 1000;
        k = 0;
        g = R - R_m;
        f_min=1000;
        while k<100%abs(f_w_new - f_w) > 0.001%(5)update lambda
            w = w_new;
            k = k + 1;
            step=1/sqrt(k+1);%%step size
            
            f_w = f_w_new;
            g = R - R_m; %%subgradient of -g(lambda)
            w_new = w -step * g; %x(k+1)=x(k)-a*g(k)
            
            
            [R,t]=t_co_opt(w_new,rho);
            

            f_w_new = sum(w_new .* R) - sum(w_new) * R_m;%%3
            if f_w_new<f_min
                f_min=f_w_new;
                t_opt=t;
                R_opt=R;
                w_opt=w_new;
                disp([' [t0, t1, t2] = ', num2str([t])])
                disp([' [R1, R2] = ', num2str([R])])
                disp([' [k] = ', num2str([k])])
                disp([' [f_w_new ] = ', num2str([f_w_new ])])
            end
            
            
            
            
        end
        
        
        g_w = - sum(w_opt .* R_opt) + sum(w_opt) * R_m;
        if g_w <= 0%% if g(lam*)<0 ,then feasible thus R_m<R_m*
            R_min = R_m;%(6)
            
        end
        
        
        
    end
    
    
end

%   co_R_common(i)=min(R);

disp([' [t0, t1, t2] = ', num2str([t])])
disp([' [R1, R2] = ', num2str([R])])
disp(['min(R1, R2) = ', num2str(min(R))])
