clc
P=10:3:19;
K_e=0.2:0.1:0.4;
for i=1:1:3
    i
%K=K_e(i);  
P0=10^(20/10);
K=K_e(i);
D10=10;
D20=10*K;
D12=10*(1-K);
%P0=10^(20/10);
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


lam_new=[100/log(1+rho(1)) ,100/log(1+rho(3))]
R_max=10;
R_min=0;

[R,t]=t_co_opt(lam,rho);
    disp([' [t0, t1, t2] = ', num2str([t])])
    disp([' [R1, R2] = ', num2str([R])])



err = 0.01;

while R_max - R_min > err
    R_max;
    R_min;
    R_m = mean([R_min, R_max]);%%(1)
    lam=lam_new;%%(2)
    
    
    
     g_lam = - sum(lam_new .* R) + sum(lam_new) * R_m;
    
  %  g_lam = - sum(lam.* R) + sum(lam) * R_m%%%%(4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if g_lam > 0 %infeasible%%(5)
        R_max = R_m
        
        
    elseif g_lam <=0 %feasible
        
        %%%%%%%%%%% min of f_lam is same max of -g_lam
        lam_new = lam;
        f_lam_new = sum(lam_new .* R) - sum(lam_new) * R_m;
        f_lam = 1000;
        a = 0;
        g = R - R_m;
        f_min=10000;
        lam_opt=lam_new;
        while a<=100%abs(f_lam_new - f_lam) > abs(0.01*f_lam)%(5)update lambda
           
            lam = lam_opt;
            a = a + 1;
            step=1/(a+1);%%step size
            
            g = R - R_m; %%subgradient of -g(lambda)
            lam_opt = abs(lam -step * g); %x(k+1)=x(k)-a*g(k)
            
            [R_opt,t_opt]=t_co_opt(lam_opt,rho);

            f_lam_new = sum(lam_opt .* R_opt) - sum(lam_opt) * R_m;%%3
            
            if f_lam_new<f_min
                f_min=f_lam_new;
                R=R_opt;
                t=t_opt;
                lam_new=lam_opt;
                disp([' [t0, t1, t2] = ', num2str([t])])
                disp([' [R1, R2] = ', num2str([R])])
                disp([' [k] = ', num2str([a])])
                disp([' [f_min] = ', num2str([f_min])])
            end
            
            
            
            
            
            
        end
        %%%%%%%%%%%%%%%%
        
        g_lam = - sum(lam_new .* R) + sum(lam_new) * R_m

        if g_lam <= 0%% if g(lam*)<0 ,then feasible thus R_m<R_m*
            R_min = R_m%(6)
            
        end
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%55
end

[R,t]=t_co_opt(lam_new,rho);

co_R_common(i)=min(R);


                disp([' [t0, t1, t2] = ', num2str([t])])
                disp([' [R1, R2] = ', num2str([R])])
disp([' [min R] = ', num2str([min(R)])])


    %%%%%%%%%%%%%%%%%%%%%%%%without user cooperation
    

    
    gam1 = eff*p*h10^2*P0/(10^((sigma)/10));%%==rho
    gam2 = eff*p*h20^2*P0/(10^((sigma)/10));
%    gam=[gam1 gam2];
 gam=[rho1_10 rho2]  ; 
        [tau_opt,R_tau] = opt_P2(gam);
 

R_common(i)=min(R_tau)




end

plot(P,co_R_common,'dr-',P,R_common,'or-');


%plot(P,co_R_common,'dr-');