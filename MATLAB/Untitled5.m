
gam=[10^(2.2),10];
k=2;
R_min = 0; R_max = 10;
%lam=[rand() rand()];
lam=[1 1];
lam_new=lam;
P=[10 0; 0 10];

err = 0.001;

while R_max - R_min > err
    R_max
    R_min
    R_m = mean([R_min, R_max]);%%(1)
    lam=lam_new;%%(2)
    [R_t, t] = t_opt(lam, gam);%%(3)
    g_lam = - sum(lam.* R_t) + sum(lam) * R_m;%%%%(4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if g_lam > 0 %infeasible%%(5)
        R_max = R_m;
        
        
    elseif g_lam <=0 %feasible
        
        %%%%%%%%%%% min of f_lam is same max of -g_lam
        lam_new = lam;
        f_lam_new = sum(lam_new .* R_t) - sum(lam_new) * R_m;
        f_lam = 1000;
        n = 1;
        g = R_t - R_m;
lam=10;

        while abs(sqrt(sum((lam -lam_new).^2))) > 0.0001% sqrt(g*P*g')> 0.1%(5)update lambda
            lam = lam_new
            g = R_t - R_m; %%subgradient of -g(lambda)
            g_=g/sqrt(g*P*g');
            step=g_*P/(n+1);%%step size
            

            f_lam = f_lam_new;
            
            lam_new = lam - step; %x(k+1)=x(k)-a*g(k)
            P_new=(P-P*g_'*g_*P*2/(n+1))*(n^2)/(n^2-1);
            [R_t, t] = t_opt(lam_new, gam);
            f_lam_new = sum(lam_new .* R_t) - sum(lam_new) * R_m%%3
            
        end
        %%%%%%%%%%%%%%%%
        
        g_lam = - sum(lam_new .* R_t) + sum(lam_new) * R_m;
        
        
        
        if g_lam <= 0%% if g(lam*)<0 ,then feasible thus R_m<R_m*
            R_min = R_m;%(6)
            
        end
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%55
end


disp([' [t0, t1, t2] = ', num2str([t])])
disp([' [R1, R2] = ', num2str([R_t])])
disp(['min(R1, R2) = ', num2str(min(R_t))])