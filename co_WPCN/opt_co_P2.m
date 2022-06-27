function [t, R] = opt_co_P2(rho)


   % wi=[100 100;20 20;15 15;7 7;1 1];%P0=0P0=5P0=10P0=15
    
   %w=[100/log(1+rho(2)) ,100/log(1+rho(3))]
   w=[0 0]
    R_max=10;
    R_min=0;
    w_new=w;
    
    
    while R_max - R_min >0.001
        R_max;
        R_min;
        
        R_m = mean([R_min, R_max]);%%(1)
        w=w_new;%%(2)
        
        [tou0 ,tou1 ,tou21, tou22, t21, t22 ,lam]=find_t_lam(rho,w_new);
        
        R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
        R1_12=tou1*log2(1+rho(2)*tou0/tou1);
        R1_20=tou21*log2(1+rho(3)*t21/tou21);
        R1_co=min(R1_10+R1_20,R1_12);
        R2_co=tou22*log2(1+rho(3)*t22/tou22);
        R=[R1_co, R2_co];
        g_w = - sum(w.* R) + sum(w) * R_m;
  
        t= [tou0 ,tou1 ,tou21, tou22, t21, t22]

        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if g_w > 0 %infeasible%%(5)
            R_max = R_m;
            
            
        elseif g_w <=0 %feasible
            
            %%%%%%%%%%% min of f_lam is same max of -g_lam
            w_new = w;
            f_w_new = sum(w_new .* R) - sum(w_new) * R_m
            f_w = 1000;
            a = 0;
            g = R - R_m;
            while abs(f_w_new - f_w) > 0.01%(5)update lambda
                w = w_new;
                a = a + 1;
                step=0.1/(a+1)^2;%%step size
                
                f_w = f_w_new;
                g = R - R_m; %%subgradient of -g(lambda)
                w_new = w -step * g; %x(k+1)=x(k)-a*g(k)
                
                [tou0 ,tou1 ,tou21, tou22, t21, t22 ,lam]=find_t_lam(rho,w_new);
               
                R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
                R1_12=tou1*log2(1+rho(2)*tou0/tou1);
                R1_20=tou21*log2(1+rho(3)*t21/tou21);
                R1_co=min(R1_10+R1_20,R1_12);
                R2_co=tou22*log2(1+rho(3)*t22/tou22);
                R=[R1_co, R2_co];
                 
                
                f_w_new = sum(w_new .* R) - sum(w_new) * R_m%%3
                
            end
            
            
            g_w = - sum(w_new .* R) + sum(w_new) * R_m;
            
            
            
            if g_w <= 0%% if g(lam*)<0 ,then feasible thus R_m<R_m*
                R_min = R_m;%(6)
                
            end
            
            
            
        end
   
    end
    
    disp([' [R1, R2] = ', num2str([R])])

end