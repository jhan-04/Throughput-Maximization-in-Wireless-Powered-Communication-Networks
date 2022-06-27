clear all
clc
m=0;
mt1=0;
mt2=0;



a=[2 2.5 3 3.5 4];

gam1=10^(2.2);
gam2=10.^([10 7 4 1 -2]./10)';
%n=1000;

%[t1,t2] = meshgrid(0.001:0.001:1);
% t1=1/n:1/n:1;%x축 열이 증가할때마다 증가
% t2=t1';
% t0=(1-t1-t2);


for i=1:length(gam2)
gam=[gam1,gam2(i)];
k=2;
R_min = 0; R_max = 10;
%lam=[rand() rand()];
lam=[1 1];
lam_new=lam;

err = 0.001;

while R_max - R_min > err
    R_max;
    R_min;
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
        a = 0;
        g = R_t - R_m;
        while abs(f_lam_new - f_lam) > 0.00001%(5)update lambda
            lam = lam_new;
            a = a + 1;
            step=0.1/(a+1);%%step size
            
            f_lam = f_lam_new;
            g = R_t - R_m; %%subgradient of -g(lambda)
            lam_new = lam -step * g; %x(k+1)=x(k)-a*g(k)
            
            [R_t, t] = t_opt(lam_new, gam);
            f_lam_new = sum(lam_new .* R_t) - sum(lam_new) * R_m;%%3
            
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
P2(i)=t(3)/t(2);
end

semilogy(a,P2,'*-')





syms z
for i=1:length(a)
A=gam1+gam2(i);

eqn=z*log(z)-z+1==A;
optz(i,1)=double(solve(eqn,z));

 end   
opt0=(optz-1)./(A+optz-1);
opt1=gam1./(A+optz-1);
opt2=gam2./(A+optz-1);
R1_op=opt1.*log2(1+gam1.*(opt0./opt1));
R2_op=opt2.*log2(1+gam2.*(opt0./opt2));
P1=opt2./opt1;





a=[2 2.5 3 3.5 4];
semilogy(a,P1,'o-')
hold on
semilogy(a,P2,'*-')

grid on
hold off