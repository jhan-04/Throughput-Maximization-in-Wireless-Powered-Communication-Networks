clear all
clc
m=0;
mt1=0;
mt2=0;

gam1=10^(2.2);
gam2=10;
%R 은 2차원에서 R다루는 r은 직선에서 계산
n=1000;

%[t1,t2] = meshgrid(0.001:0.001:1);
t1=1/n:1/n:1;%x축 열이 증가할때마다 증가
t2=t1';
t0=(1-t1-t2);
R1=(t1.*log2(1+gam1*(t0./t1))).*(t1+t2<1)+0;
R2=(t2.*log2(1+gam2*(t0./t2))).*(t1+t2<1)+0;
R=R1+R2;
R_=min(R1,R2);
hold on
contour(t1,t2,R_,'ShowText','on')


R_max=max(max(R_))
[I_t2,I_t1]=find(R_==R_max);

t1_m=t1(I_t1)
t2_m=t2(I_t2)
R1_m=R1(I_t2,I_t1)
R2_m=R2(I_t2,I_t1)
plot(t1_m, t2_m,'o')

xlabel('t1(allocated to U1)')
ylabel('t2(allocated to U2)')



%%%%%%%%%%%%%%%%%%%%


gam=[10^(2.2),10];%gam1 gam2
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

plot(t(2), t(3),'*')
hold off
