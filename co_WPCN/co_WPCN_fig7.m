clear all
clc

eff=0.67;%energy conversion efficiency ay receiver of Ui
band=10^6;%
%a=2.5;
p=0.75;%a fixed portion of harvested energy
sigma = -160 + 10*log10(10^6); % noise power,

for i=1:4
    
    K=0.2:0.2:0.8;
    
    %%%%%%%%%%%%%%%%%%%a=2.5
    a=2;
    P0=10^(20/10);
    
    D10=10;
    D20=10*K(i);
    D12=10*(1-K(i));
    h10=10^(-3)*D10^(-a);
    h20=10^(-3)*D20^(-a);
    h12=10^(-3)*D12^(-a);
    
    rho1_10=h10^2*eff*p*P0/10^((sigma)/10);
    rho1_12=h10*h12*eff*p*P0/10^((sigma)/10);
    rho2=h20^2*eff*p*P0/10^((sigma)/10);
    
    rho=[rho1_10 rho1_12 rho2];
    
%     i
%     
%     [t, R] = opt_co_P2_w(rho);
%     
%     co_R_common_a25(i)=min(R);
%     R_final_co_25(i,:)=R;
%     
%     disp([' [t0, t1, t2] = ', num2str([t])])
%     disp([' [R1, R2] = ', num2str([R])])
%     disp(['min(R1, R2) = ', num2str(min(R))])
%     
%     
%   










    %%%%%%%%%%%%%%%%%%%%%%%%without user cooperation
    %
        gam1 = eff*p*h10^2*P0/(10^((sigma)/10));%%==rho
        gam2 = eff*p*h20^2*P0/(10^((sigma)/10));
        gam=[gam1 gam2];
    
            [tau_opt,R_tau] = opt_P2(gam);
    
    
    R_common_a25(i)=min(R_tau);
    R_final(i,:)=R_tau;
    
        i







    
    %%%%%%%%%%%%%%%%%%%a=3
    
     a=3;
    P0=10^(20/10);
    
    D10=10;
    D20=10*K(i);
    D12=10*(1-K(i));
    h10=10^(-3)*D10^(-a);
    h20=10^(-3)*D20^(-a);
    h12=10^(-3)*D12^(-a);
    
    rho1_10=h10^2*eff*p*P0/10^((sigma)/10);
    rho1_12=h10*h12*eff*p*P0/10^((sigma)/10);
    rho2=h20^2*eff*p*P0/10^((sigma)/10);
    
    rho=[rho1_10 rho1_12 rho2];
    
%     i
%     
%     [t, R] = opt_co_P2_w(rho);
%     
%     co_R_common_a3(i)=min(R);
%     R_final_co_3(i,:)=R;
%     
%     disp([' [t0, t1, t2] = ', num2str([t])])
%     disp([' [R1, R2] = ', num2str([R])])
%     disp(['min(R1, R2) = ', num2str(min(R))])
%     
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%without user cooperation
    %
%         gam1 = eff*p*h10^2*P0/(10^((sigma)/10));%%==rho
%         gam2 = eff*p*h20^2*P0/(10^((sigma)/10));
%         gam=[gam1 gam2];
%     
%             [tau_opt,R_tau] = opt_P2(gam);
%     
%     
%     R_common_a3(i)=min(R_tau);
%     R_final(i,:)=R_tau;
%     
        i
    
    
end



grid on




%plot(K,co_R_common);
%co_R_common

plot(K,co_R_common_a25,'dr-',K,co_R_common_a3,'or-');
plot(K,R_common_a25,'db-',K,R_common_a3,'ob-');


grid on
legend('with User cooperation a=2.5','with user coperation a=3','without User cooperation a=2.5','without user coperation a=3')
xlabel('transmit Power,P0(dBm)')
ylabel('throuput (Mbps)')
axis([0.2,0.7,0,2])

hold off

