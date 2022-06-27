% clc
% clear all
% 

K=0.5;
D10=10;
D20=10*K;
D12=10*(1-K);

for i=5:5
P=0:5:20;
  %  i=5;
   P0=10^(P(i)/10);

    eff=0.5;%energy conversion efficiency ay receiver of Ui
    band=10^6;%
    a=2;
    p=0.5;%a fixed portion of harvested energy
    sigma = -160 + 10*log10(10^6); % noise power,
    h10=10^(-3)*D10^(-a);
    h20=10^(-3)*D20^(-a);
    h12=10^(-3)*D12^(-a);
    
    rho1_10=(h10^2)*eff*p*P0/10^((sigma)/10);
    rho1_12=h10*h12*eff*p*P0/10^((sigma)/10);
    rho2=(h20^2)*eff*p*P0/10^((sigma)/10);
    
    rho=[rho1_10 rho1_12 rho2];
  
i

[t, R] = opt_co_P2_w(rho);

   co_R_common(i)=min(R);
   R_final_co(i,:)=R;
   
disp([' [t0, t1, t2] = ', num2str([t])])
disp([' [R1, R2] = ', num2str([R])])
disp(['min(R1, R2) = ', num2str(min(R))])
    


    %%%%%%%%%%%%%%%%%%%%%%%%without user cooperation
    

    
    gam1 = eff*p*h10^2*P0/(10^((sigma)/10));%%==rho
    gam2 = eff*p*h20^2*P0/(10^((sigma)/10));
%    gam=[gam1 gam2];
%  gam=[rho1_10 rho2]  ; 
%         [tau_opt,R_tau] = opt_P2(gam);
%  
% 
% R_common(i)=min(R_tau);
% R_final(i,:)=R_tau;
%     
%     i
%     
    
end    



grid on
P=0:5:20;

plot(P,co_R_common,'dr-',P,R_common,'or-');
%axis([0, 30,0,4])
%plot(P,co_R_common-R_common);
grid on
legend('with User cooperation','without user coperation')
xlabel('transmit Power,P0(dBm)')
ylabel('throuput (Mbps)')
hold off
  
