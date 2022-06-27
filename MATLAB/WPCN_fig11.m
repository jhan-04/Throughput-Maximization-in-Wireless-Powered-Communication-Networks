clear all
clc

PA=20;
a=2.5:0.5:4;
K = 2;
eff = 0.5; % efficiency
SNR_gap  = 9.8; % SNR gap
D1 = 5; 
D2 = 10;


P_noise = -160 + 10*log10(10^6); % noise power
h_1 = 10^(-3)*D1.^(-a);
g_1 = h_1;
h_2 = 10^(-3)*D2.^(-a);
g_2 = h_2;
gamma1 = eff.*h_1.*g_1*10.^(PA./10)/(10^((SNR_gap+P_noise)/10));  
gamma2 = eff.*h_2.*g_2*10.^(PA./10)/(10^((SNR_gap+P_noise)/10));



for i=1:length(a)
  
    gam = [gamma1(i),gamma2(i)];
    %%%%%%%%%P1
    [t_p1] = Opt_P1(gam);
    
    R1_p1(i) =  t_p1(2)*log2(1+gamma1(i)*t_p1(1)/t_p1(2));
    R2_p1(i) =  t_p1(3)*log2(1+gamma2(i)*t_p1(1)/t_p1(3));
    R_sum_p1(i) = R1_p1(i) + R2_p1(i);
    R_Norm_p1(i) = R_sum_p1(i)/K;

    %%%%%%%%P2
    [t,R_t] = opt_P2(gam);
 
    R1_p2(i) =  R_t(1);
    R2_p2(i) =  R_t(2);
    R_min(i) = min(R_t);


end



hold on
plot(a,R_min,'-db',  a,R_sum_p1,':or',   a,R_Norm_p1,'-or',  a,R1_p1,'--or',  a,R2_p1,'-.or');
xlabel('PA(dBm)');
ylabel('Average Thoughput (Mbps)');
legend('Max.Common-Thoughput for(P2)','Max.Sum-Thoughput for(P1)','Normalized Max.Sum-Thoughput',...
    'Thoughput of U_1 in (P1)','Thoughput of U_2 in (P1)');
grid on;
hold off