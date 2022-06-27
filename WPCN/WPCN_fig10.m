clc
clear all


K=2;
D1=5;
D2=10;
a=2;
eff=0.5;

K = 2; % number of users
xi = 0.5; % energy harvesting efficiency
Gamma  = 9.8; % SNR gap, dB
D1 = 5; % distance from user i to the AP
D2 = 10;
alpha = 2; % channel pathloss exponent
sigma_1 = -160; % noise power spectral density, dBm/HZ;
sigma = -160 + 10*log10(10^6); % noise power, dBm;
h_1 = 10^(-3)*D1^(-alpha);
g_1 = 10^(-3)*D1^(-alpha);
h_2 = 10^(-3)*D2^(-alpha);
g_2 = 10^(-3)*D2^(-alpha);
PA = 0:5:25;

%N_PA = length(PA);
gamma1 = xi*h_1*g_1*10.^(PA./10)/(10^((Gamma+sigma)/10));
gamma2 = xi*h_2*g_2*10.^(PA./10)/(10^((Gamma+sigma)/10));

% 
% R_o1_p1 = zeros(1, N_PA); R_o2_p1 = zeros(1, N_PA); R_sum_o2_p1 = zeros(1, N_PA);
% R_sum_Norm_p1 = zeros(1, N_PA);
% R_o1_p2 = zeros(1, N_PA); R_o2_p2 = zeros(1, N_PA); R_min = zeros(1, N_PA);
for i = 1 :length(PA)
    % for P1, we can find the optimal tau based on Proposition 3.1
    gamma_p = [gamma1(i),gamma2(i)];
    
    [tau_p1] = Opt_P1(gamma_p);
    
    
    R_o1_p1(i) =  tau_p1(2)*log2(1+gamma1(i)*tau_p1(1)/tau_p1(2));
    R_o2_p1(i) =  tau_p1(3)*log2(1+gamma2(i)*tau_p1(1)/tau_p1(3));
    R_sum_o2_p1(i) = R_o1_p1(i) + R_o2_p1(i);
    R_sum_Norm_p1(i) = R_sum_o2_p1(i)/K;

    % for P2, using the iterative algorithm to obtain the optimal solution
    [tau_opt,R_tau] = opt_P2(gamma_p);
 
    R_o1_p2(i) =  R_tau(1);
    R_o2_p2(i) =  R_tau(2);
    R_min(i) = min(R_tau);
    
end
figure

plot(PA,R_min,'-db',PA,R_sum_o2_p1,':or',PA,R_sum_Norm_p1,'-or',PA,R_o1_p1,'--or',PA,R_o2_p1,'-.or','linewidth',2);
xlabel('P_A(dBm)','fontsize',12);
ylabel('Average Thoughput (Mbps)','fontsize',12);
legend('Max.Common-Thoughput for(P2)','Max.Sum-Thoughput for(P1)','Normalized Max.Sum-Thoughput',...
    'Thoughput of U_1 in (P1)','Thoughput of U_2 in (P1)');
grid on;
title('Fig. 10 in the Paper (Iterative Algorithm)');