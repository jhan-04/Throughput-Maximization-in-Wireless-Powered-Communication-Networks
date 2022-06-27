clc
clear all
DK=10;
K=2:6;
PA=20;
a=2;
eff = 0.5; % efficiency
SNR_gap  = 9.8; % SNR gap
P_noise = -160 + 10*log10(10^6);



for i=1:length(K)
 n=1: K(i);%%i
 
 
D=DK.*n./K(i);
h = 10^(-3)*D.^(-a);
g= h;
gam= eff.*h.*g*10.^(PA./10)/(10^((SNR_gap+P_noise)/10))  
  

%%%%p1
[t0,t_p1] = Opt_P1_k(gam,K(i));%%%%t0를 포함하지 않음
    R_p1 =  t_p1.*log2(1+gam.*t0./t_p1);
    R_sum_p1(i) =sum(R_p1);
    R_Norm_p1(i) = R_sum_p1(i)/K(i);%%out put
    
    
  %%%ETA
  t_equal=1/(K(i)+1);
  R_ETA=t_equal*log2(1+gam.*t_equal/t_equal); 
  R_Norm_ETA(i)=sum(R_ETA)/K(i);%%out put
  R_min_ETA(i)=min(R_ETA);
  
  
  
  %%%%%%%%%%%%%P2
  
    [t, R_t2] = opt_P2_K(gam,K(i));
    R_common(i)=min(R_t2);
    
    
    
    

end


hold on
plot( K,R_Norm_p1,'-or' ,  K,R_Norm_ETA,'-sb',  K,R_min_ETA,'-.sb',K,R_common,'-db');
xlabel('PA(dBm)');
ylabel('Average Thoughput (Mbps)');
legend('Normalized Max.Sum-Thoughput' ,'R_Norm_ETA','R_min_ETA' ,'R_common');%,'Thoughput of U_1 in (P1)','Thoughput of U_2 in (P1)');
grid on;
hold off





