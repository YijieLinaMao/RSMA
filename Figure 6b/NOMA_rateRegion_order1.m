function [Capacity]=NOMA_rateRegion_order1(weights,H,SNRdB,tolerance)

SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H) ;
h1=H(:,:,1);
h2=H(:,:,2);


%Step 1: Initialization
P_k=(SNR)/N_user;
p_1=h1'/norm(h1)*sqrt(P_k);
p_2=h2'/norm(h2)*sqrt(P_k);


loop=1;
cap_tot_past=0;
count=0;
while (loop)
    
    %T_ck and T_k calculation
    T_1_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_2_1=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    T_2_2=abs(h2*p_2)^2+1;
    I_1_1=abs(h1*p_2)^2+1;
    I_2_1=abs(h2*p_2)^2+1;
    I_2_2=1;

    
    %Step 2: MMSE Common and private Receiver update at each user
    g_1_1=inv(T_1_1)*p_1'*h1'; 
    g_2_1=inv(T_2_1)*p_1'*h2';
    g_2_2=inv(T_2_2)*p_2'*h2';

    %Step 3: MSE Common and private Weight update at each user
    E_1_1=inv(T_1_1)*I_1_1; 
    E_2_1=inv(T_2_1)*I_2_1; 
    E_2_2=inv(T_2_2)*I_2_2;
        
    W_1_1=inv(E_1_1);
    W_2_1=inv(E_2_1);
    W_2_2=inv(E_2_2);
    

    %optimization
    [cap_tot,  p_1, p_2 ]=NOMA_update_P_order1(weights,H,SNR,g_1_1,g_2_1,g_2_2,W_1_1,W_2_1,W_2_2);
    
    if abs(cap_tot-cap_tot_past)<=tolerance
        loop=0;
    else
        cap_tot_past=cap_tot;
        count=count+1;
    end
    
    if count>=2000
        break;
    end
   
end


%Calculate the capacity of each UE
T_2_2=abs(h2*p_2)^2+1;
I_2_2=T_2_2-abs(h2*p_2)^2;
T_1_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
T_2_1=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
I_1_1=abs(h1*p_2)^2+1;
I_2_1=abs(h2*p_2)^2+1;
R_1_1=real(log2(T_1_1/I_1_1));
R_2_1=real(log2(T_2_1/I_2_1));

%Private Rate

R_2=real(log2(T_2_2/I_2_2));
Capacity(1)=min(R_1_1,R_2_1);
Capacity(2)=R_2;
% Capacity


          
          
          
 
