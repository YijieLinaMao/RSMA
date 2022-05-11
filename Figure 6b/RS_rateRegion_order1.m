function [ Capacity,P_common] = RS_rateRegion_order1(weights,H,SNRdB,tolerance)

SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H) ;
h1=H(:,:,1);
h2=H(:,:,2);
H_tot=[h1' h2'];

%Step 1: Initialization--MRC+SVD
P_common=SNR*0.9;
P_private_k=(SNR-P_common)/N_user;

[U2,~,~]=svd(H_tot);
hat_p_c=U2(:,1);

p_1=h1'/norm(h1)*sqrt(P_private_k);
p_2=h2'/norm(h2)*sqrt(P_private_k);
p_c=hat_p_c*sqrt(P_common); 


loop=1;
cap_tot_past=0;
count=0;
while (loop)

    %T_ck and T_k calculation
    T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    T_c_1=abs(h1*p_c)^2+T_1;
    T_c_2=abs(h2*p_c)^2+T_2;
    I_1=T_1-abs(h1*p_1)^2;
    I_2=T_2-abs(h2*p_2)^2;

    
    %Step 2: MMSE Common and private Receiver update at each user
    g_c_1=inv(T_c_1)*p_c'*h1'; 
    g_c_2=inv(T_c_2)*p_c'*h2';
    
    g_1=inv(T_1)*p_1'*h1'; 
    g_2=inv(T_2)*p_2'*h2';

    %Step 3: MSE Common and private Weight update at each user
    E_c_1=inv(T_c_1)*T_1; 
    E_c_2=inv(T_c_2)*T_2;
    
    E_1=inv(T_1)*I_1; 
    E_2=inv(T_2)*I_2;
        
    W_c_1=inv(E_c_1);
    W_c_2=inv(E_c_2);

    W_1=inv(E_1);
    W_2=inv(E_2);
    

    %Step 4: Update P 
    [cap_tot, p_1, p_2, p_c]=RS_update_P_order1(weights,H,SNR,g_c_1,g_c_2,g_1,g_2,W_c_1,W_c_2,W_1,W_2);

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

P_common=real(trace(p_c*p_c'));

%Calculate the capacity of each UE
T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
T_c_1=abs(h1*p_c)^2+T_1;
T_c_2=abs(h2*p_c)^2+T_2;
I_1=T_1-abs(h1*p_1)^2;
I_2=T_2-abs(h2*p_2)^2;

%Common Rate
R_c_1=real(log2(T_c_1/T_1));
R_c_2=real(log2(T_c_2/T_2));

%Private Rate
R_1=real(log2(T_1/I_1));
R_2=real(log2(T_2/I_2));

%ro_c is not a optimization variable in timesharing scenario
Capacity(1)=R_1+min(R_c_1,R_c_2);
Capacity(2)=R_2;
Capacity_common(1)=R_c_1;
Capacity_common(2)=R_c_1;




          
          
          
 