function Capacity = MULP_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance)

%------------------Algorithm WMMSE Algorithm in Paper: [Weighted 
%Sum-Rate Maximization using Weighted MMSE for MIMO-BC 
%Beamforming Design]-------------------------------------------

[A_U,NT,N_user] = size(H_BC_estimate);  
SNR = 10^(SNRdB/10);
P_k=SNR/N_user;

%precoder initialization for each channel estimation
p_1=H_BC_estimate(:,:,1)'/norm(H_BC_estimate(:,:,1))*sqrt(P_k);
p_2=H_BC_estimate(:,:,2)'/norm(H_BC_estimate(:,:,2))*sqrt(P_k);

%1000 different channel errors results in 1000 real channel
for i=1:M
    H_BC_real_1(:,:,i)=H_BC_estimate(:,:,1)+H_BC_error_1(:,:,i);
    H_BC_real_2(:,:,i)=H_BC_estimate(:,:,2)+H_BC_error_2(:,:,i);
%     H_tot_real(:,:,i)=[H_BC_real_1' H_BC_real_2']'; %size: (A_U*N_user)*NT
end
    
    

loop=1;
cap_tot_past=100000;
count=0;
while(loop)

    t_1=0;
    t_2=0;

    phi_1=0;
    phi_2=0;

    f_1=zeros(NT,A_U);
    f_2=zeros(NT,A_U);

    v_1=0;
    v_2=0;

    WW_1=0;
    WW_2=0;
    
    for i=1:M
        h1=H_BC_real_1(:,:,i);
        h2=H_BC_real_2(:,:,i);
        
        T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
        T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
        
        I_1=T_1-abs(h1*p_1)^2;
        I_2=T_2-abs(h2*p_2)^2;
        
        E_1=inv(T_1)*I_1; 
        E_2=inv(T_2)*I_2;
    
        %Step 2: MMSE  Receiver update at each user
        g_1=inv(T_1)*p_1'*h1'; 
        g_2=inv(T_2)*p_2'*h2';

        %Step 3: MSE Weight update at each user    
        W_1=inv(E_1);
        W_2=inv(E_2);

        
        
        %Step 4: Calculating the averaged value for optimization 
        WW_1=WW_1+W_1;
        WW_2=WW_2+W_2;

        t_1=t_1+W_1*abs(g_1)^2;
        t_2=t_2+W_2*abs(g_2)^2;
        

        phi_1=phi_1+(W_1*abs(g_1)^2)*h1'*h1;
        phi_2=phi_2+(W_2*abs(g_2)^2)*h2'*h2;
        
        f_1=f_1+W_1*h1'*g_1';
        f_2=f_2+W_2*h2'*g_2';

        v_1=v_1+log2(W_1);
        v_2=v_2+log2(W_2);
    end

    WW_1=WW_1./M;
    WW_2=WW_2./M;
    
    t_1=t_1./M;
    t_2=t_2./M;

    phi_1=phi_1./M;
    phi_2=phi_2./M;
    
    f_1=f_1./M;
    f_2=f_2./M;
    
    v_1=v_1./M;
    v_2=v_2./M;
    

    %Step 4: Update P ---use subgradient method
    [cap_tot p_1 p_2]=MULP_update_P(M,weights,H_BC_estimate,SNR,WW_1,WW_2,t_1,t_2,...
                          phi_1,phi_2,f_1,f_2,v_1,v_2);

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
% count
% trace(p_1*p_1')+trace(p_2*p_2')



%Calculate the rate of each UE
for i=1:M
    h1=H_BC_real_1(:,:,i);
    h2=H_BC_real_2(:,:,i);

    T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+1;

    I_1=T_1-abs(h1*p_1)^2;
    I_2=T_2-abs(h2*p_2)^2;

    %Private Rate
    R_1(i)=real(log2(T_1/I_1));
    R_2(i)=real(log2(T_2/I_2));

end

Capacity(1)=sum(R_1)/M;
Capacity(2)=sum(R_2)/M;






          
          
          
 
          
          
          
 