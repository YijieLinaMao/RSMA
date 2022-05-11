function Capacity=NOMA_rateRegion_order1(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance) 
SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H_BC_estimate) ;

randn('seed',28);

for i=1:M
    H_BC_real_1(:,:,i)=H_BC_estimate(:,:,1)+H_BC_error_1(:,:,i);
    H_BC_real_2(:,:,i)=H_BC_estimate(:,:,2)+H_BC_error_2(:,:,i);
end

%Step 1: Initialization
%The user decoded first is initialized by SVD. Therefore UE1 is initialized
%by SVD. 

h1=H_BC_estimate(:,:,1);
h2=H_BC_estimate(:,:,2);
H_tot=[h1' h2'];
P_k=(SNR)/N_user;
[U2,~,~]=svd(H_tot);
p_1=U2(:,1)*sqrt(P_k);

%The user decoded secondly is initialized by MRC. Therefore UE 2 is
%initialized by MRC.
p_2=H_BC_estimate(:,:,2)'/norm(H_BC_estimate(:,:,2))*sqrt(P_k);


loop=1;
cap_tot_past=0;
count=0;
while (loop)

    t_1_1=0;
    t_2_1=0;
    t_2_2=0;

    phi_1_1=0;
    phi_2_1=0;
    phi_2_2=0;

    f_1_1=zeros(NT,A_U);
    f_2_1=zeros(NT,A_U);
    f_2_2=zeros(NT,A_U);

    v_1_1=0;
    v_2_1=0;
    v_2_2=0;

    WW_1_1=0;
    WW_2_1=0;
    WW_2_2=0;
    %T_ck and T_k calculation
    for i=1:M
        h1=H_BC_real_1(:,:,i);
        h2=H_BC_real_2(:,:,i);
        
        T_1_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
        T_2_1=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
        T_2_2=abs(h2*p_2)^2+1;
        
        I_1_1=abs(h1*p_2)^2+1;
        I_2_1=abs(h2*p_2)^2+1;
        I_2_2=1;

        
        E_1_1=inv(T_1_1)*I_1_1; 
        E_2_1=inv(T_2_1)*I_2_1; 
        E_2_2=inv(T_2_2)*I_2_2;
    
        %Step 2: MMSE Common and private Receiver update at each user
        g_1_1=inv(T_1_1)*p_1'*h1'; 
        g_2_1=inv(T_2_1)*p_1'*h2';
        g_2_2=inv(T_2_2)*p_2'*h2';

        %Step 3: MSE Weight update at each user    
        W_1_1=inv(E_1_1);
        W_2_1=inv(E_2_1);
        W_2_2=inv(E_2_2);
        
        %Step 4: Average value calculation
        WW_1_1=WW_1_1+W_1_1;
        WW_2_1=WW_2_1+W_2_1;
        WW_2_2=WW_2_2+W_2_2;
        
        t_1_1=t_1_1+W_1_1*abs(g_1_1)^2;
        t_2_1=t_2_1+W_2_1*abs(g_2_1)^2;
        t_2_2=t_2_2+W_2_2*abs(g_2_2)^2;
        
        phi_1_1=phi_1_1+(W_1_1*abs(g_1_1)^2)*h1'*h1;  %same  as used in the Hamdi's paper
        phi_2_1=phi_2_1+(W_2_1*abs(g_2_1)^2)*h2'*h2;
        phi_2_2=phi_2_2+(W_2_2*abs(g_2_2)^2)*h2'*h2;

        f_1_1=f_1_1+W_1_1*h1'*g_1_1';
        f_2_1=f_2_1+W_2_1*h2'*g_2_1';
        f_2_2=f_2_2+W_2_2*h2'*g_2_2';

        
        v_1_1=v_1_1+log2(W_1_1);
        v_2_1=v_2_1+log2(W_2_1);
        v_2_2=v_2_2+log2(W_2_2);

    end
    WW_1_1=WW_1_1./M; 
    WW_2_1=WW_2_1./M;
    WW_2_2=WW_2_2./M;
    
    t_1_1=t_1_1./M; 
    t_2_1=t_2_1./M;
    t_2_2=t_2_2./M;
    
    phi_1_1=phi_1_1./M; 
    phi_2_1=phi_2_1./M;
    phi_2_2=phi_2_2./M;
    
    f_1_1=f_1_1./M; 
    f_2_1=f_2_1./M;
    f_2_2=f_2_2./M;
    
    v_1_1=v_1_1./M; 
    v_2_1=v_2_1./M;
    v_2_2=v_2_2./M;

    %optimization
    [cap_tot  p_1 p_2 ]=NOMA_update_P_order1(H_BC_estimate,M,weights,SNR,WW_1_1,WW_2_1,WW_2_2,t_1_1,t_2_1,t_2_2,...
                                     phi_1_1,phi_2_1,phi_2_2,f_1_1,f_2_1,f_2_2,v_1_1,v_2_1,v_2_2);
    

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


for i=1:M
    h1=H_BC_real_1(:,:,i);
    h2=H_BC_real_2(:,:,i);

    
    T_1_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
    T_2_1=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
    T_2_2=abs(h2*p_2)^2+1;
    
    I_1_1=abs(h1*p_2)^2+1;
    I_2_1=abs(h2*p_2)^2+1;
    I_2_2=1;
    
    R_1_1=real(log2(T_1_1/I_1_1));
    R_2_1=real(log2(T_2_1/I_2_1));
    R_1(i)=min(R_1_1,R_2_1);
    R_2(i)=real(log2(T_2_2/I_2_2));


end

%Calculate the rate of each UE
Capacity(1)=sum(R_1)/M;
Capacity(2)=sum(R_2)/M;






          
          
          
 
