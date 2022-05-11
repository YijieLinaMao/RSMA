function [Capacity,P_common] = RS_rateRegion_order1(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance)  

SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H_BC_estimate);

randn('seed',28);

for i=1:M
    H_BC_real_1(:,:,i)=H_BC_estimate(:,:,1)+H_BC_error_1(:,:,i);
    H_BC_real_2(:,:,i)=H_BC_estimate(:,:,2)+H_BC_error_2(:,:,i);
end


%%%%%%%%%%%%%%%%%%%%%%initialization when bias is not 1%%%%%%%%%%%%%%%%%%%
%Step 1: Initialization--MRC+SVD
h1=H_BC_estimate(:,:,1);
h2=H_BC_estimate(:,:,2);
H_tot=[h1' h2'];
[U2,~,~]=svd(H_tot);
hat_p_c=U2(:,1);

P_common=SNR*0.8;
P_private_k=SNR-P_common;

[U2,~,~]=svd(H_tot);
hat_p_c=U2(:,1);

p_1=h1'/norm(h1)*sqrt(P_private_k*0.5);
p_2=h2'/norm(h2)*sqrt(P_private_k*0.5);
p_c=hat_p_c*sqrt(P_common);   


loop=1;
cap_tot_past=0;
count=0;
while (loop)
    t_c_1=0;
    t_c_2=0;
    t_1=0;
    t_2=0;
    phi_c_1=0;
    phi_c_2=0;
    phi_1=0;
    phi_2=0;
    f_c_1=zeros(NT,A_U);
    f_c_2=zeros(NT,A_U);
    f_1=zeros(NT,A_U);
    f_2=zeros(NT,A_U);
    v_c_1=0;
    v_c_2=0;
    v_1=0;
    v_2=0;
    WW_c_1=0;
    WW_c_2=0;
    WW_1=0;
    WW_2=0;
    for i=1:M
        
        h1=H_BC_real_1(:,:,i);
        h2=H_BC_real_2(:,:,i);
        
        T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+1;
        T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+1;
        T_c_1=abs(h1*p_c)^2+T_1;
        T_c_2=abs(h2*p_c)^2+T_2;
        I_1=T_1-abs(h1*p_1)^2;
        I_2=T_2-abs(h2*p_2)^2;
        
        E_c_1=inv(T_c_1)*T_1; 
        E_c_2=inv(T_c_2)*T_2;

        E_1=inv(T_1)*I_1; 
        E_2=inv(T_2)*I_2;


        %Step 2: MMSE Common and private Receiver update at each user
        g_c_1=inv(T_c_1)*p_c'*h1'; 
        g_c_2=inv(T_c_2)*p_c'*h2';

        g_1=inv(T_1)*p_1'*h1'; 
        g_2=inv(T_2)*p_2'*h2';

        %Step 3: MSE Common and private Weight update at each user

        W_c_1=inv(E_c_1);
        W_c_2=inv(E_c_2);

        W_1=inv(E_1);
        W_2=inv(E_2);
        
        WW_c_1=WW_c_1+W_c_1;
        WW_c_2=WW_c_2+W_c_2;
        WW_1=WW_1+W_1;
        WW_2=WW_2+W_2;
        
        %Step 4: Calculating the averaged value for optimization 
        t_c_1=t_c_1+W_c_1*abs(g_c_1)^2;
        t_c_2=t_c_2+W_c_2*abs(g_c_2)^2;
        t_1=t_1+W_1*abs(g_1)^2;
        t_2=t_2+W_2*abs(g_2)^2;
        
        phi_c_1=phi_c_1+(W_c_1*abs(g_c_1)^2)*h1'*h1;
        phi_c_2=phi_c_2+(W_c_2*abs(g_c_2)^2)*h2'*h2;
        phi_1=phi_1+(W_1*abs(g_1)^2)*h1'*h1;
        phi_2=phi_2+(W_2*abs(g_2)^2)*h2'*h2;
        

        f_c_1=f_c_1+W_c_1*h1'*g_c_1';
        f_c_2=f_c_2+W_c_2*h2'*g_c_2';
        f_1=f_1+W_1*h1'*g_1';
        f_2=f_2+W_2*h2'*g_2';
        
        v_c_1=v_c_1+log2(W_c_1);
        v_c_2=v_c_2+log2(W_c_2);
        v_1=v_1+log2(W_1);
        v_2=v_2+log2(W_2);
    end
    WW_c_1=WW_c_1./M;
    WW_c_2=WW_c_2./M;
    WW_1=WW_1./M;
    WW_2=WW_2./M;
    
    t_c_1=t_c_1./M;
    t_c_2=t_c_2./M;
    t_1=t_1./M;
    t_2=t_2./M;
    
    phi_c_1=phi_c_1./M;
    phi_c_2=phi_c_2./M;
    phi_1=phi_1./M;
    phi_2=phi_2./M;
    
    f_c_1=f_c_1./M;
    f_c_2=f_c_2./M;
    f_1=f_1./M;
    f_2=f_2./M;
    
    v_c_1=v_c_1./M;
    v_c_2=v_c_2./M;
    v_1=v_1./M;
    v_2=v_2./M;
    
    

    %Step 4: Update P
    [cap_tot,p_1,p_2,p_c]=RS_update_P_order1(H_BC_estimate,M,weights,SNR,WW_c_1,WW_c_2,WW_1,WW_2,t_c_1,t_c_2,t_1,t_2,...
                         phi_c_1,phi_c_2,phi_1,phi_2,f_c_1,f_c_2,f_1,f_2,v_c_1,v_c_2,v_1,v_2);

    if abs(cap_tot-cap_tot_past)<=tolerance
        loop=0;
    else
        cap_tot_past=cap_tot;
        count=count+1;
    end

    if count>=500
        loop=0;
        break;
        
    end

end

P_common=real(trace(p_c*p_c'));

%Calculate the rate of each UE
for i=1:M
    h1=H_BC_real_1(:,:,i);
    h2=H_BC_real_2(:,:,i);
        
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
    R_1(:,:,i)=real(log2(T_1/I_1))+min(R_c_1,R_c_2);
    R_2(:,:,i)=real(log2(T_2/I_2));
end


Capacity(1)=sum(R_1)/M;
Capacity(2)=sum(R_2)/M;




          
          
          
 