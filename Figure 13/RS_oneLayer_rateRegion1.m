function WSR = RS_oneLayer_rateRegion1(weights,H,SNRdB,tolerance)

SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H); 
h1=H(:,:,1);
h2=H(:,:,2);
h3=H(:,:,3);

%Step 1: Initialization
% power splitting 
alpha_3=0.2;
alpha_1=1-alpha_3;


%precoder in layer 3----use MRC with SVD
H_123=[h1' h2' h3'];
[U_123,~,~]=svd(H_123);
p_123=U_123(:,1)*sqrt(SNR*alpha_3);


%precoder in layer 1----use MRC
power_1=SNR*alpha_1/N_user; %power in layer 1 is splitted into three parts for three different precoder
p_1=h1'/norm(h1)*sqrt(power_1);
p_2=h2'/norm(h2)*sqrt(power_1);
p_3=h3'/norm(h3)*sqrt(power_1);


loop=1;
cap_tot_past=0;
count=0;
while (loop)
    
    %T_k calculation
    T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_3=abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

    T_1_123=abs(h1*p_123)^2+T_1;   
    T_2_123=abs(h2*p_123)^2+T_2;
    T_3_123=abs(h3*p_123)^2+T_3;
    
    
    %Step 2: MMSE Common and private Receiver update at each user
    g_1_123=inv(T_1_123)*p_123'*h1'; 
    g_1=inv(T_1)*p_1'*h1';
    
    g_2_123=inv(T_2_123)*p_123'*h2'; 
    g_2=inv(T_2)*p_2'*h2';
    
    g_3_123=inv(T_3_123)*p_123'*h3'; 
    g_3=inv(T_3)*p_3'*h3';
    
    g{1}=g_1_123;
    g{2}=g_1;    
    g{3}=g_2_123;
    g{4}=g_2;    
    g{5}=g_3_123;
    g{6}=g_3;

    %Step 3: MSE Common and private Weight update at each user
    E_1_123=inv(T_1_123)*(T_1_123-abs(h1*p_123)^2); 
    E_1=inv(T_1)*(T_1-abs(h1*p_1)^2); 
    
    E_2_123=inv(T_2_123)*(T_2_123-abs(h2*p_123)^2); 
    E_2=inv(T_2)*(T_2-abs(h2*p_2)^2); 
        
    E_3_123=inv(T_3_123)*(T_3_123-abs(h3*p_123)^2); 
    E_3=inv(T_3)*(T_3-abs(h3*p_3)^2); 

    W_1_123=inv(E_1_123);
    W_1=inv(E_1);
    
    W_2_123=inv(E_2_123);
    W_2=inv(E_2);
    
    W_3_123=inv(E_3_123);
    W_3=inv(E_3);
    
    W{1}=W_1_123;
    W{2}=W_1;
    
    W{3}=W_2_123;
    W{4}=W_2;
    
    W{5}=W_3_123;
    W{6}=W_3;


    %Step 4: Update P
    [cap_tot, p_1, p_2, p_3,  p_123, ro_1_123,ro_2_123,ro_3_123]=RS_oneLayer_update_P1(weights,H,SNR,g,W);

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


% Calculate the capacity of each UE
   %Private rate calculation
    T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
    T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
    T_3=abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

    E_1=inv(T_1)*(T_1-abs(h1*p_1)^2); 
    E_2=inv(T_2)*(T_2-abs(h2*p_2)^2); 
    E_3=inv(T_3)*(T_3-abs(h3*p_3)^2); 
   
    R_1=-real(log2(E_1));
    R_2=-real(log2(E_2));
    R_3=-real(log2(E_3));
    
    %Common rate calculation 
    C_1_123=-ro_1_123;
    C_2_123=-ro_2_123;
    C_3_123=-ro_3_123;

    
    % Rate of each user
    rate_1=C_1_123 + R_1;
    rate_2=C_2_123 + R_2;
    rate_3=C_3_123 + R_3;

    % WSR
    WSR=weights(1)*rate_1+weights(2)*rate_2+weights(3)*rate_3;


          
          
          
 