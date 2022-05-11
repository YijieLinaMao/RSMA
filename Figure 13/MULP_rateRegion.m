function WSR = MULP_rateRegion(weights,H,SNRdB,tolerance)

SNR = 10^(SNRdB/10);
[A_U,NT,N_user] = size(H); 
h1=H(:,:,1);
h2=H(:,:,2);
h3=H(:,:,3);

%Step 1: Initialization
% power splitting for three messages
alpha_1=0.4;
alpha_2=0.3;
alpha_3=0.3;

%Step 1: Initialization
%precoder in layer 3----use MRC
p_1=h1'/norm(h1)*sqrt(SNR*alpha_1);
p_2=h2'/norm(h2)*sqrt(SNR*alpha_2);
p_3=h3'/norm(h3)*sqrt(SNR*alpha_3);


loop=1;
cap_tot_past=0;
count=0;
while (loop)
    
    %T_k calculation
    T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;    
    T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;    
    T_3=abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;

   

    %Step 2: MMSE Common and private Receiver update at each user
    g_1=inv(T_1)*p_1'*h1';
    g_2=inv(T_2)*p_2'*h2';
    g_3=inv(T_3)*p_3'*h3';
    
    g{1}=g_1;
    g{2}=g_2;
    g{3}=g_3;

    %Step 3: MSE Common and private Weight update at each user
    E_1=inv(T_1)*(T_1-abs(h1*p_1)^2); 
    E_2=inv(T_2)*(T_2-abs(h2*p_2)^2); 
    E_3=inv(T_3)*(T_3-abs(h3*p_3)^2); 
        
    W_1=inv(E_1);
    W_2=inv(E_2);
    W_3=inv(E_3);
    
    W{1}=W_1;
    W{2}=W_2;
    W{3}=W_3;

    %optimize both power and ro_c   
    [cap_tot  p_1 p_2 p_3]=MULP_update_P(weights,H,SNR,g,W);

    if abs(cap_tot-cap_tot_past)<=tolerance
        loop=0;
    else
        cap_tot_past=cap_tot;
        count=count+1;
    end
%     count
    if count>=2000
        break;
    end
   
end

%Calculate the capacity of each UE
T_1=abs(h1*p_1)^2+abs(h1*p_2)^2+abs(h1*p_3)^2+1;
T_2=abs(h2*p_1)^2+abs(h2*p_2)^2+abs(h2*p_3)^2+1;
T_3=abs(h3*p_1)^2+abs(h3*p_2)^2+abs(h3*p_3)^2+1;


E_1=inv(T_1)*(T_1-abs(h1*p_1)^2); 
E_2=inv(T_2)*(T_2-abs(h2*p_2)^2); 
E_3=inv(T_3)*(T_3-abs(h3*p_3)^2); 



 % Rate of each user
R_1=-real(log2(E_1));
R_2=-real(log2(E_2));
R_3=-real(log2(E_3));

% WSR
WSR=weights(1)*R_1+weights(2)*R_2+weights(3)*R_3;


          
          
          
 
