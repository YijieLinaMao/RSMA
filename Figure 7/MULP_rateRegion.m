function Capacity = MULP_rateRegion(weights,H,SNRdB,tolerance)

%------------------Algorithm WMMSE Algorithm in Paper: [Weighted 
%Sum-Rate Maximization using Weighted MMSE for MIMO-BC 
%Beamforming Design]-------------------------------------------


SNR = 10^(SNRdB/10);

[A_U,NT,N_user] = size(H);  
H_tot=[];
for i_user = 1:N_user
  H_tot = cat(1,H_tot,H(:,:,i_user));  %size: (A_U*N_user)*NT
end

P_k=SNR/N_user;
for i_user=1:N_user
    h1=H(:,:,i_user);
    B(:,:,i_user)=h1'/norm(h1)*sqrt(P_k);
end

loop=1;
cap_tot_past=100000;
count=0;
while(loop)
    
    %Step 2: MMSE Receiver Calculation 
    for i_user =1 : N_user
%         B
        R(:,:,i_user)=eye(A_U);
        for j_user=1:N_user
            if j_user~=i_user
                R(:,:,i_user)=R(:,:,i_user)+H(:,:,i_user)*B(:,:,j_user)*B(:,:,j_user)'*H(:,:,i_user)';
            end
        end
        A_MMSE{i_user}=B(:,:,i_user)'*H(:,:,i_user)'*inv(H(:,:,i_user)*B(:,:,i_user)*B(:,:,i_user)'*H(:,:,i_user)'+R(:,:,i_user));
    end
    A_MMSE_blkdiag = blkdiag(A_MMSE{:});    %size: (A_U*N_user)*(A_U*N_user)

    %Step 3: MSE Weight
    for i_user=1:N_user
        E(:,:,i_user)=inv(eye(A_U)+B(:,:,i_user)'*H(:,:,i_user)'*inv(R(:,:,i_user))*H(:,:,i_user)*B(:,:,i_user));
        W{i_user}=weights(i_user)*inv(E(:,:,i_user));
    end
    W_blkdiag = blkdiag(W{:});  %size: (A_U*N_user)*(A_U*N_user)
    
    %Step 4: WMMSE transmit filter
    
    B_bar=inv(H_tot'*A_MMSE_blkdiag'*W_blkdiag*A_MMSE_blkdiag*H_tot+trace(W_blkdiag*A_MMSE_blkdiag*A_MMSE_blkdiag')*eye(NT)/SNR)*H_tot'*A_MMSE_blkdiag'*W_blkdiag;
    
    b=sqrt(SNR/trace(B_bar*B_bar'));
    B_MMSE=b*B_bar;   %size: NT*(A_U*N_user)
    for i_user=1:N_user
        B(:,:,i_user)=B_MMSE(:,(i_user-1)*A_U+1:i_user*A_U);
    end
    
    
    %loop determination
    for i_user=1:N_user
        R(:,:,i_user)=eye(A_U);
        for j_user=1:N_user
            if j_user~=i_user
                R(:,:,i_user)=R(:,:,i_user)+H(:,:,i_user)*B(:,:,j_user)*B(:,:,j_user)'*H(:,:,i_user)';
            end
        end
        E(:,:,i_user)=inv(eye(A_U)+B(:,:,i_user)'*H(:,:,i_user)'*inv(R(:,:,i_user))*H(:,:,i_user)*B(:,:,i_user));
        Capacity_loop(i_user)=real(log2(det(inv(E(:,:,i_user)))));
    end
    cap_tot=sum(Capacity_loop);
    if abs(cap_tot-cap_tot_past)<tolerance
        loop=0;
    else
        cap_tot_past=cap_tot;
        count=count+1;
    end
    
    if count>=2000
        break;
    end
end

%Calculate the rate of each UE
% trace(B(:,:,1)*B(:,:,1)')+trace(B(:,:,2)*B(:,:,2)')
% SNR

for i_user=1:N_user
    R(:,:,i_user)=eye(A_U);
    for j_user=1:N_user
        if j_user~=i_user
            R(:,:,i_user)=R(:,:,i_user)+H(:,:,i_user)*B(:,:,j_user)*B(:,:,j_user)'*H(:,:,i_user)';
        end
    end
    E(:,:,i_user)=inv(eye(A_U)+B(:,:,i_user)'*H(:,:,i_user)'*inv(R(:,:,i_user))*H(:,:,i_user)*B(:,:,i_user));
    Capacity(i_user)=real(log2(det(inv(E(:,:,i_user)))));
end



          
          
          
 