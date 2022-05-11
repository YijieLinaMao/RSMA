function Capacity = DPC_rateRegion(weights,H,SNRdB,tolerance)

%------------------Algorithm proposed in Paper: Downlink Capacity
%Evaluation of Cellular Networks With Know-Interference Cancellation
%(Algorithm 1)



SNR = 10^(SNRdB/10);

[NT,A_U,N_user] = size(H) ;

%The algorithm is based on descend weights sequence of UEs
weights=weights(:);
[sortedWeights, index_UE]=sort(weights,'descend');
H = H(:,:,index_UE);  %sorted H according to descend weights

%Step 1: Q initialization
Q = zeros(A_U,A_U,N_user);       


gradient = zeros(A_U,A_U,N_user);  %initialized gradient


loop=1;
cap_tot_past=0;
count=0;
%The algorithm iterates "Iterations" times to make Q converges
while(loop)
    
    %covariance of each UE
    for i_user=1:N_user
        HQH(:,:,i_user)=H(:,:,i_user)*Q(:,:,i_user)*H(:,:,i_user)';
    end
    
    % Calculate the gradient of each UE on Q_i
    for i_user=1:N_user
        gradient_i_part = zeros(A_U,A_U);
        for j_user=i_user:N_user-1
            gradient_i_part = gradient_i_part+(sortedWeights(j_user)-sortedWeights(j_user+1))*...
                (H(:,:,i_user)'*inv( eye(NT)+sum(HQH(:,:,1:j_user),3))*H(:,:,i_user));
        end        
        gradient(:,:,i_user) = gradient_i_part + sortedWeights(N_user)*H(:,:,i_user)'*inv(eye(NT)+sum(HQH,3))*H(:,:,i_user);
    end
    
    %Step 2 & Step 3: principal eigenvector and the corresponding principal
    %eigenvalues of the gradient
    
    for i_user=1:N_user
        [V,D]=eig(gradient(:,:,i_user));
        [eigenValue(i_user) index_principal] = max(diag(D));
        eigenVector(:,i_user)=V(:,index_principal);
    end
    
    %Step 4: optimal user selection
    [~, i_opt_UE]=max(eigenValue);

     
    % Step 5: calculating t, bisection or Matlab optimization toolbox
    % Using MATLAB optimization toolbox, reference from 'J. Lee and N. Jindal, 
    %"Symmetric capacity of a downlink MIMO channel," IEEE ISIT 2006.'

    obj_t = @(t) -obj_vvh(t,i_opt_UE,Q,sortedWeights,H,SNR,eigenVector(:,i_opt_UE));
    t_opt=fminbnd(obj_t,0,1);
    
    % Update Q
    for i_user=1:N_user
        if (i_user==i_opt_UE)
            Q(:,:,i_user)=t_opt*Q(:,:,i_user)+(1-t_opt)*SNR*eigenVector(:,i_user)*eigenVector(:,i_user)';
        else
            Q(:,:,i_user)=t_opt*Q(:,:,i_user);
        end
    end
    
    for i_user=1:N_user
        info_UE(:,:,i_user) = eye(NT);
        for j_user=1:i_user
            info_UE(:,:,i_user) = info_UE(:,:,i_user) + H(:,:,j_user)*Q(:,:,j_user)*H(:,:,j_user)';
        end
    end

    Capacity = zeros(N_user,1);
    Capacity(index_UE(1)) = log2(real(det(info_UE(:,:,1))));
    for i_user=2:N_user
        Capacity(index_UE(i_user)) = log2(real(det(info_UE(:,:,i_user)))/real(det(info_UE(:,:,i_user-1))));
    end
    
    cap_tot_1=sum(Capacity);
    if abs(cap_tot_1-cap_tot_past)<=tolerance
        loop=0;
    else
        cap_tot_past=cap_tot_1;
        count=count+1;
    end
    
    if count>=2000
        break;
    end
   
end

% Calculate R
for i_user=1:N_user
    info_UE(:,:,i_user) = eye(NT);
    for j_user=1:i_user
        info_UE(:,:,i_user) = info_UE(:,:,i_user) + H(:,:,j_user)*Q(:,:,j_user)*H(:,:,j_user)';
    end
end

Capacity = zeros(N_user,1);
Capacity(index_UE(1)) = log2(real(det(info_UE(:,:,1))));
for i_user=2:N_user
    Capacity(index_UE(i_user)) = log2(real(det(info_UE(:,:,i_user)))/real(det(info_UE(:,:,i_user-1))));
end


%
% Sub-function
%
function f=obj_vvh(t,i_opt_UE,Q,weights,H,SNR,v);

[NT,A_U,N_user] = size(H); 

for i_user=1:N_user
    if(i_user==i_opt_UE)
        Q(:,:,i_user)=t*Q(:,:,i_user)+(1-t)*SNR*v*v';
    else
        Q(:,:,i_user)=t*Q(:,:,i_user);
    end;
end;

for i_user=1:N_user
    HQH(:,:,i_user)=H(:,:,i_user)*Q(:,:,i_user)*H(:,:,i_user)';
end;

cap_part1=0;
for i_user=1:N_user-1  
    cap_part1=cap_part1+(weights(i_user)-weights(i_user+1))*log2(det(eye(NT)+sum(HQH(:,:,1:i_user),3)));
end;

cap_part2=weights(N_user)*log2(det(eye(NT)+sum(HQH(:,:,1:N_user),3)));

f=cap_part1+cap_part2;