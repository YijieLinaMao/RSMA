function [cap_tot  p_1 p_2 p_3]=NOMA_update_P(weights,H,SNR,g,W);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Then we try to use CVX to solve the problem        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    u3=weights(3);
    
    [A_U,NT,N_user] = size(H); 
    
    h1=H(:,:,1);
    h2=H(:,:,2);
    h3=H(:,:,3);
    
    g_1=g{1};
    g_2_1=g{2};
    g_2=g{3};
    g_3_1=g{4};    
    g_3_2=g{5};
    g_3=g{6};
    
    W_1=W{1};
    W_2_1=W{2};
    W_2=W{3};
    W_3_1=W{4};    
    W_3_2=W{5};
    W_3=W{6};
    
    
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_3(NT,A_U) complex
    expression constraints;
    
    %T_k calculation
    T_1=square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
    
    T_2_1=square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
    T_2=square_abs(h2*p_2)+square_abs(h2*p_3)+1;
    
    T_3_1=square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;
    T_3_2=square_abs(h3*p_2)+square_abs(h3*p_3)+1;
    T_3=square_abs(h3*p_3)+1;
    
    
    %MSE calculation
    E_1=abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
    E_2_1=abs(g_2_1)^2*T_2_1-2*real(g_2_1*h2*p_1)+1;
    E_2=abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
    
    E_3_1=abs(g_3_1)^2*T_3_1-2*real(g_3_1*h3*p_1)+1;
    E_3_2=abs(g_3_2)^2*T_3_2-2*real(g_3_2*h3*p_2)+1;
    E_3=abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;
    
    %Rate-WMMSE relationship
    Xi_1=W_1*E_1-log2(W_1);   
    Xi_2_1=W_2_1*E_2_1-log2(W_2_1);
    Xi_2=W_2*E_2-log2(W_2);
    
    Xi_3_1=W_3_1*E_3_1-log2(W_3_1);
    Xi_3_2=W_3_2*E_3_2-log2(W_3_2);
    Xi_3=W_3*E_3-log2(W_3);
    
    Xi_1_users(1)=Xi_1;
    Xi_1_users(2)=Xi_2_1;
    Xi_1_users(3)=Xi_3_1;
    
    Xi_2_users(1)=Xi_2;
    Xi_2_users(2)=Xi_3_2;
    
    %common rate part:
    WMMSE_1=max(Xi_1_users);
    WMMSE_2=max(Xi_2_users);
    WMMSE_3=Xi_3;
    
    object_func=u1*(WMMSE_1)+u2*(WMMSE_2)+u3*(WMMSE_3);

    
    
    %objective function
    minimize(object_func)
    
    
    
    %constraints    
    constraints(19)=trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)-SNR;

    
    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end

    cap_tot=object_func;
    
    
end