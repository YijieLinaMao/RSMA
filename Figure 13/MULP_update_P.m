function [cap_tot  p_1 p_2 p_3]=MULP_update_P(weights,H,SNR,g,W);


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
    g_2=g{2};
    g_3=g{3};

    
    W_1=W{1};
    W_2=W{2};
    W_3=W{3};

    
    
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_3(NT,A_U) complex
    expression constraints;
    
    %T_k calculation
    T_1=square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
    
    T_2=square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
    
    T_3=square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;

    
    %MSE calculation
    E_1=abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
    E_2=abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
    E_3=abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;
    
    %Rate-WMMSE relationship
    Xi_1=W_1*E_1-log2(W_1);   
    Xi_2=W_2*E_2-log2(W_2);
    Xi_3=W_3*E_3-log2(W_3);

    
    %common rate part:
    WMMSE_1=Xi_1;
    WMMSE_2=Xi_2;
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