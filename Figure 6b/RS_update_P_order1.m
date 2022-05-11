function [cap_tot, p_1, p_2, p_c ]=RS_update_P_order1(weights,H,SNR,g_c_1,g_c_2,g_1,g_2,W_c_1,W_c_2,W_1,W_2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Then we try to use CVX method to solve the problem     % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    
    [A_U,NT,N_user] = size(H); 
    h1=H(:,:,1);
    h2=H(:,:,2);
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_c(NT,A_U) complex
%     variable ro_c(N_user)  
    expression constraints;
    
    %objective function
    %common rate part:
    T_c_1=square_abs(h1*p_c)+square_abs(h1*p_1)+square_abs(h1*p_2)+1;
    T_c_2=square_abs(h2*p_c)+square_abs(h2*p_1)+square_abs(h2*p_2)+1;
    
    E_c_1=square_abs(g_c_1)*T_c_1-2*real(g_c_1*h1*p_c)+1;
    E_c_2=square_abs(g_c_2)*T_c_2-2*real(g_c_2*h2*p_c)+1;
    
    c_1=(W_c_1*E_c_1-log2(W_c_1));
    c_2=(W_c_2*E_c_2-log2(W_c_2));
    c=max(c_1,c_2);
    
    %we assume common rate is for user 1, therefore, c_2=0
    object_common=u1*c;
    
    %Private rate part
    T_1=square_abs(h1*p_1)+square_abs(h1*p_2)+1;
    T_2=square_abs(h2*p_1)+square_abs(h2*p_2)+1;
    E_1=abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
    E_2=abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
    object_private=u1*(W_1*E_1-log2(W_1))+u2*(W_2*E_2-log2(W_2));
    object_func=object_common+object_private;
    
    minimize(object_func)
    
    
    
    %constraints 
    constraints=trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_c'*p_c)-SNR;

    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end
    cap_tot=object_func;
    
    
end