function [cap_tot, p_1, p_2  ]=NOMA_update_P_order1(weights,H,SNR,g_1_1,g_2_1,g_2_2,W_1_1,W_2_1,W_2_2)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Then we try to use CVX  to solve the problem       % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    
    [A_U,NT,N_user] = size(H); 
    h1=H(:,:,1);
    h2=H(:,:,2);
    
    cvx_begin quiet
    
   
    variable p_2(NT,A_U) complex
    variable p_1(NT,A_U) complex
    expression constraints(1,2*N_user);
    
    T_1_1=square_abs(h1*p_1)+square_abs(h1*p_2)+1;
    T_2_1=square_abs(h2*p_1)+square_abs(h2*p_2)+1;
    
    E_1_1=square_abs(g_1_1)*T_1_1-2*real(g_1_1*h1*p_1)+1;
    E_2_1=square_abs(g_2_1)*T_2_1-2*real(g_2_1*h2*p_1)+1;
    
    r_1_1=(W_1_1*E_1_1-log2(W_1_1));
    r_2_1 =(W_2_1*E_2_1-log2(W_2_1));
    
    %objective function
    object_common=u1*max(r_1_1,r_2_1);


    T_2_2=square_abs(h2*p_2)+1;
    E_2_2=abs(g_2_2)^2*T_2_2-2*real(g_2_2*h2*p_2)+1;
    object_private=u2*(W_2_2*E_2_2-log2(W_2_2));
    object_func=object_common+object_private;
    
    minimize(object_func)
    
    
    
    %constraints    
    constraints=trace(p_2'*p_2)+trace(p_1'*p_1)-SNR;
   
    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end
    cap_tot=object_func;
    
    
end