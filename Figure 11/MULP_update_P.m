function [cap_tot,p_1,p_2]=MULP_update_P(M,weights,H_BC_estimate,SNR,W_1,W_2,t_1,t_2,...
                          phi_1,phi_2,f_1,f_2,v_1,v_2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Then we try to use CVX to solve the problem        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    [A_U,NT,N_user] = size(H_BC_estimate);  
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
 
    expression constraints;

    private_1=quad_form(p_1,phi_1)+quad_form(p_2,phi_1)+t_1-2*real(f_1'*p_1)+W_1-v_1;
    private_2=quad_form(p_1,phi_2)+quad_form(p_2,phi_2)+t_2-2*real(f_2'*p_2)+W_2-v_2;

    object_1=u1*private_1;
    object_2=u2*private_2;
    

    object_func=object_1+object_2;
    minimize(object_func)

    
    %constraints 
    constraints=trace(p_1'*p_1)+trace(p_2'*p_2)-SNR;

    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end

    cap_tot=object_func;
    
    
end