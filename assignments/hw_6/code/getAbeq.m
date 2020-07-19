function [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 2.1 p,v,a constraint in start
    n_constraint = length(start_cond);
    Aeq_start = zeros(n_constraint, n_all_poly);
    for k=1:n_constraint
        Aeq_start(k,1:k)=prod(n_order-k+2:n_order) * getRowVec(k) * ts(1)^(2-k);
    end
    
    beq_start = start_cond';
    
    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    n_constraint = length(end_cond);
    Aeq_end = zeros(n_constraint, n_all_poly);
    for k=1:n_constraint
        Aeq_end(k,n_all_poly-k+1:n_all_poly)=prod(n_order-k+2:n_order) * getRowVec(k) * ts(end)^(2-k);
    end
    
    beq_end = end_cond';
    
    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1,1);
    
    for k = 1:n_seg-1
        Aeq_con_p(k,k*(n_order+1)) = 1 * ts(k)^(1);
        Aeq_con_p(k,k*(n_order+1)+1) = -1 * ts(k+1)^(1);
    end

    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    Aeq_con_v =  zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1,1);
    
    for k = 1:n_seg-1
        Aeq_con_v(k,k*(n_order+1)-1:k*(n_order+1)) = n_order * [-1, 1] * ts(k)^(0);
        Aeq_con_v(k,k*(n_order+1)+1:k*(n_order+1)+2) = -1 * n_order * [-1, 1] * ts(k+1)^(0);
    end    

    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1,1);
    
    for k = 1:n_seg-1
        Aeq_con_a(k,k*(n_order+1)-2:k*(n_order+1)) = n_order * (n_order - 1) * [1, -2, 1] * ts(k)^(-1);
        Aeq_con_a(k,k*(n_order+1)+1:k*(n_order+1)+3) = -1 * n_order * (n_order - 1) * [1, -2, 1] * ts(k+1)^(-1);
    end  

    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start; Aeq_end; Aeq_con];
    beq = [beq_start; beq_end; beq_con];
end

function row = getRowVec(k)
    row = zeros(1, k);
    for i=1:k
        row(i) = (-1)^(i+k) * nchoosek(k-1,i-1);
    end
end