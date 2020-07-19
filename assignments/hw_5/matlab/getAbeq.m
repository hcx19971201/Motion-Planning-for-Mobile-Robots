function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    %
    %
    %
    %
    n_coef = n_order+1;
    
    Aeq_start(1,1:n_coef) = getRowVec(0,n_order,0);
    Aeq_start(2,1:n_coef) = getRowVec(0,n_order,1);
    Aeq_start(3,1:n_coef) = getRowVec(0,n_order,2);
    Aeq_start(4,1:n_coef) = getRowVec(0,n_order,3);
    
    beq_start = start_cond';
    
    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    %
    %
    %
    
    Aeq_end(1,n_coef*(n_seg-1)+1:n_coef*n_seg) = getRowVec(ts(end),n_order,0);
    Aeq_end(2,n_coef*(n_seg-1)+1:n_coef*n_seg) = getRowVec(ts(end),n_order,1);
    Aeq_end(3,n_coef*(n_seg-1)+1:n_coef*n_seg) = getRowVec(ts(end),n_order,2);
    Aeq_end(4,n_coef*(n_seg-1)+1:n_coef*n_seg) = getRowVec(ts(end),n_order,3);
    
    beq_end = end_cond';
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    %
    %
    %
    
    for i=1:n_seg-1
        Aeq_wp(i,n_coef*i+1:n_coef*(i+1)) = getRowVec(0,n_order,0);
        beq_wp(i) = waypoints(i+1);
    end
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    %
    %
    %
    
    for i=1:n_seg-1
        p1 = getRowVec(ts(i),n_order,0);
        p2 = getRowVec(0,n_order,0);
        
        Aeq_con_p(i,n_coef*(i-1)+1:n_coef*(i+1))=[p1,-p2];
    end
    
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %
    
    for i=1:n_seg-1
        v1 = getRowVec(ts(i),n_order,1);
        v2 = getRowVec(0,n_order,1);
        
        Aeq_con_v(i,n_coef*(i-1)+1:n_coef*(i+1))=[v1,-v2];
    end

    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %
    
    for i=1:n_seg-1
        a1 = getRowVec(ts(i),n_order,2);
        a2 = getRowVec(0,n_order,2);
        
        Aeq_con_a(i,n_coef*(i-1)+1:n_coef*(i+1))=[a1,-a2];
    end
    
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    
    for i=1:n_seg-1
        j1 = getRowVec(ts(i),n_order,3);
        j2 = getRowVec(0,n_order,3);
        
        Aeq_con_j(i,n_coef*(i-1)+1:n_coef*(i+1))=[j1,-j2];
    end
    
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end


function row = getRowVec(T,n_order,k)
    row = zeros(1,n_order+1);
    for i=k+1:n_order+1
        row(i) = prod(i-k:i-1)*T^(i-k-1);
    end
end