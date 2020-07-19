function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % STEP 3.2.1 p constraint
    Aieq_p = [];
    bieq_p = [];
    
    for k=1:n_seg
        Aieq_p = blkdiag(Aieq_p, ts(k)*eye(n_order+1));
        bieq_p = [bieq_p; ones(n_order+1,1)*corridor_range(k,2)];
    end
    
    for k=1:n_seg
        bieq_p = [bieq_p; -ones(n_order+1,1)*corridor_range(k,1)];
    end
    
    Aieq_p = [Aieq_p;-Aieq_p];

    %#####################################################
    % STEP 3.2.2 v constraint   
    n_ctr = n_order;
    Aieq_v = [];

    for k = 1:n_seg
        A_k = zeros(n_ctr, n_order+1);
        for n = 1:n_ctr
            A_k(n, n:n+1) = n_order * [-1, 1] * ts(k)^(0);
        end
        Aieq_v = blkdiag(Aieq_v, A_k);
    end
    Aieq_v = [Aieq_v;-Aieq_v];
    bieq_v = ones(n_ctr*n_seg*2,1)* v_max;

    %#####################################################
    % STEP 3.2.3 a constraint   
    n_ctr = n_order-1;
    Aieq_a = [];
    
    for k = 1:n_seg
        A_k = zeros(n_ctr, n_order+1);
        for n = 1:n_ctr
            A_k(n, n:n+2) = n_order * (n_order-1) * [1, -2, 1] * ts(k)^(-1);
        end
        Aieq_a = blkdiag(Aieq_a, A_k);
    end
    Aieq_a = [Aieq_a;-Aieq_a];
    bieq_a = ones(n_ctr*n_seg*2,1)* a_max;
    
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a];
    bieq = [bieq_p; bieq_v; bieq_a];
%     Aieq = Aieq_p;
%     bieq = bieq_p;
end