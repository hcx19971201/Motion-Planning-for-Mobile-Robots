function M = getM(n_seg, n_order, ts)
    M = [];
    
    for k = 1:n_seg
        M_k = [];
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        %
        %
        %
        %
        M_k = zeros(n_order+1);
        
        for i = 1:(n_order+1)/2
            M_k(i,i) = factorial(i-1);
        end
        
        T = ts(k);
        for i = (n_order+1)/2+1:(n_order+1)
            for j = (i-(n_order+1)/2):(n_order+1)
                M_k(i, j) = prod(j-(i-(n_order+1)/2-1):j-1)*T^(j-1-(i-(n_order+1)/2-1));
            end
        end
        
        M = blkdiag(M, M_k);
    end
end