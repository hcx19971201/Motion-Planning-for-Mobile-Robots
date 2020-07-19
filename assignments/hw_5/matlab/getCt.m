function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    %
    %
    %
    %
    
    %% continuity constraints
    C0 = eye(4);
    Ck = eye(4);
    
    fblock = [1;0;0;0;1;0;0;0];
    Cwf = [];
    for i = 1:(n_seg-1)
        Cwf = blkdiag(Cwf,fblock);
    end
    Cf = blkdiag(C0,Cwf,Ck);
    
    %% free variable
    Cp = [];
    pblock = [[0,0,0];eye(3);[0,0,0];eye(3)];
    for i = 1:(n_seg-1)
        Cp = blkdiag(Cp, pblock);
    end
    Cp = [zeros(4,3*(n_seg-1));Cp;zeros(4,3*(n_seg-1))];
    
    Ct = [Cf,Cp];
end