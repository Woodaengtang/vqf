function est_param = rls_mag(mag)
    persistent N L init_flag rls_flag count P_init Y_init P_cov output

    if isempty(init_flag)
        count = 1;
        N = 9;
        L = 100;
        rls_flag = false;
        P_init = NaN([L, N]);
        output = -1;
    end
    
    input = input_vec(mag);
    if rls_flag
        
        
    else
        P_init(count, :) = input;
        Y_init(count) = output;
        count = count + 1;
        if count > 100
            y = -1*ones([L, 1]);
            est_param = (P_init'*P_init) \ (P_init'*y);
            P_cov = eye(N) / (P_init'*P_init);
        end
    end
    
end

function input = input_vec(mag)
    input = [mag(1)^2, mag(2)^2, mag(3)^2, 2*mag(1)*mag(2), 2*mag(1)*mag(3), 2*mag(2)*mag(3), mag(1), mag(2), mag(3)];
end

