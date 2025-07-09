function param = RLS_fcn(mag_raw)
    persistent P
    persistent Y
    persistent is_first
    persistent count
    persistent size
    persistent init_flag
    
    if isempty(is_first)
        P = NaN([5, 5]);
        Y = NaN([5, 1]);
        count = 1;
        size = 5;
        init_flag = true;
    end
    
    input = [mag_raw.mx*mag_raw.my, mag_raw.my.^2, mag_raw.mx, mag_raw.my, 1];
    output = -mag_raw.mx.^2;
    if count > size
        if init_flag
            param = P \ Y;
            P = I / (P'*P);
            init_flag = false;
        end
        e_k = output - input * param;
        P = P - (P * input') / (1 + input * P * input') * (input * P);
        param = param + P * input' * e_k;
    else
        P(count, :) = input;
        Y(count) = -mag_raw.mx.^2;
        param = NaN(5, 1);
    end
    count = count + 1;
end