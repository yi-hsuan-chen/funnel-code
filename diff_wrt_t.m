function dXdt   = diff_wrt_t(i,X,dt,n)
    switch i
        case 1
            % use FORWARD difference here for the first point
            dXdt        = (-3*X{i}+4*X{i+1}-X{i+2})/(2*dt);
        case n
            % use BACKWARD difference here for the last point
            dXdt        = (3*X{i}-4*X{i-1}+X{i-2})/(2*dt);
        otherwise
            % use CENTRAL difference
            dXdt        = (X{i+1}-X{i-1})/(2*dt);
    end
end