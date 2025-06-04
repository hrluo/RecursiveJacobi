function [structure_array] = Reveal_recursive_structure(n, n_threshold, f)
    % structure_array = [recursion depth, block size, b], where b = n^f
    
    % Initialization
    recdepth = 1;
    b = round(n^f);
    structure_array = [recdepth, n, b];

    if 2*b <= n_threshold
        return;
    end
    if b >= round(n / 2)
        return;
    end
    
    while true
        n = 2*b;
        b = round(n^f);
        recdepth = recdepth + 1;
        structure_array = [structure_array; recdepth, n, b];

        if n <= n_threshold
            return;
        end
        if b >= round(n / 2)
            return;
        end
    end
end