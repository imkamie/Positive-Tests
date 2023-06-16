function [A, dy] =  fun_robertson(t, y)
    A = Kmatrix_robertson(y);

    dy = A*y;
end

function kY = Kmatrix_robertson(Y)
    % Laplacian matrix for robertson example
    kY = [-0.04 (10^4)*Y(3) 0;
        0.04 -3*(10^7)*Y(2)-(10^4)*Y(3) 0;
        0 3*(10^7)*Y(2) 0];
end