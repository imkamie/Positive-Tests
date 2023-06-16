function [Ay, J] =  jac_robertson(t, y)
    A = Kmatrix_robertson(y);
    
    Ay = matrix_derivative_robertson();

    J = matrix_jac_robertson(t,y) + A;
end

function kY = Kmatrix_robertson(Y)
    % Laplacian matrix for robertson example
    kY = [-0.04 (10^4)*Y(3) 0;
        0.04 -3*(10^7)*Y(2)-(10^4)*Y(3) 0;
        0 3*(10^7)*Y(2) 0];
end

function Ay = matrix_derivative_robertson()  
    daY_dY2 = [0 0 0;
        0 -3*(10^7) 0;
        0 3*(10^7) 0];

    daY_dY3 = [0 10^4 0;
        0 -10^4 0;
        0 0 0];
    
    Ay = daY_dY2 + daY_dY3;
end

function daY = matrix_jac_robertson(t, Y)  
    daY_dY2 = [0 0 0;
        0 -3*(10^7) 0;
        0 3*(10^7) 0];

    daY_dY3 = [0 10^4 0;
        0 -10^4 0;
        0 0 0];
    
    daY = daY_dY2*Y(2) + daY_dY3*Y(3);
end