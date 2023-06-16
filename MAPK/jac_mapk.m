function [Ay, J] =  jac_mapk(t, y)
    A = Kmatrix_mapk(y);

    Ay = matrix_derivative_mapk();

    J = matrix_jac_mapk(t,y) + A;
end

function [alpha,k] = mapk_parameters
    alpha=1;
    k=[100/3; 1/3; 50; 0.5; 10/3; 0.1; 7/10];
end

function A = Kmatrix_mapk(Y)
    % Laplacian matrix for MAPK example
    [alpha,k] = mapk_parameters;
    A = [(-k(7)-k(1)*Y(2)) 0 0 k(2) 0 k(6); 
        0 (-k(1)*Y(1)) k(5) 0 0 0; 
        0 0 (-k(3)*Y(1)-k(5)) k(2) k(4) 0; 
        (1-alpha)*(k(1)*Y(2)) (alpha*k(1)*Y(1)) 0 (-k(2)) 0 0; 
        0 0 (k(3)*Y(1)) 0 (-k(4)) 0; 
        k(7) 0 0 0 0 (-k(6))];
end

function Ay = matrix_derivative_mapk()
    [alpha,k] = mapk_parameters;

    daY_dY1 = [0 0 0 0 0 0; 
        0 (-k(1)) 0 0 0 0; 
        0 0 (-k(3)) 0 0 0; 
        0 (alpha*k(1)) 0 0 0 0; 
        0 0 (k(3)) 0 0 0; 
        0 0 0 0 0 0];
    
    daY_dY2 = [(-k(1)) 0 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0; 
        (1-alpha)*(k(1)) 0 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0];
    
    Ay = daY_dY1 + daY_dY2;
end

function daY = matrix_jac_mapk(t, Y)
    [alpha,k] = mapk_parameters;

    daY_dY1 = [0 0 0 0 0 0; 
        0 (-k(1)) 0 0 0 0; 
        0 0 (-k(3)) 0 0 0; 
        0 (alpha*k(1)) 0 0 0 0; 
        0 0 (k(3)) 0 0 0; 
        0 0 0 0 0 0];
    
    daY_dY2 = [(-k(1)) 0 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0; 
        (1-alpha)*(k(1)) 0 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0];
    
    daY = daY_dY1*Y(1) + daY_dY2*Y(2);
end
