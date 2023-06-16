function [Ay, J] =  jac_stratospheric(t, y)
    A = Kmatrix_stratospheric(y);
    
    Ay = matrix_derivative_stratospheric(t);

    J = matrix_jac_stratospheric(t,y) + A;
end

function [k] = stratospheric_parameters(t)
    Tr=4.5;
    Ts=19.5;
    Tl=mod((t/3600),24);

    if((Tr<=Tl) && (Tl<=Ts))
        Ttmp = (2*Tl-Tr-Ts)/(Ts-Tr);
        sigma=(1/2)+(1/2)*cos(pi*(abs(Ttmp))*(Ttmp));
    else
        sigma=0;
    end

    k=[(2.643E-10)*(sigma^3), 8.018E-17, (6.120E-04)*sigma, 1.576E-15, (1.070E-03)*(sigma^2), 7.110E-11, 1.200E-10, 6.062E-15, 1.069E-11, (1.289E-02)*sigma];

end

function aY = Kmatrix_stratospheric(Y, t)
    % Laplacian matrix for stratospheric example
    [k] = stratospheric_parameters(t);

    aY = [(-k(6)-k(7)*Y(3)) 0 k(5) 0 0 0;
        k(6) (-k(2)*Y(4)-k(4)*Y(3)-k(9)*Y(6)) k(3) 2*k(1) 0 k(10);
        0 (1/3)*k(2)*Y(4) (-k(3)-k(5)-k(4)*Y(2)-k(7)*Y(1)-k(8)*Y(5)) (2/3)*k(2)*Y(2) 0 0;
        (1/2)*k(7)*Y(3) (k(4)*Y(3)+(1/2)*k(9)*Y(6)) k(3)+k(5)+k(4)*Y(2)+k(7)*Y(1)+k(8)*Y(5)+(1/2)*k(7)*Y(1) (-k(1)-k(2)*Y(2)) 0 (1/2)*k(9)*Y(2);
        0 0 0 0 (-k(8)*Y(3)) k(10)+k(9)*Y(2);
        0 0 0 0 k(8)*Y(3) (-k(10)-k(9)*Y(2))];
end

function Ay = matrix_derivative_stratospheric(t)  
    [k] = stratospheric_parameters(t);
    daY_dY1 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(7) 0 0 0;
        0 0 (3/2)*k(7) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    daY_dY2 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(4) (2/3)*k(2) 0 0;
        0 0 k(4) -k(2) 0 (1/2)*k(9);
        0 0 0 0 0 k(9);
        0 0 0 0 0 -k(9)];

    daY_dY3 = [-k(7) 0 0 0 0 0;
        0 -k(4) 0 0 0 0;
        0 0 0 0 0 0;
        (1/2)*k(7) k(4) 0 0 0 0;
        0 0 0 0 -k(8) 0;
        0 0 0 0 k(8) 0];

    daY_dY4 = [0 0 0 0 0 0;
        0 -k(2) 0 0 0 0;
        0 (1/3)*k(2) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    daY_dY5 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(8) 0 0 0;
        0 0 k(8) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    daY_dY6 = [0 0 0 0 0 0;
        0 -k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 (1/2)*k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    Ay = daY_dY1 + daY_dY2 + daY_dY3 + daY_dY4 + daY_dY5 + daY_dY6;
end

function daY = matrix_jac_stratospheric(t, Y)  
    [k] = stratospheric_parameters(t);
    daY_dY1 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(7) 0 0 0;
        0 0 (3/2)*k(7) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    daY_dY2 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(4) (2/3)*k(2) 0 0;
        0 0 k(4) -k(2) 0 (1/2)*k(9);
        0 0 0 0 0 k(9);
        0 0 0 0 0 -k(9)];

    daY_dY3 = [-k(7) 0 0 0 0 0;
        0 -k(4) 0 0 0 0;
        0 0 0 0 0 0;
        (1/2)*k(7) k(4) 0 0 0 0;
        0 0 0 0 -k(8) 0;
        0 0 0 0 k(8) 0];

    daY_dY4 = [0 0 0 0 0 0;
        0 -k(2) 0 0 0 0;
        0 (1/3)*k(2) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    daY_dY5 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(8) 0 0 0;
        0 0 k(8) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    daY_dY6 = [0 0 0 0 0 0;
        0 -k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 (1/2)*k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    daY = daY_dY1*Y(1) + daY_dY2*Y(2) + daY_dY3*Y(3) + daY_dY4*Y(4) + daY_dY5*Y(5) + daY_dY6*Y(6);
end
