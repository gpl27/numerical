function [x] = bissecao(f, intervalo, TOL, N)
    a = intervalo(1);
    b = intervalo(2);
    for i = 1:N
        x = (a + b)/2;
        e = abs(b - a)/2;
        fx = f(x);
        if (fx == 0 || e < TOL)
            x = return(x);
        elseif (f(a)*fx < 0)
            b = x;
        elseif (fx*f(b) < 0)
            a = x;
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

function [x] = regulafalsi(f, x0, x1, TOL, N)
    for i = 1:N
        fx0 = f(x0);
        fx1 = f(x1);
        x2 = (x1*fx0 - x0*fx1)/(fx0 - fx1);
        fx2 = f(x2);
        e = abs(fx1 - fx0)/2;
        if (fx2 == 0 || e < TOL)
            x = return(x2);
        elseif (fx0*fx2 < 0)
            x1 = x2;
        elseif (fx2*fx1 < 0)
            x0 = x2;
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

function [x] = pontofixo(g, x0, delta, N)
    for i = 1:N
        x1 = g(x0);
        e = abs(x1 - x0);
        if (x0 == x1 || e < delta)
            x = return(x1);
        else
            x0 = x1;
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

function [x] = newton(f, fl, xn, delta, N)
    for i = 1:N
        xn1 = xn - (f(xn)/fl(xn));
        e = abs(xn1 - xn);
        if (f(xn1) == 0 || e < delta)
            x = return(xn1);
        else
            xn = xn1;
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

function [x] = secante(f, x0, x1, delta, N)
    // N deve ser maior ou igual a 2
    for i = 1:N
        f0 = f(x0);
        f1 = f(x1);
        x2 = x1 - f1*((x1 - x0)/(f1 - f0));
        e = abs(x2 - x1);
        if (f(x2) == 0 || e < delta)
            x = return(x2)
        else
            x0 = x1;
            x1 = x2;
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

// Exemplo solucao de sistemas
// 10x - 9y = 1
// -9x + 10y = 1

// Solucao usando Jacobi
x = 0; y = 0;
for i = 1:20
    xn = (1+9*y)/10;
    yn = (1+9*x)/10;
    x = xn; y = yn;
end
disp('Jacobi:')
x
y

// Solucao usando Gauss-Seidel
x = 0; y = 0;
for i = 1:20
    x = (1+9*y)/10;
    y = (1+9*x)/10;
end
disp('Gauss-Seidel:')
x
y

// Solucao usando SOR e Jacobi
x = 0; y = 0; w = 4/3;
for i = 1:20
    xn = (1 - w)*x + w*((9*y + 1)/10);
    yn = (1 - w)*y + w*((9*x + 1)/10);
    x = xn; y = yn;
end
disp('SOR e Jacobi:')
x
y

// Solucao usando SOR e Gauss-Seidel
x = 0; y = 0; w = 4/3;
for i = 1:20
    xn = (1 - w)*x + w*((9*y + 1)/10);
    yn = (1 - w)*y + w*((9*xn + 1)/10);
    x = xn; y = yn;
end
disp('SOR e Gauss-Seidel:')
x
y
