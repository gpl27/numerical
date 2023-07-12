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

function [x] = bairstow(a, alpha0, beta0, TOL, N)
    for j = 1:N
        n = length(a);
        b(1) = a(1);
        b(2) = a(2) + alpha0*b(1);
        for i = 3:n
            b(i) = a(i) + alpha0*b(i-1) + beta0*b(i-2);
        end
        c(1) = b(1);
        c(2) = b(2) + alpha0*c(1);
        for i = 3:n
            c(i) = b(i) + alpha0*c(i-1) + beta0*c(i-2);
        end
        // Calcular alpha1 e beta1
        delta = inv([c(n-2), c(n-3); c(n-1), c(n-2)])*[-b(n-1);-b(n)];
        alpha0 = alpha0 + delta(1);
        beta0 = beta0 + delta(2);
        if (b(n) <= TOL && b(n-1) <= TOL)
            x(1) = (alpha0 + sqrt(alpha0**2 + 4*beta0))/2;
            x(2) = (alpha0 - sqrt(alpha0**2 + 4*beta0))/2;
            x = return(x)
        end
    end
    error('Número máximo de iterações excedido!')
endfunction

// Jacobi
x = 0; y = 0;
for i = 1:20
    xn = (1+9*y)/10;
    yn = (1+9*x)/10;
    x = xn; y = yn;
end
// Gauss-Seidel
x = 0; y = 0;
for i = 1:20
    x = (1+9*y)/10;
    y = (1+9*x)/10;
end

// Gauss-Seidel
x = 0; y = 0; w = 4/3;
for i = 1:20
    xn = (1 - w)*x + w*((9*y + 1)/10);
    yn = (1 - w)*y + w*((9*xn + 1)/10);
    x = xn; y = yn;
end

