// Ajuste de Curvas
function melhor_ajuste(x, y)
    for i = 1:(length(x) - 1)
        deltax(i) = x(i+1) - x(i);
        deltay(i) = y(i+1) - y(i);
        deltalnx(i) = log(x(i+1)) - log(x(i));
        deltalny(i) = log(y(i+1)) - log(y(i));
        deltax2(i) = x(i+1)**2 - x(i)**2;
        deltay2(i) = y(i+1)**2 - y(i)**2;
        deltax3(i) = x(i+1)**3 - x(i)**3;
        deltay3(i) = y(i+1)**3 - y(i)**3;
    end
    reta = deltay ./ deltax;
    parab = deltay2 ./ deltax2;
    cubica = deltay3 ./ deltax3;
    pot = deltalny ./ deltax;
    expon = deltalny ./ deltalnx;
    
    disp(reta, parab, cubica, expon, pot);
endfunction

// y = a1*x**p + a2*x**p-1 + ... + an
function [c] = regressao(x, y, p)
    n = p + 1;
    for i = 1:n
        for j = 1:n
            if (i == 1 && j == 1)
                R(i, j) = length(x);
            end
            R(i, j) = sum(x**(i + j - 2));
        end
        b(i) = sum(x**(i-1) .* y);
    end
    c = return(inv(R)*b);
endfunction

// y = a1*e**(-a2*x)
function [c] = regressao_exp(x, y)
    R = [length(x), sum(x);sum(x), sum(x**2)];
    b = [sum(log(y)); sum(x.*log(y))];
    c = inv(R)*b;
    c(1) = exp(c(1));
    c = return(inv(R)*b);
endfunction

// y = a*x**b
function [c] = regressao_pot(x, y)
    R = [length(x), sum(log(x)); sum(log(x)), sum(log(x)**2)];
    b = [sum(log(y)); sum(log(x).*log(y))];
    c = return(inv(R)*b);
endfunction

// Interpolacao polinomial
// 	simplesmente montar sistema com n+1 equacoes para polinomio de grau n

// Temos ainda a possibilidade de usar:
// 	interpolacao por polinomios ortogonais 
// (Legendre, Lagrange, Newton) quando os dados sao mais precisos
//	interpolacao por spline cubico (funcao cubica entre um ponto e outro)


// Derivacao e Integracao Numerica
// Devolve fi = f([x-2h, x-h, x, x+h, x+2h])
// Ou seja, f aplicada nas redondezas de x
function [fi] = func_para_tab(f, x, h)
    X = (x - 2*h):h:(x + 2*h);
    for i = 1:length(X)
        fi(i) = f(X(i));
    end
    fi = return(fi);
endfunction

// Aplica f aos x
function [fi] = calcula_pontos(f, x)
    for i = 1:length(x)
        fi(i) = f(x(i));
    end
    fi = return(fi);
endfunction

// O(h**2)
// Recebe f = [..., x-2h, x-h, x, x+h, x+2h, ...];
function [c] = derivada1(f, h)
    m = ceil(length(f)/2);
    c = return((f(m+1) - f(m-1))/(2*h));
endfunction
function [c] = derivada2(f, h)
    m = ceil(length(f)/2);
    c = return((f(m+1) - 2*f(m) + f(m-1))/(h**2));
endfunction
function [c] = derivada3(f, h)
    m = ceil(length(f)/2);
    c = return((f(m+2) - 2*f(m+1) + 2*f(m-1) - f(m-2))/(2*h**3));
endfunction
function [c] = derivada4(f, h)
    m = ceil(length(f)/2);
    c = return((f(m+2) - 4*f(m+1) + 6*f(m) - 4*f(m-1) + f(m-2))/(h**4));
endfunction

f = [1.7, 1.869, 2.037];
disp(derivada1(f, 0.1));

function [c] = trapezios(f, h)
    T = (h/2)*(f(1) + 2*sum(f(2:length(f)-1)) + f(length(f)))
    c = return(T);
endfunction
v = [4.2, 7.5, 9.0, 10.5, 7.0];
disp(trapezios(v, 0.1));

function [c] = simpson(f, h)
    S = (h/3)*(f(1) + 4*(sum(f(2:2:(length(f)-1)))) + 2*(sum(f(3:2:(length(f)-1)))) + f(length(f)));
    c = return(S);
endfunction
deff('y=f(x)', 'y = 1/(1+x**2)');
A = calcula_pontos(f, 0:0.1:1);
disp(simpson(A, 0.1));

// Existe ainda o metodo da Quadratura de Gauss
// que possui ordem 2n-1 para n pontos


// Equacoes diferenciais
// Euler O(h)
deff('yn = yl(t, y)', 'yn = 2 - t + 3*y');
h = 0.1; x = 0:h:0.2; y(1) = 1;
for i = 2:length(x)
   y(i) = y(i-1) + h*yl(x(i-1), y(i-1)); 
end
disp(y);

// Problema de exemplo:
// y em [0,1] com h = 0.1
//    yl = -x*y
//    y(0) = 1

// Runge-Kutta 2 - O(h**2)
deff('yn = yl(x, y)', 'yn = -x*y');
x = 0:0.1:1; h = 0.1; y(1) = 1;
for i = 2:length(x)
    k1 = h*yl(x(i-1), y(i-1));
    k2 = h*yl(x(i-1)+h, y(i-1)+k1);
    y(i) = y(i-1) + ((k1 + k2)/2);
end
disp(y);

// Runge-Kutta 4 - O(h**4)
deff('yn = yl(x, y)', 'yn = -x*y');
x = 0:0.1:1; h = 0.1; y(1) = 1;
for i = 2:length(x)
    k1 = yl(x(i-1), y(i-1));
    k2 = yl(x(i-1)+h/2, y(i-1)+(h/2)*k1);
    k3 = yl(x(i-1)+h/2, y(i-1)+(h/2)*k2);
    k4 = yl(x(i-1)+h, y(i-1)+h*k3);
    y(i) = y(i-1) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end
disp(y);
y = [];

// Adams / Predição Correção
// Problema de Exemplo:
// y(0.4) com h = 0.1
//    yl = -x*y
//    y(0) = 1

// Primeiro usar Runge-Kutta 4 para encontrar estimativas de [yn-3, yn]
deff('yn = yl(x, y)', 'yn = -x*y');
h = 0.1; x = 0:h:0.3; y(1) = 1; n = length(x);
for i = 2:n
    k1 = yl(x(i-1), y(i-1));
    k2 = yl(x(i-1)+h/2, y(i-1)+(h/2)*k1);
    k3 = yl(x(i-1)+h/2, y(i-1)+(h/2)*k2);
    k4 = yl(x(i-1)+h, y(i-1)+h*k3);
    y(i) = y(i-1) + (h/6)*(k1 + 2*k2 + 2*k3 + k4); 
end

// Usar formula de Adams-Bashforth para encontrar uma estimativa de yn+1 e portanto de fn+1 (predicao)
y(n+1) = y(n) + (h/24)*(55*yl(x(n), y(n)) - 59*yl(x(n-1), y(n-1)) + 37*yl(x(n-2), y(n-2)) - 9*yl(x(n-3), y(n-3)));
disp(y(n+1));

// Usar formula de Adams-Moulton para obter um valor mais preciso de yn+1 (correcao)
y(n+1) = y(n) + (h/24)*(9*yl(x(n)+h, y(n+1)) + 19*yl(x(n), y(n)) - 5*yl(x(n-1), y(n-1)) + yl(x(n-2), y(n-2)));
disp(y);

// Transforma-se equacoes diferenciais de ordens maiores em sistemas de equacoes diferenciais de primeira ordem
// Sistemas de equacoes diferenciais
// Problema de Exemplo:
//    xll(t) + 4*xl(t) + 5*x(t) = 0
//    x(0) = 3
//    xl(0) = -5
// substituindo xl(t) = y(t)
//    xl = y
//    yl = -5*x - 4*y
//    x(0) = 3
//    y(0) = -5
deff('xl = f(t, x, y)', 'xl = y');
deff('yl = g(t, x, y)', 'yl = -5*x - 4*y');
h = 0.1; t = 0:h:5; x(1) = 3; y(1) = -5; n = length(t);
for i = 2:n
    k1 = f(t(i-1), x(i-1), y(i-1));
    l1 = g(t(i-1), x(i-1), y(i-1));
    k2 = f(t(i-1)+h/2, x(i-1)+(h*k1)/2, y(i-1)+(h*l1)/2);
    l2 = g(t(i-1)+h/2, x(i-1)+(h*k1)/2, y(i-1)+(h*l1)/2);
    k3 = f(t(i-1)+h/2, x(i-1)+(h*k2)/2, y(i-1)+(h*l2)/2);
    l3 = g(t(i-1)+h/2, x(i-1)+(h*k2)/2, y(i-1)+(h*l2)/2);
    k4 = f(t(i-1)+h, x(i-1)+h*k3, y(i-1)+h*l3);
    l4 = g(t(i-1)+h, x(i-1)+h*k3, y(i-1)+h*l3);
    x(i) = x(i-1) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    y(i) = y(i-1) + (h/6)*(l1 + 2*l2 + 2*l3 + l4);
end
disp(x(1), x(11), x(n));


// Otimizacao
function [r] = aurea(f, a, b, TOL)
    phi = (1+sqrt(5))/2;
    L = b - a;
    while L > TOL
        atil = b - L/phi; btil = a + L/phi;
        fatil = f(atil); fbtil = f(btil);
        if fatil < fbtil
            b = btil;
        elseif fatil > fbtil
            a = atil;
        else
            a = atil; b = btil;
        end
        L = b - a;
    end
    minimizador = (a + b)/2;
    r = return([a, b])
endfunction

