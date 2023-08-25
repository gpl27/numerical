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


