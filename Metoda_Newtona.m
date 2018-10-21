function [ X ] = Metoda_Newtona( f1, f2, f3, J, x0, iMax, eps )
% Funkcja rozwi�zuj�ca uk�ad 3 r�wna� nieliniowych postaci f1(x,y,z)=0,
%  f2(x,y,z)=0, f3(x,y,z)=0 metod� Newtona.

% Dane wej�ciowe:
%  f1, f2, f3 - funkcje zmiennych x, y, z
%  J - jakobian np. J= @(x, y,z)[2*x 2*y 2*z; 2*x 2*y 0; 1 0 1];
%  x0 - przybli�enie pocz�tkowe
%  iMax - maksymalna liczba iteracji
%  eps - dok�adno�� kolejnych przyblizen
% Dane wyj�ciowe:
%  Wektor X=[x, y, z] b�d�cy rozwi�zanie uk�adu r�wna�.

X = x0; % przybli�enie pocz�tkowe
pX = X; % poprzednie przybli�enie (na pocz�tku r�wne przybli�eniu pocz�tkowemu)
F = zeros(3, 1);
i=0;
%iter=zeros(iMax, 1);
while ( (norm(X-pX) > eps && i<iMax) || i==0)
    i=i+1;
    iter(i)=norm(X-pX);
    F(1)=f1(X(1), X(2), X(3));
    F(2)=f2(X(1), X(2), X(3));
    F(3)=f3(X(1), X(2), X(3));

    dF = J(X(1), X(2), X(3));

    if ( det(dF)==0 ) 
        disp('podaj inny punkt startowy');
        X = zeros(1, 3);
        break;
    else
        dF = -inv( dF );
        pX = X; 
        X = X + (dF * F)';
    end
end
iter'
if ( i==iMax ) 
        disp('Przekroczono maksymalna liczbe iteracji');
        X = zeros(1, 3);
end
        
end