function[roznice, bladmax] = Interpolacja_Hermitea(x, f, df, m)
%Interpolacja funkcjami sklejanymi kubicznymi hermitowskimi z przestrzeni
% S3(Dn, 1). Tablicowanie roznic s'(x)-f'(x) w m rownoodleglych punktach.
% Obliczenie b³êdu maksymalnego w tych punktach.

%Dane wejsciowe:
% x - wektor wezlow interpolacji
% f - funkcja interpolowana
% df - funkcja pochodna funkcji interpolowanej
% m -ilosc rownoodleglych punktow, w ktorych obliczamy wartosci funkcji
% interpolujacej

%Dane wyjsciowe
% roznice - tablica roznic s'(x)=f'(x), gdzie pierwsza kolumna to x, a
% druga to roznica s'(x)-f'(x)
% bladmax - maksymalny blad w tych punktach

n = length(x); %ilosc wezlow
a=x(1); % a - poczatek przedzia³u, b - koniec
b=x(n);
xi = linspace(a, b, m); %tworzymy wektor punktow w ktorych bedziemy szukac wartosci funkcji interpolujacej
x(n+1:n+2) = 0; %dla uproszczenia dalszych obliczen powiekszamy wektor x
x = x([end:end 1:end-1]);
bladmax = 0;

for j=1:length(xi)
    s(j)=0;  %s(x)
    sp(j)=0; %s'(x)
    for i=2:n+1  %petla obliczajaca s(x) oraz s'(x)
       if (xi(j)>=x(i-1) && xi(j)<=x(i) && xi(j)>=a && xi(j)<=b && i~=2)
           sum1=(((xi(j)-x(i-1))*(xi(j)-x(i-1))*(-2*xi(j)+3*x(i)-x(i-1))) ./ ((x(i)-x(i-1)).*(x(i)-x(i-1)).*(x(i)-x(i-1))));
           sum2=(((xi(j)-x(i-1))*(xi(j)-x(i-1))*(xi(j)-x(i))) ./ ((x(i)-x(i-1)).*(x(i)-x(i-1))));
           
           sum3=((2*(xi(j)-x(i-1))*(-2*xi(j)+3*x(i)-x(i-1))-(2*((xi(j)-x(i-1))*(xi(j)-x(i-1))))) ./ ((x(i)-x(i-1)).*(x(i)-x(i-1)).*(x(i)-x(i-1))));
           sum4=(((xi(j)-x(i-1))*(xi(j)-x(i-1))+2*(xi(j)-x(i-1))*(xi(j)-x(i))) ./ ((x(i)-x(i-1)).*(x(i)-x(i-1))));
       elseif (xi(j)>=x(i) && xi(j)<=x(i+1)&& xi(j)>=a && xi(j)<=b  && i~=n+1)
           sum1=(((xi(j)-x(i+1))*(xi(j)-x(i+1))*(-2*xi(j)+3*x(i)-x(i+1))) ./ ((x(i)-x(i+1)).*(x(i)-x(i+1)).*(x(i)-x(i+1))));
           sum2=(((xi(j)-x(i+1))*(xi(j)-x(i+1))*(xi(j)-x(i))) ./ ((x(i+1)-x(i)).*(x(i+1)-x(i))));
           
           sum3=((2*(xi(j)-x(i+1))*(-2*xi(j)+3*x(i)-x(i+1))-(2*((xi(j)-x(i+1))*(xi(j)-x(i+1)))) ) ./ ((x(i)-x(i+1)).*(x(i)-x(i+1)).*(x(i)-x(i+1))));
           sum4=(((xi(j)-x(i+1))*(xi(j)-x(i+1))+2*(xi(j)-x(i+1))*(xi(j)-x(i))) ./ ((x(i)-x(i+1)).*(x(i)-x(i+1))));
       else
           sum1=0;
           sum2=0;
           sum3=0;
           sum4=0;
       end
       s(j) = s(j) + (f(x(i)).*sum1) + (df(x(i)).*sum2);
       sp(j) = sp(j) + (f(x(i))*sum3) + (df(x(i))*sum4);
       
       y(i) = f(x(i)); %wektor wartosci funkcji w wezlach
    end
    fp(j) = df(xi(j)); %f'(x)
    
    if ( abs( sp(j)-fp(j) ) > bladmax )
        bladmax = abs( sp(j)-fp(j) );
    end
    
end
roznice(:, 1)=xi';
roznice(:, 2)=abs(fp-sp)';


%rysowanie wykresow funkcji interpolujacej i interpolowanej
fplot(f,[a b], 'r')
hold on
plot(xi, s, 'b+')
hold on
plot(xi, s)
hold on
x(1)=[];
x(n+1)=[];
y(1)=[];
plot(x, y, 'r*')
hold on
legend('funkcja interpolowana f(x)', 'm rownoodleglych punktow', 'funkcja interpolujaca s(x)', 'wezly funkcji interpolowanej' );

end