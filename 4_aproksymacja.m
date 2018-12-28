function[wyniki, blad] = aproksymacja(f, n, m)
%Tresc: Aproksymacja œredniokwadratowa ci¹g³a w obszarze D = [0,1]×[0,1] w bazie funkcji: 1,?2,?2,?? . Ca³kowanie z³o¿onymi formu³ami trapezów ze wzglêdu na ka¿d¹ zmienn¹. Tablicowanie funkcji, przybli¿enia i b³êdu w ?2 punktach obszaru D oraz obliczenie b³êdu œredniokwadratowego w tych punktach.
%Dane wejsciowe:
% f(x,y) - funkcja dw?ch zmiannych
% n - w n*n punktach przybli?amy blad
% m - ilo?? podprzedzia??w przedzia?u [0,1] calkowania
%Dane wyjsciowe:
% wyniki - macierz zawierajaca n*n wierszy bedacych kolejnymi punktami (x,y) i 3
% kolumny: f(x,y),  f*(x,y),  |f(x,y)-f*(x,y)|
% blad - blad sredniokwadratowy w n*n punktach

%baza:
baza = {@(x,y)(1), @(x,y)(x*x), @(x,y)(y*y), @(x,y)(x*y)};

G = zeros(4,4); %maciez Grama

for i=1:4
   for j=1:i
       G(i, j) = calkowanie_trapezy(m, m, @(x, y)(baza{i}(x,y)*baza{j}(x, y)));
       G(j, i) = G(i, j);
   end
   d(i,1)=calkowanie_trapezy(m, m, @(x, y)(baza{i}(x,y)*f(x, y)));
end

a=G\d;

f1 = @(x, y)(a(1)*baza{1}(x,y)+a(2)*baza{2}(x,y)+a(3)*baza{3}(x,y)+a(4)*baza{4}(x,y));

wyniki = zeros(n*n, 3);
p=linspace(0, 1, n);
q=linspace(0, 1, n);
pom=1;
blad = 0;
for i=1:n
  for j=1:n
       wyniki(pom, 1) = f(p(i), q(j));
       wyniki(pom, 2) = f1(p(i), q(j));
       wyniki(pom, 3) = abs( wyniki(pom, 2) - wyniki(pom, 1));
       pom=pom+1;
       blad = blad + (f(p(i), q(j))-f1(p(i), q(j)))^2;
  end
end
blad = blad/n;
blad = sqrt(blad);

end