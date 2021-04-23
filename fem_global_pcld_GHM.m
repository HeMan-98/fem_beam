%%
%--Paramters--%
l = 0.5;
n = 50;
le = l/n;
E1 = 70.3e9;
E3 = 70.3e9;
G2 = 7.9e4;
al1 = 1.047;
w1 = 4943.06;
z1 = 3911.89;
h1 = 0.0045;
h2 = 0.0002;
h3 = 0.0005;
A1 = h1*0.038;
A2 = h2*0.038;
A3 = h3*0.038;
I1 = 0.038*h1^3/12;
I3 = 0.038*h3^3/12;
rho1 = 2750;
rho2 = 1099.5;
rho3 = 2750;
%%
%--FEM Matrix Formulations--%
syms x;
P1 = [1 0;1 le];
N12 = [1 x]/P1;
P2 = [1 0 0 0;0 1 0 0;1 le le^2 le^3;0 1 2*le 3*le^2];
N34 = [1 x x^2 x^3]/P2;
N34x = diff(N34,x);
N1 = [N12(1); 0; 0; 0; N12(2); 0; 0; 0];
N1x = diff(N1,x);
Ke1 = double(E1*A1*int(N1x*N1x.',x,[0 le]));
N2 = [0; N12(1); 0; 0; 0; N12(2); 0; 0];
N2x = diff(N2,x);
Ke3 = double(E3*A3*int(N2x*N2x.',x,[0 le]));
N3 = [0; 0; N34(1); N34(2); 0; 0; N34(3); N34(4)];
N3xx = diff(diff(N3,x),x);
Kb1 = double(E1*I1*int(N3xx*N3xx.',x,[0 le]));
Kb3 = double(E3*I3*int(N3xx*N3xx.',x,[0 le]));
N4 = [0; 0; N34x(1); N34x(2); 0; 0; N34x(3); N34x(4)];
N6 = 1/h2*((N2-N1)+(h1/2+h3/2+h2)*N4);
N5 = 1/2*((N1+N2)+((h3-h1)/2)*N4);
Me1 = double(rho1*A1*int(N1*N1.',x,[0 le]));
Me2 = double(rho2*A2*int(N5*N5.',x,[0 le]));
Me3 = double(rho3*A3*int(N2*N2.',x,[0 le]));
Mb1 = double(rho1*A1*int(N3*N3.',x,[0 le]));
Mb2 = double(rho2*A2*int(N3*N3.',x,[0 le]));
Mb3 = double(rho3*A3*int(N3*N3.',x,[0 le]));
%Element matrices%
Kv = double(A2*int(N6*N6.',x,[0 le]));
Ke = Ke1+Ke3+Kb1+Kb3;
Me = Me1+Me2+Me3+Mb1+Mb2+Mb3;
%Assembly%
ndf = 4;
KE = zeros(4+n*ndf);
KV = zeros(4+n*ndf);
M = zeros(4+n*ndf);
for i = 1:n
    for j = 1:ndf+4
        for k = 1:ndf+4
            KE(ndf*(i-1)+j,ndf*(i-1)+k) = KE(ndf*(i-1)+j,ndf*(i-1)+k) + Ke(j,k);
            KV(ndf*(i-1)+j,ndf*(i-1)+k) = KV(ndf*(i-1)+j,ndf*(i-1)+k) + Kv(j,k);
            M(ndf*(i-1)+j,ndf*(i-1)+k) = M(ndf*(i-1)+j,ndf*(i-1)+k) + Me(j,k);
        end
    end
end
%GHM%
[R,S] =eigs(KV,101);
S = G2*S;
R = R*S;
zer = zeros(length(M),length(S));
zero = zeros(length(M));
K = [KE+G2*(1+al1)*KV -al1*R;
    -al1*R.' al1*S];
D = [zero zer
    zer.' (2*al1*z1/w1)*S];
M = [M zer;
    zer.' (al1/w1^2)*S];
%BCs%
K([1,2,3,4],:) = [];
K(:,[1,2,3,4]) = [];
D([1,2,3,4],:) = [];
D(:,[1,2,3,4]) = [];
M([1,2,3,4],:) = [];
M(:,[1,2,3,4]) = [];
%%
%--Analysis--%
[Vectors,Values]=eig(K,M);
wn=sqrt(diag(Values));
Freq=wn/(2*pi);
F(ndf*n+101,1) = 0;
F(ndf*n-1,1) = 1;
w = wn(2,1);
S = complex(K-M*w^2,-D*w);
Z = S\F;
%figure;
plot(abs(Z(ndf-1:ndf:ndf*n-1,1)).*sign(real(Z(ndf-1:ndf:ndf*n-1,1))),'red','lineWidth',1);