%%
%--Parameters--%
N=50; % Number of elements
rho=2750; % Density
E=70.3e9; % Young Modulus
A=0.038*0.0050; % Cross section area
I=(0.038*0.0050^3)/12; % Second moment of inertia
Le=0.5; % Length
L=Le/N;
%%
% FEM Formulation
syms x
P = [1 0 0 0;0 1 0 0;1 L L^2 L^3;0 1 2*L 3*L^2];
Nn = [1 x x^2 x^3]*inv(P);
Nnxx = diff(diff(Nn,x),x);
Ke = double(E*I*int(Nnxx.'*Nnxx,x,[0 L]));
Me = double(rho*A*int(Nn.'*Nn,x,[0 L]));
% Assembly
K = zeros(2*N+2);
M = zeros(2*N+2);
for i = 1:N
    for j = 1:4
        for k = 1:4
            K(2*(i-1)+j,2*(i-1)+k) = K(2*(i-1)+j,2*(i-1)+k) + Ke(j,k);
            M(2*(i-1)+j,2*(i-1)+k) = M(2*(i-1)+j,2*(i-1)+k) + Me(j,k);
        end
    end
end
% Boundary conditions (cantilever beam)
K([1,2],:)=[]; 
K(:,[1,2])=[];
M([1,2],:)=[];
M(:,[1,2])=[];
%%
% Eigenvectores and eigenvalues
[V,D]=eig(K,M);

omega=sqrt(D);
f=omega/(2*pi);
% finite element and exact solution

ft(1)=((1.8751^2)/(2*pi*Le^2))*sqrt((E*I)/(rho*A)); % principal mode
ft(2)=((4.6941^2)/(2*pi*Le^2))*sqrt((E*I)/(rho*A));
ft(3)=((7.8548^2)/(2*pi*Le^2))*sqrt((E*I)/(rho*A));
modes(1)=f(1,1);
modes(2)=f(2,2);
modes(3)=f(3,3);
modes;
ft;

% Vibration modes (first, second and third modes)

for i=1:(size(V,1)/2)
    X(i,1)=V(2*i-1,1);
    X(i,2)=V(2*i-1,2);
    X(i,3)=V(2*i-1,3);
end

%figure;
%hold on;
%plot(L*(1:N),X(:,1),'red','lineWidth',3);
%plot(L*(1:N),X(:,2),'blue','lineWidth',3);
%plot(L*(1:N),X(:,3),'green','lineWidth',3);
%legend('1st Mode','2nd Mode','3rd Mode');
%Beam Displacement%
F(2*N,1) = 0;
F(2*N-1,1) = 1;
w = omega(2,2);
Z = (K-M*w^2)\F;
%figure;
plot(L*(1:N),Z(1:2:2*N-1,1),'blue','lineWidth',1);
%FRF%
%{
for i = 1:400
    if i <=100
        w(i) = omega(1,1)*i/100;
    elseif i<=200
        w(i) = omega(1,1)+(omega(2,2)-omega(1,1))*(i-100)/100;
    else
        w(i) = omega(2,2)+(omega(3,3)-omega(2,2))*(i-200)/100;
    end
    Z = (K-M*w(i)^2)\F;
    B(i,1) = log10(max(abs(Z(1:2:2*N-1,1))));
end
%figure;
plot(w,B(:,1),'lineWidth',1);
%}