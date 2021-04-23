for o=1:50
    %--Paramters--%
    l = 0.5;
    n = 50;
    le = l/n;
    E1 = 70.3e9;
    G2 = 7.9e4;
    al1 = 1.047;
    w1 = 4943.06;
    z1 = 3911.89;
    h1 = 0.0050;
    h2 = 0.0002;
    A1 = h1*0.038;
    A2 = h2*0.038;
    I1 = 0.038*h1^3/12;
    rho1 = 2750;
    rho2 = 1099.5;
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
    N3 = [0; 0; N34(1); N34(2); 0; 0; N34(3); N34(4)];
    N3xx = diff(diff(N3,x),x);
    Kb1 = double(E1*I1*int(N3xx*N3xx.',x,[0 le]));
    N4 = [0; 0; N34x(1); N34x(2); 0; 0; N34x(3); N34x(4)];
    N6 = 1/h2*(2*(N2-N1)+(h1+h2)*N4);
    Me1 = double(rho1*A1*int(N1*N1.',x,[0 le]));
    Me2 = double(rho2*A2*int(N2*N2.',x,[0 le]));
    Mb1 = double(rho1*A1*int(N3*N3.',x,[0 le]));
    Mb2 = double(rho2*A2*int(N3*N3.',x,[0 le]));
    %Element matrices%
    Kv = double(A2*int(N6*N6.',x,[0 le]));
    Ke = Ke1+Kb1;
    Me = Me1+Me2+Mb1+Mb2;
    %GHM+Rearrangement%
    [R,S] = eigs(Kv,3);
    S = G2*S;
    R = R*S;
    zer = zeros(length(Me),length(S));
    zero = zeros(length(Me));
    Men = [Me(1:4,1:4) zer(1:4,:) Me(1:4,5:8);
        zer(1:4,:).' (al1/w1^2)*S zer(5:8,:).';
        Me(5:8,1:4) zer(5:8,:) Me(5:8,5:8)];
    Den = [zero(1:4,1:4) zer(1:4,:) zero(1:4,5:8);
        zer(1:4,:).' (2*al1*z1/w1)*S zer(5:8,:).';
        zero(5:8,1:4) zer(5:8,:) zero(5:8,5:8)];
    Ken = [Ke(1:4,1:4)+G2*(1+al1)*Kv(1:4,1:4) -al1*R(1:4,:) Ke(1:4,5:8)+G2*(1+al1)*Kv(1:4,5:8);
        -al1*R(1:4,:).' al1*S -al1*R(5:8,:).';
        Ke(5:8,1:4)+G2*(1+al1)*Kv(5:8,1:4) -al1*R(5:8,:) Ke(5:8,5:8)+G2*(1+al1)*Kv(5:8,5:8)];
    %Assembly%
    nn=o;
    nm=n-nn;
    ll=nn*le;
    ndf = 4 + length(S);
    K = zeros(4+nn*ndf+2*nm);
    D = zeros(4+nn*ndf+2*nm);
    M = zeros(4+nn*ndf+2*nm);
    for i = 1:nn
        for j = 1:ndf+4
            for k = 1:ndf+4
                K(ndf*(i-1)+j,ndf*(i-1)+k) = K(ndf*(i-1)+j,ndf*(i-1)+k) + Ken(j,k);
                D(ndf*(i-1)+j,ndf*(i-1)+k) = D(ndf*(i-1)+j,ndf*(i-1)+k) + Den(j,k);
                M(ndf*(i-1)+j,ndf*(i-1)+k) = M(ndf*(i-1)+j,ndf*(i-1)+k) + Men(j,k);
            end
        end
    end
    Me=Me1+Mb1;
    Ke=Ke1+Kb1;
    Ke([1,2,5,6],:) = [];
    Ke(:,[1,2,5,6]) = [];
    Me([1,2,5,6],:) = [];
    Me(:,[1,2,5,6]) = [];
    for i = 1:nm
        for j = 1:4
            for k = 1:4
                K(ndf*nn+2+2*(i-1)+j,ndf*nn+2+2*(i-1)+k) = K(ndf*nn+2+2*(i-1)+j,ndf*nn+2+2*(i-1)+k) + Ke(j,k);
                M(ndf*nn+2+2*(i-1)+j,ndf*nn+2+2*(i-1)+k) = M(ndf*nn+2+2*(i-1)+j,ndf*nn+2+2*(i-1)+k) + Me(j,k);
            end
        end
    end
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
    F(nn*ndf+2*nm,1) = 0;
    F(nn*ndf+2*nm-1,1) = 1;
    w = wn(3,1);
    S = complex(K-M*w^2,-D*w);
    Z = S\F;
    B(o)=max(abs(Z([ndf-1:ndf:ndf*nn-1,ndf*nn+2-1:2:ndf*nn+2*nm-1],1)));
    hh(o)=ll;
    clearvars -except B hh;
end
plot(hh(10:50),B(10:50),'black','lineWidth',1)