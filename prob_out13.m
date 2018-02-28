function P13 = prob_out13(H,t1,t3,rho)

[Ix,Iy,Iz]=prodop(1/2,3);

J = 1;
t = 1/2/J;
Iz1 = Iz(:,:,1);
Iz2 = Iz(:,:,2);
Iz3 = Iz(:,:,3);


%UT = expm(-1i*2*pi*J*Iz1*Iz2*t);

U1 = expm(-1i*H*t1);
U2 = expm(-1i*H*t3);

UZ1 = expm(-1i*pi/2*Iz(:,:,1));
UZ2 = expm(-1i*pi/2*Iz(:,:,2));
UZ3 = expm(-1i*pi/2*Iz(:,:,3));

UY1 = expm(-1i*pi/2*Iy(:,:,1));
UY2 = expm(-1i*pi/2*Iy(:,:,2));
UY3 = expm(-1i*pi/2*Iy(:,:,3));

UX1 = expm(-1i*pi*Ix(:,:,1));
UX2 = expm(-1i*pi*Ix(:,:,2));
UX3 = expm(-1i*pi*Ix(:,:,3));

Cnot13 = UY2'*UZ2'*UZ1...
    *UX3'*UT(H,t/4)*UX2*UX1*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UX1*UT(H,t/4)...
    *UY2;

H = Iz1*Iz2
t = t3
    OP5 = Cnot13...
    *UY1'*UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)...
    *UY1;

t = t3-t1
    OP6 = UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)*OP5

P13 = diag(OP6*rho*OP6')
end