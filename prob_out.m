function P = prob_out(H,t1,t2,rho)

[Ix,Iy,Iz]=prodop(1/2,3);

J = 1;
t = 1/2/J;
Iz1 = Iz(:,:,1);
Iz2 = Iz(:,:,2);
Iz3 = Iz(:,:,3);


%UT = expm(-1i*2*pi*J*Iz1*Iz2*t);

U1 = expm(-1i*H*t1);
U2 = expm(-1i*H*t2);

UZ1 = expm(-1i*pi/2*Iz(:,:,1));
UZ2 = expm(-1i*pi/2*Iz(:,:,2));
UZ3 = expm(-1i*pi/2*Iz(:,:,3));


UY1 = expm(-1i*pi/2*Iy(:,:,1));
UY2 = expm(-1i*pi/2*Iy(:,:,2));
UY3 = expm(-1i*pi/2*Iy(:,:,3));

UX1 = expm(-1i*pi*Ix(:,:,1));
UX2 = expm(-1i*pi*Ix(:,:,2));
UX3 = expm(-1i*pi*Ix(:,:,3));

Cnot12 = UY2'*UZ2'*UZ1...
    *UX3'*UT(H,t/4)*UX2*UX1*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UX1*UT(H,t/4)...
    *UY2;

H = Iz1Iz2
for t = t1
    Cnot12...
    *UY1'*UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)...
    *UY1;

for t = t2-t1
    UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)

P12 = diag(U2'*U1*Cnot12*U1'*rho*U1*Cnot12'*U1'*U2);
%(U2'*U1*Cnot*U1'*rho*U1*Cnot'*U1'*U2);

Cnot23 = UY2'*UZ2'*UZ1...
    *UX1'*UT(H,t/4)*UX2*UX3*UT(H,t/4)*UX1*UT(H,t/4)*UX2*UX3*UT(H,t/4)...
    *UY2;


H = Iz2Iz3
for t = t2
    Cnot23...
   *UY1'*UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)...
    *UY1;

for t = t3-t2
    UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)

P23 = diag(U2'*U1*Cnot23*U1'*rho*U1*Cnot23'*U1'*U2)

Cnot13 = UY2'*minusUZ2*UZ1...
    *UX2*UT(H,t/4)*UX1*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX1*UX3*UT(H,t/4)...
    *UY2;

H = Iz1Iz3
for t = t3
    Cnot13...
    *UY1'*UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)...
    *UY1;

for t = t3-t1
    UX2'*UT(H,t/4)*UX3*UT(H,t/4)*UX2*UT(H,t/4)*UX3*UT(H,t/4)

P13 = diag(U2'*U1*Cnot13*U1'*rho*U1*Cnot13'*U1'*U2)

end