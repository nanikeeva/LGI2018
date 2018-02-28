function aP = aprob_out(H,t1,t2,rho)

[Ix,Iy,Iz]=prodop(1/2,2);

J = 1;
t = 1/2/J;
Iz1 = Iz(:,:,1);
Iz2 = Iz(:,:,2);

U1 = expm(-1i*H*t1);
U2 = expm(-1i*H*t2);

UZ1 = expm(-1i*pi/2*Iz(:,:,1));
minusUZ2 = expm(1i*pi/2*Iz(:,:,2));
UY2 = expm(-1i*pi/2*Iy(:,:,2));
minusUY2 = expm(1i*pi/2*Iy(:,:,2));

UT = expm(-1i*2*pi*J*Iz1*Iz2*t);

Cnot = minusUY2*UZ1*minusUZ2*UT*UY2;

X =kron([0 1; 1 0],eye(2));
aCnot = X*Cnot*X;
aP = diag(U2'*U1*aCnot*U1'*rho*U1*aCnot'*U1'*U2);