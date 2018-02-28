%% FileName : generate_free_evolH

%% Description
%
% This program calculates the Free Evolution hamiltonian or the Hamiltonian
% in absence of Radio Frequency Pulse

function [H_free_evo,D_hilbert] = generate_free_evolH(spinlist,spinNumbers,v,J,useother)

% Ix=gra.Ix; Iy=gra.Iy; Iz=gra.Iz;
[Ix,Iy,Iz,~,~,~,sIHz] = prodop(spinNumbers,spinlist);

%% Generating Zeeman Hamiltonian
% We are generating $Ho = -\sum_i 2\pi v_i I_z^i$
D_hilbert =1;
for spi=1:length(spinNumbers)
    D_hilbert=D_hilbert*(2*spinNumbers(spi)+1)^spinlist(spi);
end

Ho=zeros(D_hilbert);
for k=1:length(v)
    H_off = -2*pi*v(k)*Iz(:,:,k);
    Ho=Ho+H_off;
end

%% Generating J-coupling Hamiltonian
% We are generating $Hj = \sum_{ij}^{j>i} 2\pi J (I_x^i I_x^j + I_y^i I_y^j + I_z^i I_z^j)$
Hj=zeros(D_hilbert);

for j=1:length(spinlist)
    spinlistsum(j)=sum(spinlist(1:j));
end

for k=1:sum(spinlist)
    for n=1:sum(spinlist)
        H_coup = 2*pi*J(k,n)*(Iz(:,:,k)*Iz(:,:,n));
        Hj=Hj+H_coup;
    end
end
if useother==1
for l=1:length(spinlist)
    if(l==1)
        for k=1:spinlist(l)
            for n=1:spinlist(l)
                H_coup = 2*pi*J(k,n)*(Ix(:,:,k)*Ix(:,:,n) + Iy(:,:,k)*Iy(:,:,n));
                Hj=Hj+H_coup;
            end
        end
    else
        for k=spinlist(l-1)+1:spinlistsum(l)
            for n=spinlist(l-1)+1:spinlistsum(l)
                H_coup = 2*pi*J(k,n)*(Ix(:,:,k)*Ix(:,:,n) + Iy(:,:,k)*Iy(:,:,n));
                Hj=Hj+H_coup;
            end
        end
    end
end
end
%% Total Internal Hamiltonian
H_free_evo=Ho+Hj;

