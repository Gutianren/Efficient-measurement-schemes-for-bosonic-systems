function [e,phi] = AnyPhi(d,n,H)
   % We use GHZ here
   phi = zeros(d^n,1);
   %phi(1) = 1;
   %phi(1) = 1/sqrt(2);
   %phi(2) = 1/sqrt(2);
   
   for i = 1:d
       phi((i-1)*(d^n-1)/(d-1)+1) = 1/sqrt(d);
   end

   % the Hamiltonian we use the one in GroundEnergy
   e = conj(phi') * H * phi;
   fprintf('%f',e);
end