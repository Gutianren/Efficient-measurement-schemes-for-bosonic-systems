function [y] = GroundEnergyConverge(dmax, file)
    for d = 1:dmax+1
        y(d) = GetGroundEnergy(d,file);
        fprintf('%f\n',y(d));
    end
%     eg = y(end);
%     y = abs(eg - y);
%     x = 2:dmax;
%     semilogy(x,y(x),x,y(x)/eg);
%     legend('Absolute','Relative');
%     xlabel('d');
%     ylabel('\Delta E_g');
%     title('Convergence of ground energy');
%     grid minor;
end