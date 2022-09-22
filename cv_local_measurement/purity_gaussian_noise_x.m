function [p_noise, gs_noise] = purity_gaussian_noise_x(gs, e_a, sigma_a)
    % Adding to the ground state a gaussian shift exp(ipa)
    % returning p_noise denoting the probability and gs_noise(:,i) denoting
    % the noised state
    % The states are written in the energy space
    % for numerics, we set the range of deviation at most 3*sigma
    
    global dx
    
    a_max = 3*sigma_a;
    a = -a_max+ea:dx:a_max+ea;
    pdf_noise = normpdf(a,e_a,sigma_a);
    p_noise = pdf_noise * dx;
    na = floor(a/dx);
    
    
    
    
    
    
    
end