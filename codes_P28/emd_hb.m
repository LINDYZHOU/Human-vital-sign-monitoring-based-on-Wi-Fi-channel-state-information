function [output] = emd_hb(input,no_imf)
% For EMD process
output =[];
for i = 1:size(input,2)       %each subcarrier
    [imf_h,~] = emd(input(:,i),'Interpolation','pchip','display',0);
    %Name-'interpolation', consisting of  either 'spline' or 'pchip'.
    % 'spline', if x is a smooth signal
    % 'pchip', if x is a nonsmooth signal
%     output=[output,imf_h(:,no_imf)+ imf_h(:,no_imf+1)];    
    output=[output,imf_h(:,no_imf)];    
end
end