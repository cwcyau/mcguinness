function p = lognormpdf(x, m, s)
%
% calculate the probability density function for a normal distribution
%
% by christopher yau (may 2006)
%

% dim = size(s, 1);
% 
% if dim == 1
%     
    p = -0.5*((x - m)./s).^2 - 0.5*log(2*pi*s.^2);
%     
% else
%     
%     p = -0.5*( (x - m)'*inv(s)*(x - m) ) - 0.5*dim*log(2*pi) - 0.5*log(det(s));
%     
% end
