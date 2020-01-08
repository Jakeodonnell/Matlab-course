%   dfun.m
function df = dfun(a,t)
dfa1 = exp(-a(2)*t); 
dfa2 = -a(1)*t.*exp(-a(2)*t); 
dfa3 = exp(-a(3)*t); 
dfa4 = -a(3)*t.*exp(-a(4)*t);
df = [dfa1 dfa2 dfa3 dfa4]

%dfa1 = t./(1+a(2)*t);
%dfa2 = -(a(1)*t.^2)./((a(2)*t + 1).^2);
%df = [dfa1 dfa2];