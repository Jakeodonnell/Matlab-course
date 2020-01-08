%   dfun.m
function df = dfunB(a,t)


dfa1 = t./(1+a(2)*t);
dfa2 = -(a(1)*t.^2)./((a(2)*t + 1).^2);
df = [dfa1 dfa2];