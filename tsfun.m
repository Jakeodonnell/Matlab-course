function n = tsfun(yrensad)
load tsoutlier.mat
%y = randn(12)
n = length(y);
k = input('Give half window, k ');
figure(1)
plot(y,'b');
hold on
%beräkna median värden ymed i fönster
for i = 1:n
    n1 = max(1,i-k);
    n2 = min(n,i+k);
    ymed(i) = median(y(n1:n2));
    %ymed(i) = mean(y(n1:n2))
    
end
%beräkna MAD(i)
outlier = zeros(1,n);
for i = 1:n
    n1 = max(1,i-k);
    n2 = min(n,i+k);
    mad(i) = 1.486 * median(abs(y(n1:n2)-ymed(i)));
    if (abs(y(i)-ymed(i))) > 3*mad(i)
        outlier(i) = 1;
        
    end
end
    
%plot outlier med rött kryss.
yrensad = y;
yrensad(outlier(:) == 1) = []; %outlier(:) alla element i outlier
for i = 1:n %n:-1:1  
    %börja baklänges, annars tar man bort i och i+1 ersätter i., dvs man kan bara ta bort ett element.
    if outlier(i) == 1
        
        plot(i,y(i),'xr')
        
    end
end
figure(2)
plot(yrensad,'k')
%ylim([-20 25]) % for tsoutlier
title('Outliers removed')



