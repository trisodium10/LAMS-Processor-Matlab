function b = smooth_std(a,n)
a = a(:,end:-1:1);
b = zeros(size(a));
for k = 1:size(a,2)
    if k <= n
        b(:,k) = nanstd(a(:,1:k).');
    else
        b(:,k) = nanstd(a(:,(k-n):k).');
    end
end
b = b(:,end:-1:1);