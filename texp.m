function Te = texp
c = 299792458;
h = 0.7;
a = 0.25;
pix = 1:16;
[i,j] = ind2sub([4,4],pix);
m = zeros(16);
for mu = 1:16
    for mb = 1:16
        m(mu,mb) = sqrt((i(mu)-i(mb))^2 + (j(mu)-j(mb))^2);
    end
end
Te = sqrt((m.*a).^2 +h^2)./c .*1e9;
save('texp.mat','Te')
end