function Te = texp
SPEEDOFLIGHT = 299792458;
global HEIGHT PIXELSIZE;
pix = 1:16;
[a,b] = ind2sub([4,4],pix);
m = zeros(16);
for mu = 1:16
    for mb = 1:16
        m(mu,mb) = sqrt((a(mu)-a(mb))^2 + (b(mu)-b(mb))^2);
    end
end
Te = sqrt((m*PIXELSIZE).^2 +HEIGHT^2)./SPEEDOFLIGHT .*1e9;
save('texp.mat','Te')
end