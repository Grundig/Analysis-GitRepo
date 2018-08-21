function direction = hitdirection

pix = 1:16;
[i,j] = ind2sub([4,4],pix);
m = zeros(16);
for mu = 1:16
    for mb = 1:16
        m(mu,mb) = atan((i(mu)-i(mb))/(j(mu)-j(mb)));
    end
end

direction = (rad2deg(m)+40)

end