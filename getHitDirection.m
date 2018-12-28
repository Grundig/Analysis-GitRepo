function hitDirection = getHitDirection
global HEIGHT PIXELSIZE DETECTOR_NORTH;

pix = 1:16;
[i,j] = ind2sub([4,4],pix);
hitDirection = struct();
hitDirection.direction = zeros(16);
hitDirection.inclination = zeros(16);

for mu = 1:16
    for mb = 1:16
        x = i(mu)-i(mb);
        y = j(mu)-j(mb);
        if mu == mb
            %for "vertical" muons we give no direction
            hitDirection.direction(mu,mb) = NaN;
        else
            %NOTE: detector pixel grid is left-handed by design (same as compass angle)
            hitDirection.direction(mu,mb) = rad2deg(atan2(y,x));
        end
        % angle of attack measured from normal to earth surface
        hitDirection.inclination(mu,mb) = rad2deg(atan(sqrt(x^2 + y^2)*PIXELSIZE)/HEIGHT);
    end
end


% correct for the real detector position
% and translate direction range from -180:180 to 0:360
hitDirection.direction = mod(hitDirection.direction - DETECTOR_NORTH, 360); 

