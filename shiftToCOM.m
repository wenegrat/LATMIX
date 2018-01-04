function out = shiftToCOM(input, dye, x)


thetai = x./max(x).*2*pi; % Convert to angle
xi = cos(thetai);
zi = sin(thetai);
nt = length(dye); % Dye should already be integrated vertically

for i=1:nt
    xii = xi.*dye(:,i).';
   zii = zi.*dye(:,i).';
   tbar = atan2(-mean(zii), -mean(xii)) + pi;
   ZCM(i) = max(x).*tbar./(2*pi);
   ind = find(x>=ZCM(i), 1, 'first');
   indshift = ind -length(x)./2;
   out(:,i) = circshift(input(:,i), -indshift, 1); % This properly centers the dye (note not applicable in buoyancy space).

end

end