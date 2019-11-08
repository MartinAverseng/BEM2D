function[] = animateWave(gridx,gridy,k,amplitude)

nit_per_lambda = 20;
dt = pi/(nit_per_lambda*k);

figure;
vals = amplitude;
l = max(abs(amplitude(:)));
while true
    vals = exp(-1i*k*dt)*vals;
    imagesc(gridx,gridy,real(vals))
    caxis([-l l])
    colormap('gray')
    axis xy;
    axis equal;
    drawnow;
end


end


