function captured_frames = makevideo(y, t, lbeam, wbeam, img_noise)
% function makevideo fabricates a video of the motion of the beam

% INPUTS:
% y:                displacement of the beam
% t:                time vector
% lbeam and wbeam:  dimensions of the beam, used to adjust the window
% img_noise:        flag for adding noise ( 1: noise is added, 0: no noise )

% OUTPUTS:
% captured_frames:  synthetic video of the beam motion

% F. T. Colombo, September 2024

ny = size(y,2);

x = 0:lbeam/ny:lbeam;
fs = 1/(t(2)-t(1));
yt0 = ceil(t(1)*fs);
loops = length(t);

ax = 2;
f1 = figure;
if ~yt0
    plot(x, [0 y(1,:)], 'LineWidth',15, 'Color','k')
    ax_vert = ax*max(abs(y(:)));
else
    plot(x, [0 y(yt0,:)], 'LineWidth',15, 'Color','k')
    ax_vert = ax*max(max(abs(y(yt0:end,:))));
end

axis equal
axis([x(1) x(end) -ax_vert ax_vert])
ax = gca;
ax.NextPlot = 'replaceChildren';
ax.Visible = 'off';

M = getframe;
Mx = size(M.cdata,1);
My = size(M.cdata,2);

captured_frames = zeros(Mx,My,loops,'uint8');

lw = ceil((My-1)*wbeam/lbeam);

close(f1)
    
figure,

axis equal
axis([x(1) x(end) -ax_vert ax_vert])
ax = gca;
ax.NextPlot = 'replaceChildren';
ax.Visible = 'off';

for j = 1:loops
    plot(x, [0 y(yt0+j,:)], 'LineWidth',lw, 'Color','k')
    drawnow
    M = getframe;
    
    if ~img_noise
        captured_frames(:,:,j) = rgb2gray(M.cdata);
    else
        captured_frames(:,:,j) = imnoise(rgb2gray(M.cdata),'salt & pepper',0.02);
    end
    
    
end

end