addpath imgs
%% Extract curve data from BW image boundary
img = imread("HighNote.jpg");
img_bw = img(:,:,1) < 127;
B=bwboundaries(img_bw,8);%'noholes'
figure
imshow(img)
hold on 
for Bdr_i=1:length(B)
plot(B{Bdr_i}(:,2),B{Bdr_i}(:,1),"LineWidth",3)
end
%% Extract curve from SVG files
%https://www.flaticon.com/packs/birds-silhouette
%SVG = loadsvg("bird.svg",0.1,false);%building.svg
imgnm = "BDay";
SVG = loadsvg(strcat(imgnm,".svg"),0.2,false);
%imread("bird.svg");
figure;hold on;set(gca,"YDir","reverse");axis image equal
for i = 1:length(SVG)
plot(SVG{i}(:,1),SVG{i}(:,2))
fprintf("%d\n",length(SVG{i}))
%pause
end
%% Post processing, connect the SVG separate parts
Xseq = []; Yseq = []; interp_step = 5;%interp_pnts = 3;
for i = 1:length(SVG)
    Xseq = [Xseq; SVG{i}(:,1)];
    Yseq = [Yseq; SVG{i}(:,2)];
    csr  = [Xseq(end), Yseq(end)];
    if i ~= length(SVG)
        csr_new = [SVG{i+1}(1, :)];
    else
        csr_new = [SVG{1}(1, :)];
    end
    interp_pnts = floor(norm(csr_new - csr) / interp_step);
    Xseq = [Xseq; linspace(csr(1),csr_new(1),interp_pnts+2)'];
    Yseq = [Yseq; linspace(csr(2),csr_new(2),interp_pnts+2)'];
end
figure;hold on;set(gca,"YDir","reverse");axis image equal
plot(Xseq,Yseq)
N = length(Xseq)
%%
Xcoor = Xseq;
Ycoor = Yseq;
% Xcoor = SVG{1}(:,1);
% Ycoor = SVG{1}(:,2);
% Xcoor = B{1}(:,2);
% Ycoor = B{1}(:,1);
Zcoor = Xcoor + j * Ycoor;
N = length(Xcoor);

figure
plot(Xcoor,Ycoor)
axis image equal
set(gca,"YDir","reverse")
Ucoef = fft(Zcoor);
UAmp = abs(Ucoef);
UAng = angle(Ucoef);
figure
subplot(211)
plot(UAmp(2:end-1))
subplot(212)
plot(UAng(2:end-1))
%%
comp_n = 600;
FFTidx = [2:comp_n + 1,N-comp_n+1:N];
%1/length(UAmp)
theta = [0:1/N:1];%[0:length(UAmp)-1]/length(UAmp);%[0:0.001:1]%
freq = [0:length(UAmp)-1]';
Xfit = sum(UAmp(FFTidx) .* cos(2*pi*freq(FFTidx) * theta + UAng(FFTidx)), 1)/length(UAmp);
Yfit = sum(UAmp(FFTidx) .* sin(2*pi*freq(FFTidx) * theta + UAng(FFTidx)), 1)/length(UAmp);
figure
plot(Xfit,Yfit)
axis image equal
set(gca,"YDir","reverse")
%%
Save = false;
window_W = 100;
window_H = 50;%window_W;
comp_n = 600;

N = length(UAmp);
idxes = reshape([2:comp_n+1;N+1-[1:comp_n]],1,[]); % select the index of Freq component
figure("Position",[200,300,1200,600]);%axis equal;set(gca,"YDir","reverse")
if Save
v = VideoWriter(sprintf('%s_FFT_anime',imgnm),'MPEG-4');%'Motion JPEG AVI'
v.FrameRate = 20;open(v);
end
%Frms = {};Fi=1;

% Precompute the coordinates
t_list = 0:1/N:1;
z_displ_mat = Ucoef(idxes).*exp(1i*2*pi*freq(idxes)*t_list )/N; % compute complex displacement vector
z_list_mat = cumsum([zeros(1, length(t_list)); z_displ_mat],1); % cum sum the displacement vectors
z_store = z_list_mat(end,:);
rads = UAmp(idxes); % Amplitude for each cycle
for i = 1:length(t_list)
cla;hold on;axis equal;set(gca,"YDir","reverse");xticks([]);yticks([]);%axis off;
set(gca,'Color','k','position',[0 0 1 1],'units','normalized')

z_list = z_list_mat(:, i); zcur = z_list_mat(end, i); z_traj = z_store(1, 1:i);
% plot the epi circle
drawCircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N)
% viscircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N,'Color',[0.3,0.3,0.3]);
% that is too slow

% Line plot ot arrow plot 
plot(real(z_list),imag(z_list),"LineWidth",0.75, 'Color', [1,1,1])
%quiver(real(z_list(1:end-1)),imag(z_list(1:end-1)),real(z_displ),imag(z_displ), 0, ...
%    'Color', [1,1,1],'MaxHeadSize',0.05)

% Plot some history
plot(real(z_traj),imag(z_traj),"LineWidth",2, 'Color', [0.9290, 0.6940, 0.1250])
%plot(Xfit,Yfit,"LineWidth",2,"Color","black");

% Set the field of View!
xlim(real(zcur) + [-window_W, window_W]);
ylim(imag(zcur) + [-window_H, window_H]);
drawnow
%Frms{fi}=getframe(gcf);fi=fi+1;
if Save, writeVideo(v,getframe(gcf));end
%drawnow
end
if Save, close(v);end

%%
%% Allows zooming in and out, change the camera zooming by a predefined curve.
Save = false;
window_W_fin = 1000;window_W_init = 10; % too small is not good looking
window_H_fin = 500; window_H_init = 5; % window_W;
comp_n = 600;

N = length(UAmp);
idxes = reshape([2:comp_n+1;N+1-[1:comp_n]],1,[]); % select the index of component
figure("Position",[200,300,1200,600]);%axis equal;set(gca,"YDir","reverse")
if Save
v = VideoWriter(sprintf('%s_FFT_zoom_anime',imgnm),'MPEG-4');%'Motion JPEG AVI'
v.FrameRate = 20;open(v);
end
%Frms = {};Fi=1;
z_store = [];
t_list = 0:1/N:1;
winW_list = logspace(log10(window_W_init), log10(window_W_fin), length(t_list) + 1);
winH_list = logspace(log10(window_H_init), log10(window_H_fin), length(t_list) + 1);

% Precompute the coordinates
z_displ_mat = Ucoef(idxes).*exp(1i*2*pi*freq(idxes)*t_list )/N; % compute complex displacement vector
z_list_mat = cumsum([zeros(1, length(t_list)); z_displ_mat],1); % cum sum the displacement vectors
z_store = z_list_mat(end,:);
rads = UAmp(idxes); % Amplitude for each cycle

for i = 1:length(t_list)
cla;hold on;axis equal;set(gca,"YDir","reverse");xticks([]);yticks([]);%axis off;
set(gca,'Color','k','position',[0 0 1 1],'units','normalized')
z_list = z_list_mat(:, i); zcur = z_list_mat(end, i); z_traj = z_store(1, 1:i);
% plot the epi circle
drawCircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N)
% viscircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N,'Color',[0.3,0.3,0.3]);
% that is too slow

% Line plot ot arrow plot 
plot(real(z_list),imag(z_list),"LineWidth",0.75, 'Color', [1,1,1])
%quiver(real(z_list(1:end-1)),imag(z_list(1:end-1)),real(z_displ),imag(z_displ), 0, ...
%    'Color', [1,1,1],'MaxHeadSize',0.05)

% Plot some history
plot(real(z_traj),imag(z_traj),"LineWidth",2, 'Color', [0.9290, 0.6940, 0.1250])
%plot(Xfit,Yfit,"LineWidth",2,"Color","black");

% Set the field of View!
xlim(real(zcur) + [-winW_list(i), winW_list(i)]);
ylim(imag(zcur) + [-winH_list(i), winH_list(i)]);
drawnow
%Frms{fi}=getframe(gcf);fi=fi+1;
if Save, writeVideo(v,getframe(gcf));end
%drawnow
end
if Save, close(v);end
% Goto https://www.youcompress.com/ for video compression
%%
function w_list = camera_zoom_curv(t_list, winit, wend, mode)
if nargin == 3
    mode = "exp";
end
switch mode
    case "exp"
        w_list = logspace(log10(winit), log10(wend), length(t_list) + 1);
    case "staircase"
        w_list = logspace(log10(winit), log10(wend), length(t_list) + 1);
end
end