%% Extract curve from SVG files
%https://www.flaticon.com/packs/birds-silhouette
%SVG = loadsvg("bird.svg",0.1,false);%building.svg
imgnm = "BDayLove";
SVG = loadsvg(strcat(imgnm,".svg"),0.2,false);
%imread("bird.svg");
figure;hold on;set(gca,"YDir","reverse");axis image equal
for i = 1:length(SVG)
plot(SVG{i}(:,1),SVG{i}(:,2))
fprintf("%d\n",length(SVG{i}))
%pause
end
%%
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
%
Xcoor = Xseq;
Ycoor = Yseq;
% Xcoor = SVG{1}(:,1);
% Ycoor = SVG{1}(:,2);
% Xcoor = B{1}(:,2);
% Ycoor = B{1}(:,1);
Zcoor = Xcoor + j * Ycoor;
N = length(Xcoor);
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
Save = true;
t_list = [0:1/N:0.08, ...
          floor(0.08*N)/N:2/N:0.25, ...
          floor(0.25*N)/N:4/N:0.66, ...
          floor(0.66*N)/N:8/N:1];
t_full = 0:1/N:1;
i_list = floor(t_list * N) + 1;
%t_list = [0:1/N:1];
comp_n = 800;
ar = 16/9;
window_W_fin = 400; window_W_init = 5; % too small is not good looking
%window_H_fin = 150; window_H_init = 5; % window_W;
window_H_fin = window_W_fin / ar; window_H_init = window_W_init / ar; % window_W;

winW_list = camera_zoom_curv(t_list, window_W_init, window_W_fin, 'staircase', 11);%  logspace(log10(window_W_init), log10(window_W_fin), length(t_list) + 1);
winH_list = camera_zoom_curv(t_list, window_H_init, window_H_fin, 'staircase', 11);%  logspace(log10(window_H_init), log10(window_H_fin), length(t_list) + 1);
if Save
    v = VideoWriter(sprintf('%s_FFT_zoom_anime',imgnm),'Motion JPEG AVI');%'Motion JPEG AVI''MPEG-4'
    v.FrameRate = 15;open(v);
end
N = length(UAmp);
idxes = reshape([2:comp_n+1; 
                 N+1-[1:comp_n]],1,[]); % select the index of component, interlace the components, make the f, N-f pair neighbor
% Precompute the coordinates
z_displ_mat = Ucoef(idxes).*exp(1i*2*pi*freq(idxes)*t_full )/N; % compute complex displacement vector
z_list_mat = cumsum([zeros(1, length(t_full)); z_displ_mat],1); % cum sum the displacement vectors
z_store = z_list_mat(end,:);
rads = UAmp(idxes)/N; % Amplitude for each cycle
%
%figure("Position",[200,300,1200,600]);%axis equal;set(gca,"YDir","reverse")
figure("Position",[2100          110        1280         720]);%axis equal;set(gca,"YDir","reverse")
for idx = 1:length(t_list)-3 %0001:0300%length(t_list)-5%1:length(t_list)-8%1:length(t_list)
i = i_list(idx);
cla;hold on;axis equal;set(gca,"YDir","reverse");xticks([]);yticks([]);%axis off;
set(gca,'Color','k','position',[0 0 1 1],'units','normalized')
z_list = z_list_mat(:, i); zcur = z_list_mat(end, i); z_traj = z_store(1, 1:i);
% plot the epi circle
drawCircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads])
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
xlim(real(zcur) + [-winW_list(idx), winW_list(idx)]);
ylim(imag(zcur) + [-winH_list(idx), winH_list(idx)]);
drawnow
%Frms{fi}=getframe(gcf);fi=fi+1;
if Save, writeVideo(v,getframe(gcf));end
%drawnow
end
if Save, close(v);end