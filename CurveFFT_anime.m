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
imgnm = "BDayLove";
SVG = loadsvg(strcat(imgnm,".svg"),0.2,false);
%imread("bird.svg");
figure;hold on;set(gca,"YDir","reverse");axis image equal
for i = 1:length(SVG)
plot(SVG{i}(:,1),SVG{i}(:,2))
fprintf("%d\n",length(SVG{i}))
pause
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
N = length(Xcoor)
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
comp_n = 700;
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
window_W = 100;
window_H = 50;%window_W;
comp_n = 800;

N = length(UAmp);
idxes = reshape([2:comp_n+1;N+1-[1:comp_n]],1,[]); % select the index of component
figure("Position",[200,300,1200,600]);%axis equal;set(gca,"YDir","reverse")
if Save
v = VideoWriter(sprintf('%s_FFT_anime',imgnm),'MPEG-4');%'Motion JPEG AVI'
v.FrameRate = 15;open(v);
end
%Frms = {};Fi=1;
z_store = [];
for t = 0:1/N:1
cla;hold on;axis equal;set(gca,"YDir","reverse");xticks([]);yticks([]);%axis off;
set(gca,'Color','k','position',[0 0 1 1],'units','normalized')
z_displ = Ucoef(idxes).*exp(1i*2*pi*t*freq(idxes) )/N; % compute complex displacement vector
z_list = cumsum([0; z_displ]);zcur=z_list(end); % cum sum the 
z_store = [z_store, zcur]; 

% plot the epi circle
rads = UAmp(idxes);
drawCircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N)
% viscircles([real(z_list(1:end-1)),imag(z_list(1:end-1))], [rads]/N,'Color',[0.3,0.3,0.3]);
% that is too slow

% Line plot ot arrow plot 
plot(real(z_list),imag(z_list),"LineWidth",0.75, 'Color', [1,1,1])
%quiver(real(z_list(1:end-1)),imag(z_list(1:end-1)),real(z_displ),imag(z_displ), 0, ...
%    'Color', [1,1,1],'MaxHeadSize',0.05)

% Plot some history
plot(real(z_store),imag(z_store),"LineWidth",2, 'Color', [0.9290, 0.6940, 0.1250])
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
% Goto https://www.youcompress.com/avi/ for video compression