function w_list = camera_zoom_curv(t_list, winit, wend, mode, varargin)
if nargin == 3
    mode = "exp";
end
switch mode
    case "exp"
        w_list = logspace(log10(winit), log10(wend), length(t_list) + 1);
    case "lin"
        w_list = linspace((winit), (wend), length(t_list) + 1);
    case "staircase"
        lvl_n = 11;
        if length(varargin) == 1, lvl_n = varargin{1}; end
        w_levels = logspace(log10(winit), log10(wend), lvl_n);
        w_list = zeros(1,length(t_list));
        stairL = length(t_list) / lvl_n;
        lmts = floor(((0:lvl_n-1)' + [0.2,0.8]) * stairL);lmts(1,1)=1;lmts(end,2)=length(t_list);
        for l = 1:lvl_n
            w_list(lmts(l,1):lmts(l,2)) = w_levels(l);
            if  l ~= lvl_n
            w_list(lmts(l,2):lmts(l+1,1)) = logspace(log10(w_levels(l)), log10(w_levels(l + 1)), lmts(l+1,1) - lmts(l,2) + 1);
            end
        end
end
end