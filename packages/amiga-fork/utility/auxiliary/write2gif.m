
function write2gif(handle, step, args)
    arguments
        handle 
        step (1,1) {mustBeInteger(step)}
        args.FileName string {mustEndWithGif} = 'mygif.gif';
        args.DelayTime = 0.1;
    end

    frame = getframe(handle);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    
    if step == 1
        imwrite(imind, cm, args.FileName, 'gif', 'DelayTime', args.DelayTime, 'Loopcount',inf);
    else
        imwrite(imind, cm, args.FileName, 'gif', 'DelayTime', args.DelayTime, 'WriteMode','append');
    end
end

function mustEndWithGif(s)
    if ~endsWith(s, ".gif", 'IgnoreCase', true)
        error("Input must be a string ending with '.gif'.");
    end
end
