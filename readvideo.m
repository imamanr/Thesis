function [movdr nFrames vidHeight vidWidth] = readvideo(filename,sizeIm)

Vid = mmreader ( filename );
nFrames = Vid.NumberOfFrames ;
vidHeight = Vid.Height;
vidWidth = Vid.Width;

mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
movd(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 1, 'double'),...
           'colormap', []);
movdr(1:nFrames) = ...
    struct('cdata', zeros(sizeIm, sizeIm, 1, 'double'),...
           'colormap', []);
       
       
    for k = 1 : 300 
        mov(k).cdata = read(Vid, k);
        movd(k).cdata = im2double(mov(k).cdata(:,:,1));
        movdr(k).cdata = imresize(movd(k).cdata,[sizeIm sizeIm],'bicubic');
    end
