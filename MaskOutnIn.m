function [MeasMaskAvg MeasOutAvg MaskSegAvg] = MaskOutnIn(MaskSeg,dm,sizeIm)
     
        dem = 1;
        segm = 1;
        skip = 1;
        MaskSegAvg(1:8*dm) = struct('seg',zeros(sizeIm,sizeIm,'double'),'out',zeros(sizeIm,sizeIm,'double'),'mask',zeros(sizeIm,sizeIm,'double'));
        NewMaskAvg(1:8*dm) = struct('seg',zeros(sizeIm,sizeIm,'double'),'out',zeros(sizeIm,sizeIm,'double'));
        l = 0;
        MeasMaskAvg = zeros(dm*8,sizeIm.^2);
        MeasOutAvg = zeros(dm*8,sizeIm.^2);

        for j = 1:dm             %Calculates the masks of directions
                                 %Maks each segments into 8 separte
                                 %directions
%                 if (nargin == 4)
%                     if(size(find( j == dmB),2) ~= 0 ) %zero when not already done = true for condition
%                         skip = 0;                     %one when already done
%                     end
%                 end
%                 
%                 if (skip ~= 0 || nargin == 3)       % Donot skip if skip is 1
                    
                    c = MaskSeg(j).y;
                    r = MaskSeg(j).x;
                 if (size(c,1)~=0 && size(r,1)~=0)
                    min_Xy = [c(min(find(r==min(r)))) min(r)];             %#ok<MXFND>
                    min_Yx = [min(c) r(min(find(c==min(c))))];             %#ok<MXFND>
                    max_Xy = [ c(max(find(r==max(r)))) max(r)];            %#ok<MXFND> % Finding Max and min for Segment i
                    max_Yx = [ max(c) r(max(find(c==max(c))))];            %#ok<MXFND>

                    avgy=0;
                    avgx=0;
                    point=0;
                    segm = max(max(MaskSeg(j).ses));
                    for ig=1:sizeIm
                        for ih=1:sizeIm
                            if (MaskSeg(j).ses(ig,ih) == segm)
                                avgx = avgx + ig;
                                avgy = avgy + ih;
                                point = point + 1;
                            end
                        end
                    end

                    Cen(1,1) = avgx./point;
                    Cen(1,2) = avgy./point;

                    Radius = [sqrt((Cen(1,1)-min_Xy(1,2)).^2+(Cen(1,2)-min_Xy(1,1)).^2)...
                              sqrt((Cen(1,1)-max_Xy(1,2)).^2+(Cen(1,2)-max_Xy(1,1)).^2)...
                              sqrt((Cen(1,1)-min_Yx(1,2)).^2+(Cen(1,2)-min_Yx(1,1)).^2)...
                              sqrt((Cen(1,1)-min_Yx(1,2)).^2+(Cen(1,2)-min_Yx(1,1)).^2)];

                    MaxRad = max(Radius);
                    MaskSeg(j).Rad(1,1) = MaxRad;
                    MaskSeg(j).Rad(1,2) = Cen(1,1);
                    MaskSeg(j).Rad(1,3) = Cen(1,2);

                    th = 0:pi/50:2*pi;
                    xunit = MaxRad * cos(th) + Cen(1,1);
                    yunit = MaxRad * sin(th) + Cen(1,2);

                    HorLinex1 = MaxRad * cos(0) + Cen(1,1);
                    HorLiney1 = MaxRad * sin(0) + Cen(1,2);


                    HorLinex2 = MaxRad * cos(th(1,50)) + Cen(1,1);
                    HorLiney2 = MaxRad * sin(th(1,50)) + Cen(1,2);

                    VerLinex1 = MaxRad * cos(th(1,25)) + Cen(1,1);
                    VerLiney1 = MaxRad * sin(th(1,25)) + Cen(1,2);


                    VerLinex2 = MaxRad * cos(th(1,75)) + Cen(1,1);
                    VerLiney2 = MaxRad * sin(th(1,75)) + Cen(1,2);

                    TiltLine1x1 = MaxRad * cos(th(1,14)) + Cen(1,1);
                    TiltLine1y1 = MaxRad * sin(th(1,14)) + Cen(1,2);             

                    TiltLine2x1 = MaxRad * cos(th(1,39)) + Cen(1,1);
                    TiltLine2y1 = MaxRad * sin(th(1,39)) + Cen(1,2);             

                    TiltLine1x2 = MaxRad * cos(th(1,64)) + Cen(1,1);
                    TiltLine1y2 = MaxRad * sin(th(1,64)) + Cen(1,2);             

                    TiltLine2x2 = MaxRad * cos(th(1,89)) + Cen(1,1);
                    TiltLine2y2 = MaxRad * sin(th(1,89)) + Cen(1,2); 

                    masks(1:8) = struct ( 'seg',zeros(sizeIm,sizeIm,'double'));
                    in=0;
                    fi = pi/4;

                    for f=1:8

                        for th = in:pi/100:fi

                            for ra = 0:.1:MaxRad+2


                                x = floor(ra * cos(th) + MaskSeg(j).Rad(1,2));
                                y = floor(ra * sin(th) + MaskSeg(j).Rad(1,3));
                                if (x>sizeIm)
                                    x=sizeIm;
                                end
                                if (y>sizeIm)
                                    y=sizeIm;
                                end
                                if (x <= 0)
                                    x=1;
                                end
                                if(y <= 0)
                                    y=1;
                                end
                                masks(f).seg(x,y) = 1; %#ok<AGROW>

                            end
                        end
                        in = fi;
                        fi = in + pi/4;
            %           masks(f).seg && MaskSeg(f).Seg
                        MaskSegAvg(f+l).seg = MaskSeg(j).ses.* masks(f).seg; 
                        MaskSegAvg(f+l).out = masks(f).seg.*(1-MaskSeg(j).ses./max(max(MaskSeg(j).ses)));
                        MaskSegAvg(f+l).mask = masks(f).seg;
                    end    
                    mh=1;    

                    for ma = 1:8
                        if (sum(sum(MaskSegAvg(mh).seg)) ~= 0)
                            MeasMaskAvg(mh+l,:) = MaskSegAvg(mh+l).seg(:);
                            MeasOutAvg(mh+l,:) = MaskSegAvg(mh+l).out(:);
                            mh = mh + 1;
                        end
                    end

                    mh=1;
                    l = l + 8;                            
                    dem = dem + 1;
                
                %skip = 1;
                 end
                 
        end

    
