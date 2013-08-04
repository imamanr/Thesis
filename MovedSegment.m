function [dm MaskSeg dmV Est_thresh] = MovedSegment(movd_seg,movdr,ReconsIm, k,sizeIm,nbSegments,p1)
dm = 1;
dmV = [];
MaskforMeas = zeros(sizeIm);
Est_thresh = 1;
Seg(1) = ...
    struct('Mean',zeros(nbSegments,1,1,'double'),'Var', zeros(nbSegments,1,'double'));

MaskSeg(1:nbSegments) = ... 
    struct('seg',zeros(sizeIm,sizeIm,'double'),'ses',zeros(sizeIm,sizeIm,'double'), ...
            'x',zeros(sizeIm,1,'double'),'y',zeros(sizeIm,1,'double'),'Rad',zeros(1,3,'double'));
 
diff = zeros(nbSegments,1);
MaxHist = nbSegments;

for i=1:nbSegments
    
    [r,c] = find(movd_seg(1,p1).cdata == i);        
    SegFirstFrameMean =   mean(diag((ReconsIm(k).cdata(r,c)),0)); %movdr changed to ReconsIm
    SegSecondframeMean = mean(diag((movdr(k+1).cdata(r,c)),0));   %Averages from net frame
    diff(i,1) = abs(SegFirstFrameMean - SegSecondframeMean);
    
end

%diff = diff./max(diff);

if(sum(diff)~=0)
    
    hist_diff = hist(diff, MaxHist)';
    thresh_est = threshold(hist_diff,MaxHist)+1;

    scale = min(diff):max(diff)/nbSegments:max(diff);
    if (thresh_est>size(scale,2))
        thresh_est = size(scale,2);
    end
    if (thresh_est<=0)
            thresh_est = 2;%2 for rest of videos 
    end


    if (size(scale,2)~=0)
        Est_thresh = scale(1,thresh_est);

        for j=1:nbSegments      % Calculates the Segment that moved. and output is MaskSeg structure having all
                                    % Segements that moved, dm is the count

                [r,c] = find(movd_seg(p1).cdata == j);
                if ( size(r,1)> 5)
                    
                    Seg(k).Mean(j,1) = mean(diag((ReconsIm(k).cdata(r,c)),0));
                    Seg(k).Var(j,1) = var(diag((ReconsIm(k).cdata(r,c)),0));

                    Seg(k+1).Mean(j,1) = mean(diag((movdr(k+1).cdata(r,c)),0));
                    Seg(k+1).Var(j,1)  = mean(diag((movdr(k+1).cdata(r,c)),0));
                    Maskdum = zeros(sizeIm,sizeIm);


                    if (abs(Seg(k).Mean(j,1)-Seg(k+1).Mean(j,1)) > Est_thresh)
                        for p=1:size(r)
                            MaskforMeas(r(p,1),c(p,1)) = j./j;                
                            Maskdum(r(p,1),c(p,1)) = j./j;                              
                        end
                    end

                    if (sum(sum(Maskdum)) ~= 0)
                        MaskSeg(dm).seg = edge(Maskdum,0.01); %ceil(abs(imfilter(
                        MaskSeg(dm).ses = Maskdum;
                        MaskSeg(dm).x = r;                    %should be opposite i.e.    
                        MaskSeg(dm).y = c;                    %MaskSeg(dm).x = r  
                        dm = dm + 1;                          %MaskSeg(dm).y = c    fixed at later in the program after line 307
                        dmV = [dmV j];
                    end

                        clear r c;
                end
        end
    end
end
