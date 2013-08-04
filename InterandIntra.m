%intercluster and intracluster ratio [FinalVar FinalCon]

function [Cont FinalVar] = InterandIntra(ReconsIm, movd_seg, nbSegments, sizeIm,p1,k)

%similarity measure is Euclidean distance between
%average of the cluster and each pixel 
Maskdum = zeros(sizeIm,sizeIm);   
MaskSeg(1:nbSegments) = ... 
        struct('ses',zeros(sizeIm,sizeIm,'double'),'seg',zeros(sizeIm,sizeIm,'double'),'x',zeros(sizeIm,1,'double'), ...
                'y',zeros(sizeIm,1,'double'), 'AvgDis',zeros(sizeIm,1,'double'),'contrast',zeros(1,1,'double'), ...
                'var',zeros(1,1,'double'),'variance',zeros(1,1,'double'),'mean',zeros(1,1,'double'));
NeighSeg(1:nbSegments)= ...
        struct('cdata',zeros(sizeIm,sizeIm,'double'),'NoSeg',zeros(sizeIm,1,'double'),'distoutseg',zeros(1,1,'double'),'Len',zeros(1,1,'double'));
    segm = 1;
    for j=1:nbSegments      % Calculates the nonzero indices of each segment
                            % Segements that moved, dm is the count
         [r,c] = find(movd_seg(p1).cdata == j);
         if ( sum(r)~=0)
             AvgofSeg = mean(diag((ReconsIm(k).cdata(r,c)),0));
             distAvgInClus = zeros(size(r,1),1);
             for p=1:size(r,1)                
                 Maskdum(r(p,1),c(p,1)) = ReconsIm(k).cdata(r(p,1),c(p,1)); 
                 MaskSeg(j).seg(r(p,1),c(p,1)) = segm;
                 distAvgInClus(p,1) = sqrt(abs(Maskdum(r(p,1),c(p,1)).^2- AvgofSeg.^2));
             end
             MaskSeg(j).ses = Maskdum;
             MaskSeg(j).seg = edge(MaskSeg(j).seg,.01);
             MaskSeg(j).x = r;                    
             MaskSeg(j).y = c;
             MaskSeg(j).AvgDis = distAvgInClus;
             segm = segm + 1; 
             Maskdum = zeros(sizeIm,sizeIm);
         end
    end

    for j=1:nbSegments
        u = 1;
         [r,c] = find(MaskSeg(j).seg == 1);       
         for l=1:nbSegments   
                InSec =(MaskSeg(l).seg + MaskSeg(j).seg);
                 if( j~=l && logical(sum(find(InSec == 2))))
                     NeighSeg(j).cdata = NeighSeg(j).cdata + MaskSeg(l).ses;
                     NeighSeg(j).NoSeg(u,1) = l;
                     u = u + 1;
                 end
                 NeighSeg(j).Len = u-1;
                 %NeighSeg(j).NoSeg = NeighSeg(j).NoSeg(1:u,1);
             InSec = zeros(sizeIm,sizeIm);
         end
    end
    

%Ratio of intercluster and intracluster similarity measure
 
for j=1:nbSegments
    
     l = NeighSeg(j).Len;                              % Selecting the main segment
     r = MaskSeg(j).x;
     c = MaskSeg(j).y;
     if(sum(r)~=0)
        RatioInOut = zeros(1,l);
        distAvgInClus = MaskSeg(j).AvgDis;
        [m n] = find(MaskSeg(j).seg == 1);

    %     for d=1:size(m,1)
    %         pm = m(d,1);
    %         pn = m(d,1)-1;
    %        
    %         if (pm > sizeIm)
    %             pm = sizeIm;
    %         end
    %         if (pn>sizeIm)               %Boundary Conditions
    %             pn = sizeIm;
    %         end
    %         if (pm <= 0)
    %            pm = 1;
    %         end
    %         if(pn <= 0)
    %            pn = 1;
    %         end
    % 
    %         MaskSeg(j).contrast = MaskSeg(j).contrast + abs(movdr(1).cdata(pm,n(d,1))-movdr(1).cdata(pn,n(d,1)));
    %     end
    %     
    %     if size(m,1)~=0
    %         MaskSeg(j).contrast = MaskSeg(j).contrast./size(m,1);
    %     end
    %     MaskSeg(j).var = mean(distAvgInClus);
        MaskSeg(j).variance = var(diag((ReconsIm(k).cdata(r,c)),0));
        MaskSeg(j).mean = mean(diag((ReconsIm(k).cdata(r,c)),0));
        if(l~=0)

            for p=1:l
                Contrast = 0;
                idx = NeighSeg(j).NoSeg(p,1);                  %Selecting the 1st Neighbour and so on
                v = MaskSeg(idx).x;
                w = MaskSeg(idx).y;
                MeanOutClus = mean(diag((ReconsIm(k).cdata(v,w)),0));
                VarOutClus = var(diag(ReconsIm(k).cdata(v,w),0));   
                Num = sqrt((MaskSeg(j).mean - MeanOutClus).^2+(MaskSeg(j).variance - VarOutClus).^2);
                den  = sqrt((MaskSeg(j).mean + MeanOutClus).^2+(MaskSeg(j).variance + VarOutClus).^2);
                if (den ~= 0)
                    Contrast = Num./den;
                end

                RatioInOut(1,p) = Contrast;       %choose average if max is too much precise 
            end

    %       Ratio between var of current and neighbour to the differnce between averages                
    %     
         NeighSeg(j).distoutseg = max(RatioInOut);
        end
    end
end

SortedVar = zeros(nbSegments,1);
% % SortedCon = zeros(nbSegments,1);
% % 
for i=1:nbSegments
    SortedVar(i,1) = MaskSeg(i).variance;
% %     SortedCon(i,1) = MaskSeg(i).contrast
end
% % 
% % SortedVar = sort(SortedVar,'ascend');
% % SortedCon = sort(SortedCon,'descend');
% % 
% % percenV = floor(.7*nbSegments);
FinalVar = sum(SortedVar)./nbSegments;
% % percenC = floor(.5*nbSegments);
% % FinalCon = sum(SortedCon(1:percenC))./percenC;

% maxR = NeighSeg(1).distoutseg;
% % for j=1:nbSegments
% % 
% %     if (maxR < NeighSeg(j).distoutseg)
% %         maxR = NeighSeg(j).distoutseg;
% %     end
% % end
meanRatio = 0;
for i =1:nbSegments
    meanRatio = meanRatio + (NeighSeg(i).distoutseg);
end

Cont = meanRatio./nbSegments;
