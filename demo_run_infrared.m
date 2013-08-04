clc
clc
clear all
pathname = pwd;
filename = strcat(pathname,'/foreman.avi');

BndryThres = .23;
global RowThres nbSegments sizeIm
RowThres = .01;
segm = 1;
dm = 1;
sizeIm = 64;
p1=0;    
k=1;%1;
Re_seg = 1;
DetChangeSeg = 0;
DB_ratio = .1;
VarA = [];
ContA = [];
nbSegments = 5;
b=1;
NoiseFloor = 0;
SelectSeg(1:8) = struct( 'cdata',zeros(sizeIm,sizeIm,'double'));
SegMov = zeros(nbSegments,1);
MaskforMeas = zeros(sizeIm);
H = fspecial('gaussian',[2 2],1);
load('irdp07')

%%%%%%%%%%%%%%%%%%%%%%Incldue for Infrared Videos%%%%%%%%%%%%%%%%
nFrames = size(mov,2);                                         %
movdr(1:nFrames) = ...                                        %
        struct('cdata', zeros(sizeIm, sizeIm, 1, 'double'),...%
               'colormap', []);                               %
k = 1;                                                        %
for i=25:nFrames   % start at 250 for infrared video          %
    movdr(k).cdata = mov(i).cdata./255;                       %
    k = k + 1;                                                %
end                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fix_Seq(1:6) = ...
        struct('data',zeros(sizeIm, sizeIm, 1, 'uint8'));
ReconsIm(1:nFrames) = ... 
        struct('cdata', zeros(sizeIm, sizeIm, 1, 'double'),...
               'colormap', []);
MeasMatrix(1:nFrames) = ... 
        struct('cdata', zeros(sizeIm, sizeIm, 1, 'uint8'),'meas',zeros(1,1,1,'double'), ...
               'colormap', []);           

CorrMotn(1:nFrames-200) = ...
    struct('Self',zeros(sizeIm,sizeIm,1,'double'),'Crx', zeros(sizeIm,sizeIm,1,'double'));
movd_seg(1:nFrames) = ...
    struct('cdata',zeros(sizeIm,sizeIm,1,'double'));


for i=1:6
    Fix_Seq(i).data = sqrt(6).* randn(sizeIm,sizeIm);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% First frame reconstruction %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [r,c] = find(movdr(k).cdata > 0.22);
     NoMeas = floor(log10(sizeIm*sizeIm)*size(r,1)*1/.8);%10*percen
     NoMeas = .6*4096;%size(r,1);
     NewMeasMat_ber = zeros(NoMeas,sizeIm*sizeIm);
     NewMeasMat = ones(64,64);
     for j = 1:NoMeas
         Ber = sqrt(6)*randn(sizeIm);
         NewMeasMat_ber(j,:) = (NewMeasMat(:).* Ber(:))';
     end
         NewMeasVector = NewMeasMat_ber * movdr(k).cdata(:);
         
    opts.mu = 2^8;
    opts.beta = 2^5;
    opts.tol = 1E-3;
    opts.maxit = 300;
    opts.TVnorm = 1;
    opts.nonneg = true;

    [U, out] = TVAL3(NewMeasMat_ber,NewMeasVector,sizeIm,sizeIm,opts);
    
    ReconsIm(k).cdata = reshape(U,sizeIm,sizeIm);
    MeasMatrix(k).cdata = NewMeasMat;
    MeasMatrix(k).meas = NoMeas;
    clear NoMeas r c 
    sumt = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction using adaptive Sampling
%Segment the previous frame and run inside loop for all segments
while (k <= nFrames) 
    DiffBoundary = sum(abs(movdr(k).cdata(1,:) - movdr(k+1).cdata(1,:))) + ...
              sum(abs(movdr(k).cdata(:,1) - movdr(k+1).cdata(:,1))) + ...
              sum(abs(movdr(k).cdata(sizeIm,:) - movdr(k+1).cdata(sizeIm,:))) + ...
              sum(abs(movdr(k).cdata(:,sizeIm) - movdr(k+1).cdata(:,sizeIm)));
      
      if (Re_seg == 1)
          p1 = p1+1;
          b = 1;
          VarA = [];
          ContA = [];
          nbSegments = 40;
          [movd_seg(p1).cdata,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(ReconsIm(k).cdata,nbSegments);
          [dm MaskSeg dmV thresh] = MovedSegmentManual(movd_seg,movdr,movdr, k,sizeIm,nbSegments,p1);  %%*added ReconsIm
          dm = dm - 1;
      end
     
    if ( DiffBoundary < BndryThres)         
        c = 0;
        dmMov = [];
        if (DetChangeSeg == 1)
            [dc MaskSegCh dmCh threshCh] = MovedSegmentManual(movd_seg,movdr,movdr, k,sizeIm,nbSegments,p1);
            for i=1:size(dmCh,2)
                if(sum(dmCh(1,i)~=dmV)==size(dmV,2))
                    c = c + 1;
                   dmMov(1,c) = dmCh(1,i);
                end
            end
        end                                                

        CombineMask (1:dm+c) = struct('ses',zeros(sizeIm,sizeIm,'double'), ... 
                                    'x',zeros(sizeIm,1,'double'),'y',zeros(sizeIm,1,'double'));                
        if (c~=0)
            for i=1:c
                imx = find(dmMov(1,i) == dmCh);
                CombineMask(1,i).ses = MaskSegCh(1,imx).ses;
                CombineMask(1,i).x = MaskSegCh(1,imx).x;
                CombineMask(1,i).y = MaskSegCh(1,imx).y;
            end
        end
        
        for i=c+1:dm+c                    
            CombineMask(i).ses = MaskSeg(i-c).ses;
            CombineMask(i).x = MaskSeg(i-c).x;
            CombineMask(i).y = MaskSeg(i-c).y;            
        end
        dmV = [dmMov dmV];
        dm = size(dmV,2);
        
        [MeasMaskAvg MeasOutAvg MaskSegAvg] = MaskOutnIn(CombineMask,dm,sizeIm); %Calculates In and outside mask for each segment


        MeasAvgSeg = MeasMaskAvg * ReconsIm(k).cdata(:);
        MeasAvgSeg1 = MeasMaskAvg * movdr(k+1).cdata(:);

        MeasOutSeg = MeasOutAvg * ReconsIm(k).cdata(:);
        MeasOutSeg1 = MeasOutAvg * movdr(k+1).cdata(:);

        DiffOut = abs(MeasOutSeg1 - MeasOutSeg);
        k1=1;
        for u = 0:8:(dm*8)-8
            DiffOut(k1:k1+7,1)=DiffOut(k1:k1+7,1)./sum(DiffOut(k1:k1+7,1));
            k1 = k1 + 8;
        end
        MotionMat = CalcMotionMatrix ( DiffOut, dm ); % Calculate the motion matrix
        
        %%%%%%%%%%%%%%%%%%%%%% NEW MEASUREMENT MASK%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [NewMeasMat MaskSeg] = DispMeasureMaskNew( MotionMat, CombineMask, MaskSegAvg, movdr, sizeIm,dm,MeasOutSeg1,k);
             
    else if ( DiffBoundary > BndryThres)

        [dmB BorderMask] = MovedBorderSeg(movdr, movdr, movd_seg,sizeIm,nbSegments,k,RowThres,p1,thresh);%%*added ReconsIm

        % calculate the motion vector of these segments Calcualte the measurement matrix for these segments
        % Make measurements for this area on the outer side                                                  
            
        dmbS = size(dmB,2)-1;
        c1=0;
        dmBor = [];
        for i=1:size(dmB,2)
            if(sum(dmB(1,i)~=dmV)==size(dmV,2))
               c1 = c1 + 1;
               dmBor(1,c1) = dmB(1,i);
            end
        end
            
        dm = c1 + size(dmV,2);
        dmV = [dmBor dmV];
        
        c2 = 0;
        dmMov = [];

        if (DetChangeSeg == 1)
            [dc MaskSegCh dmCh threshCh] = MovedSegmentManual(movd_seg,movdr,movdr, k,sizeIm,nbSegments,p1);

            for i=1:size(dmCh,2)
                if(sum(dmCh(1,i)~=dmV)==size(dmV,2))
                    c2 = c2 + 1;
                   dmMov(1,c2) = dmCh(1,i);
                end
            end
        end 
        
        dmV = [dmMov dmV];
        dm = size(dmV,2);
        
        
        CombineMask (1:dm) = struct('ses',zeros(sizeIm,sizeIm,'double'), ... 
                                    'x',zeros(sizeIm,1,'double'),'y',zeros(sizeIm,1,'double'));
       if (c2~=0)
            for i=1:c2
                imx = find(dmMov(1,i) == dmCh);
                CombineMask(1,i).ses = MaskSegCh(1,imx).ses;
                CombineMask(1,i).x = MaskSegCh(1,imx).x;
                CombineMask(1,i).y = MaskSegCh(1,imx).y;
            end
       end
        
       if(c1~=0)
           
            for i=1+c2:c1+c2
                ibx = find(dmB==dmBor(1,i-c2));
                CombineMask(i).ses = BorderMask(ibx).ses;
                CombineMask(i).x = BorderMask(ibx).x;
                CombineMask(i).y = BorderMask(ibx).y;
            end
       end
        for i=c1+c2+1:dm
            
            CombineMask(i).ses = MaskSeg(i-(c1+c2)).ses;
            CombineMask(i).x = MaskSeg(i-(c1+c2)).x;
            CombineMask(i).y = MaskSeg(i-(c1+c2)).y;
            
        end
       
        [MeasMaskAvg MeasOutAvg MaskSegAvg] = MaskOutnIn(CombineMask,dm,sizeIm); %Calculates In and outside mask for each segment

        %FinalMeasMat = 
        MeasAvgSeg = MeasMaskAvg * ReconsIm(k).cdata(:);
        MeasAvgSeg1 = MeasMaskAvg * movdr(k+1).cdata(:);

        MeasOutSeg = MeasOutAvg * ReconsIm(k).cdata(:);
        MeasOutSeg1 = MeasOutAvg * movdr(k+1).cdata(:);

        DiffOut = abs(MeasOutSeg1 - MeasOutSeg);        
        k1=1;
        for u = 0:8:(dm*8)-8
            DiffOut(k1:k1+7,1)=DiffOut(k1:k1+7,1)./sum(DiffOut(k1:k1+7,1));
            k1 = k1 + 8;
        end
        MotionMat = CalcMotionMatrix ( DiffOut, dm ); % Calculate the motion matrix
       
        %%%%%%%%%%%%%%%%%%%%%% NEW MEASUREMENT MASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [NewMeasMat MaskSeg] = DispMeasureMaskNew( MotionMat, CombineMask, MaskSegAvg, ReconsIm, sizeIm,dm,MeasOutSeg1,k);
        end
    end
    
        
%% Reconstruction from adaptivly acquired samples

if(sum(sum(isnan(NewMeasMat)==0)))
   U = zeros(sizeIm,sizeIm);
   NoMeas = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if(sum(sum(NewMeasMat)~=0))         
         N = log10(sum(sum(NewMeasMat)));
         MaxIm = max(max(ReconsIm(k).cdata));
         deltaVal = 0.22;
         [r c] = find(ReconsIm(k).cdata.*NewMeasMat > (MaxIm - deltaVal));
         sparsty = size(r,1);
         C = 1.5;
         Err = .8;
         Nt = sum(sum(NewMeasMat));
         NoMeas = .4 * Nt;

         upBound = 0.5;
         LwBound = 0.3;
         if (NoMeas > upBound*sum(sum(NewMeasMat)))
             NoMeas = upBound*sum(sum(NewMeasMat));
         end
         if (NoMeas < LwBound*sum(sum(NewMeasMat)))
             NoMeas = LwBound*sum(sum(NewMeasMat));
         end

         NewMeasMat_ber = zeros(NoMeas,sizeIm*sizeIm);

         for j = 1:NoMeas
             Ber = sqrt(6)*randn(sizeIm);
             NewMeasMat_ber(j,:) = (NewMeasMat(:).* Ber(:))';
         end
             Noise_var = .01;
             regp = 4.5;
             s2 = sizeIm*sizeIm;
             NewMeasVector = NewMeasMat_ber * movdr(k+1).cdata(:);% + Noise_var.*randn(NoMeas,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RECONSTRUCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        opts.mu = 2^8;
        opts.beta = 2^5;
        opts.tol = 1E-3;
        opts.maxit = 300;
        opts.TVnorm = 1;
        opts.nonneg = true;
        tic 

        [U, out] = TVAL3(NewMeasMat_ber,NewMeasVector,sizeIm,sizeIm,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RECONSTRUCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        t = toc;
        sumt = t + sumt;
         U = U.*NewMeasMat;
     end
     
    temp =  ReconsIm(k).cdata;

    [r,c] = find(NewMeasMat ~=0 );
    for i=1:size(r,1)
        temp(r(i,1),c(i,1)) = U(r(i,1),c(i,1));
    end
    
    NoiseLevel = sum(sum(sqrt(abs(temp.^2 - movdr(k+1).cdata.^2))))./(sizeIm*sizeIm);
    DetChangeSeg = 1;
    PercenCoverThres = .7; % .1 for simulations
    if (sum(sum(NewMeasMat)) < PercenCoverThres*sizeIm*sizeIm)
        Re_seg = 0;
    else
        Re_seg = 1;
    end
        
    ReconsIm(k+1).cdata=temp;
    MeasMatrix(k+1).cdata = NewMeasMat;
    MeasMatrix(k+1).meas = NoMeas;
    k = k + 1;
    display(k)
    figure(2)
    colormap(gray)
    title 'Reconstructed Image using 36% Measurements'
    colorbar
    imagesc(temp)
    close all
    clear NewMeasMat U out temp NewMeasMat_ber NewMeasVector Ber MotionMat DiffOut ... 
        MeasAvgSeg MeasAvgSeg1 MeasOutSeg MeasOutSeg1 MeasMaskAvg MeasOutAvg MaskSegAvg ...
        BorderMask CombineMask VarA ContA Cont Var NoMeas
        
 end
 
end
%% Video generation

Videocmp(1:nFrames) = ...
struct('cdata', zeros(sizeIm,2* sizeIm, 3, 'uint8'),'colormap',[]);
 for i=1:nFrames
     comparison(i).cdata(:,1:64)= movdr(i).cdata(:,:);
     comparison(i).cdata(:,64+1:2*64)= ReconsIm(i).cdata;
 end
for i=1:nFrames
     Videocmp(i).cdata(:,:,1)=comparison(i).cdata.*255;
     Videocmp(i).cdata(:,:,2)=comparison(i).cdata.*255;
     Videocmp(i).cdata(:,:,3)=comparison(i).cdata.*255;
 end
 for i=1:nFrames
     Videocmp(i).cdata(:,:,1)=comparison(i).cdata.*255;
     Videocmp(i).cdata(:,:,2)=comparison(i).cdata.*255;
     Videocmp(i).cdata(:,:,3)=comparison(i).cdata.*255;
 end
movie2avi(Videocmp,'InsertFramesVid1','compression','None')
