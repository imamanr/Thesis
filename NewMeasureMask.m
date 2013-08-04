function NewMeasMat = NewMeasurMask( MotionMat, MaskSeg, sizeIm)

    NewMask(1:dm*8) = struct('seg', zeros(sizeIm,sizeIm,'double'));
    NewAvgMeas = zeros(7,1); 
    diffAvgNew = zeros(7,1);
    NewMeasMat = zeros(sizeIm);
    NewMaskAvgseg = zeros(sizeIm);
    NewMaskAvgout = zeros(sizeIm);
    NewTempMask(j).seg = zeros(sizeIm,sizeIm);
    multx = 0;
    multy = 0;
    Temp1x = 0;
    Temp1y = 0;
    Temp2x = 0;
    Temp2y = 0;
    NewTempMask(1:nbSegments) = struct('seg',zeros(sizeIm,sizeIm,'double'));
    ns = 1;
    TempM = zeros(sizeIm);
    for j = 1:dm
         for d = 1:8
             ns = ns + 1;
             nsc = int2str(d);
             temx = strcat('xd',nsc);
             temxm = strcat('xm',nsc);
             temy = strcat('yd',nsc);
             temym = strcat('ym',nsc);         
             Temp1y = getfield(MotionMat,{1 j},temy); %#ok<GFLD>
             Temp2y = ceil(1*(getfield(MotionMat,{1 j},temym))); %#ok<GFLD>
             if ( sum(Temp1y) ~= 0)                     
                if (Temp1y(1,2) == 1) 
                    multy = 1;
                else if (Temp1y(1,1) == 1 )
                        multy = -1;
                    end
                end                                    
             end

             Temp1x = ceil(1*(getfield(MotionMat,{1 j},temx)));  %#ok<GFLD>            
             Temp2x = ceil(1*(getfield(MotionMat,{1 j},temxm))); %#ok<GFLD>         

            if ( sum(Temp1x) ~= 0)

                 r = MaskSeg(j).x;
                 c = MaskSeg(j).y;
                 if (Temp1x(1,2) == 1) 
                    multx = 1;
                else if (Temp1x(1,1) == 1)
                        multx = -1;
                    end
                 end
            end 
            t=0;
            inc = 1;
            if (multx == 1 || multx == -1 || multy == 1 || multy == -1)

                while ( inc <= 7)

                  for i = 1:size(r,1)                  
                          x = ceil(r(i,1)+ multx*inc*Temp2x);
                          y = ceil(c(i,1)+ multy*inc*Temp2y);
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
                          NewTempMask(j).seg(x,y) = MaskSeg(j).ses(r(i,1),c(i,1));                  
                  end 

                   NewTempMask(j).seg = NewTempMask(j).seg(1:sizeIm,1:sizeIm);               
                   NewMaskAvgout = MaskSegAvg(ns-1).mask .*(1-NewTempMask(j).seg./max(max(NewTempMask(j).seg)));           
                   NewAvgMeas(inc,1) = sum(sum(NewMaskAvgout .* movdr(k).cdata));
                   diffAvgNew(inc,1) = abs(NewAvgMeas(inc,1) - MeasOutSeg1(ns-1));

                   inc = inc + 1;
                   NewMaskAvgseg = zeros(sizeIm);
                   NewMaskAvgout = zeros(sizeIm);
                   NewTempMask(j).seg = zeros(sizeIm);
                end

                incF = min(find( diffAvgNew==min(diffAvgNew)));                               %#ok<MXFND>

                for i = 1:size(r,1)                  
                    x = ceil(r(i,1)+ multx*incF*Temp2x);
                    y = ceil(c(i,1)+ multy*incF*Temp2y);
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
                    TempM(x,y) = MaskSeg(j).ses(r(i,1),c(i,1));                  
                end
                TempM = TempM(1:sizeIm,1:sizeIm);
                NewMask(j).seg = NewMask(j).seg + TempM;
                NewMask(j).seg = ceil(NewMask(j).seg./max(max(NewMask(j).seg))); 
                TempM = zeros(sizeIm);
                NewTempMask(j).seg = zeros(sizeIm,sizeIm);
                diffAvgNew = zeros(7,1);
                multx = 0;
                multy = 0;
                Temp1x = 0;
                Temp1y = 0;
                Temp2x = 0;
                Temp2y = 0;
            end
         end    
         NewMask(j).seg(isnan(NewMask(j).seg)) = 1;
         %NewMask(j).seg =  NewMask(j).seg./max(max( NewMask(j).seg));
         Temp = NewMask(j).seg(1:sizeIm,1:sizeIm);
         NewMeasMat = NewMeasMat + Temp;
         NewMeasMat = ceil(NewMeasMat./max(max(NewMeasMat)));
         Temp = zeros(sizeIm);     
    end   
