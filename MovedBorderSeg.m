function [dmB BorderMask] = MovedBorderSeg(movdr, ReconsIm, movd_seg,sizeIm,nbSegments,k,RowThres,p1,threh)

        dmB = [];
MaskforMeas = zeros(sizeIm);
BorderMask(1:nbSegments) = struct('ses',zeros(sizeIm,sizeIm,'double'), ... 
                                   'x',zeros(sizeIm,1,'double'),'y',zeros(sizeIm,1,'double'));

SegI(1:2) = ...
    struct('Mean',zeros(nbSegments,1,1,'double'));
     
                Diff_FirstCol = sum(abs(ReconsIm(k).cdata(:,1) - movdr(k+1).cdata(:,1)));
                Diff_LastCol  = sum(abs(ReconsIm(k).cdata(:,sizeIm) - movdr(k+1).cdata(:,sizeIm)));
                Diff_FirstRow = sum(abs(ReconsIm(k).cdata(1,:) - movdr(k+1).cdata(1,:)));
                Diff_LastRow  = sum(abs(ReconsIm(k).cdata(sizeIm,:) - movdr(k+1).cdata(sizeIm,:)));
                mm = 1;
                
                for d = 1:nbSegments
                    [r,c] = find(movd_seg(p1).cdata == d);
                    if(size(r,1)>5)
                        SegI(1).Mean(d,1) = mean(diag((ReconsIm(k).cdata(r,c)),0));
                        SegI(2).Mean(d,1) = mean(diag((movdr(k+1).cdata(r,c)),0));  %Calcaulating averages of next frame
                        tempF_col = find(c==1);
                        tempF_row = find(r==1);
                        tempL_col = find(c==sizeIm);
                        tempL_row = find(r==sizeIm);

                        if (size(tempF_col,1) == 0)
                            C_condf = false;
                        else 
                            C_condf = true;
                        end

                        if (size(tempF_row,1) == 0)
                            R_condf = false;
                        else 
                            R_condf = true;
                        end
                        if (size(tempL_col,1) == 0)
                            C_condl = false;
                        else 
                            C_condl = true;
                        end

                        if (size(tempL_row,1) == 0)
                            R_condl = false;
                        else 
                            R_condl = true;
                        end

                        if ( Diff_FirstCol > RowThres && C_condf && abs(SegI(1).Mean(d,1)-SegI(2).Mean(d,1)) > threh)
                            for p=1:size(r)
                                BorderMask(mm).ses(r(p,1),c(p,1)) = movd_seg(p1).cdata(r(p,1),c(p,1))./d;%Orignal Segmentation
                            end

                            BorderMask(mm).x = r;
                            BorderMask(mm).y = c;                         
                            mm = mm + 1;
                            dmB = [dmB d];

                            else if ( Diff_LastCol > RowThres && C_condl && abs(SegI(1).Mean(d,1)-SegI(2).Mean(d,1)) > threh)
                                     %enable all segments along this Col
                                     for p=1:size(r)
                                         BorderMask(mm).ses(r(p,1),c(p,1)) = movd_seg(p1).cdata(r(p,1),c(p,1))./d;%Orignal Segmentation
                                     end
                                     BorderMask(mm).x = r;
                                     BorderMask(mm).y = c;                                 
                                     mm = mm + 1;
                                     dmB = [dmB d];                                                                                                    

                                 else if ( Diff_FirstRow > RowThres && R_condf && abs(SegI(1).Mean(d,1)-SegI(2).Mean(d,1)) > threh)
                                        %enable all segments along this row               
                                        for p=1:size(r)
                                            BorderMask(mm).ses(r(p,1),c(p,1)) = movd_seg(p1).cdata(r(p,1),c(p,1))./d; %Orignal Segmentation
                                        end

                                        BorderMask(mm).x = r;
                                        BorderMask(mm).y = c;                                     
                                        mm = mm + 1;
                                        dmB = [dmB d];               

                                    else if ( Diff_LastRow > RowThres && R_condl && abs(SegI(1).Mean(d,1)-SegI(2).Mean(d,1)) > threh)
                                    %enable all segments along this row                
                                            for p=1:size(r)
                                                BorderMask(mm).ses(r(p,1),c(p,1)) = movd_seg(p1).cdata(r(p,1),c(p,1))./d;  %Original Segmentation
                                            end

                                            BorderMask(mm).x = r;
                                            BorderMask(mm).y = c; 
                                            mm = mm + 1;
                                            dmB = [dmB d];
                                        end
                                     end
                                end
                        end
                    end
                 end
                
                
                
                % Calculate the average of consolidated segment Mask
                
                
                % Keep segments above the threshold
