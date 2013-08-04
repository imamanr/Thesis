function MotionMat = CalcMotionMatrix ( DiffOut, dm )  
thrsed = .15;
MotionMat(1:dm) = ...
        struct('xd1',zeros(2,1,'double'),'yd1', zeros(1,2,'double'), ... 
               'xm1',zeros(1,1,'double'),'ym1', zeros(1,1,'double'), ...
               'xd2',zeros(2,1,'double'),'yd2', zeros(1,2,'double'), ... 
               'xm2',zeros(1,1,'double'),'ym2', zeros(1,1,'double'), ...
               'xd3',zeros(2,1,'double'),'yd3', zeros(1,2,'double'), ... 
               'xm3',zeros(1,1,'double'),'ym3', zeros(1,1,'double'), ...
               'xd4',zeros(2,1,'double'),'yd4', zeros(1,2,'double'), ... 
               'xm4',zeros(1,1,'double'),'ym4', zeros(1,1,'double'), ...
               'xd5',zeros(2,1,'double'),'yd5', zeros(1,2,'double'), ... 
               'xm5',zeros(1,1,'double'),'ym5', zeros(1,1,'double'), ... 
               'xd6',zeros(2,1,'double'),'yd6', zeros(1,2,'double'), ... 
               'xm6',zeros(1,1,'double'),'ym6', zeros(1,1,'double'), ... 
               'xd7',zeros(2,1,'double'),'yd7', zeros(1,2,'double'), ... 
               'xm7',zeros(1,1,'double'),'ym7', zeros(1,1,'double'), ... 
               'xd8',zeros(2,1,'double'),'yd8', zeros(1,2,'double'), ... 
               'xm8',zeros(1,1,'double'),'ym8', zeros(1,1,'double'));
     np = 1;
     for  pt = 1:dm

          for ns = 1:8

              if(DiffOut(np,1)> thrsed)
                  nsc = int2str(ns);
                  temx = strcat('xd',nsc);
                  temxm = strcat('xm',nsc);
                  temy = strcat('yd',nsc);
                  temym = strcat('ym',nsc);
                  if (ns == 1 || ns == 4 || ns == 5 || ns == 8)
                     yInc = 1;%DiffOut(np,1)*2;
                     xInc = 1;%DiffOut(np,1)*1;
                  end
                  if (ns == 2 || ns == 3 || ns == 6 || ns == 7)
                     yInc = 1;%DiffOut((ns+pt-1),1)*1;
                     xInc = 1;%DiffOut((ns+pt-1),1)*2;
                  end

                  if ( ns == 1 || ns == 8)
                      PosY = 0;
                      if (ns == 1) ,PosX = 1; end
                      if (ns == 8), PosX = 0; end
                  end
                  if ( ns == 2 || ns == 3)
                      PosX = 1;
                      if (ns == 2) ,PosY = 0; end
                      if (ns == 3), PosY = 1; end
                  end 
                  if ( ns == 4 || ns == 5)
                      PosY = 1;
                      if (ns == 4) ,PosX = 1; end
                      if (ns == 5), PosX = 0; end
                  end
                  if ( ns == 6 || ns == 7)
                      PosX = 0;
                      if (ns == 6) ,PosY = 1; end
                      if (ns == 7), PosY = 0; end
                  end

                  MotionMat = setfield(MotionMat,{1,pt},temy,[1-PosY PosY]);   %Determines wether neg or pos direction on yaxis
                  MotionMat = setfield(MotionMat,{1,pt},temx,[1-PosX PosX]);   %Determines wether neg or pos direction on xaxis
                  MotionMat = setfield(MotionMat,{1,pt},temxm,xInc);
                  MotionMat = setfield(MotionMat,{1,pt},temym,yInc);
              end
              np = np + 1;
          end    
     end
