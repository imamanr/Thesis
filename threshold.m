function Value_Thresh = threshold(hist,MAX_HIST)
    
    big_peak = 1;
    start_one = 0;
    smooth = 0;
    delete = 0;
    delete_first = 0;
    do_invert = 0;
    
    pk = findstart(hist,MAX_HIST);
    st = findfirst(hist,MAX_HIST);
    fi = findend(hist,MAX_HIST);
    d1 = pk - st;
    d2 = fi - pk;
    if ((d1 < 0) || (d2 < 0)) 
        printf('Error:histogram have strange peak')
        exit
    end
    
    if (d1>d2)
        do_invert = 1;
    end
    
    if (do_invert) 
        for i = 1:MAX_HIST/2
            tmp = hist(i,1);
            hist(i,1) = hist(MAX_HIST-1-i,1);
            hist(MAX_HIST-1-i,1) = tmp;
        end
    end
    
    pk2 = findsecond(hist,MAX_HIST,pk);
%    ratio = hist(pk,1) / hist(pk2,1);
%     if (ratio > 10) 
%         sprintf('WARNING: ratio of largest to second largest histogram bin %d\n',st);
%         sprintf('++++++++ maybe you should delete the largest histogram bin using the -D option');
%     end

    if (big_peak)
        st = findstart(hist,MAX_HIST);
    elseif (start_one)
        st = 1;       
    else
        st = 1;      
        sprintf('starting from peak at position %d\n',st);
    end
    
    Value_Thresh = findcorner2(hist,st,MAX_HIST);
    
    
    %//* invert threshold back again */
    if (do_invert) 
        Value_Thresh = MAX_HIST - Value_Thresh; 
    end
    
%    /* find largest peak */








