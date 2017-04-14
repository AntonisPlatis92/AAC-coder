%Funtion to define the frameType of a Frame
%OLS=ONLY_LONG_SEQUENCE
%LSS=LONG_START_SEQUENCE
%LPS=LONG_STOP_SEQUENCE
%ESH=EIGHT_SHORT_SEQUENCE
function frameType = SSC(frameT, nextFrameT, prevFrameType,i)
if (strcmp(prevFrameType,'LSS'))   %Check if previous is LSS
    frameType='ESH';
    return;
elseif (strcmp(prevFrameType,'LPS')) %Check if previous is LPS
    frameType='OLS';
    return
else      %Prev is ESH or OLS
    b=[0.7548 -0.7548];  %Numurator of filter 
    a=[1 -0.5095]; %Denumenator of the filter
    hop=128;  
    index=1+448+128; %Starting index
    esh_1=false; %Check if we found at least one ESH in channel 1 & 2
    esh_2=false;
    for l=1:8
        filtered_frame_1=filter(b,a,nextFrameT(index:index+127,1)); %Filter the subFrames
        filtered_frame_2=filter(b,a,nextFrameT(index:index+127,2));
        s_1(l)=sumsqr(filtered_frame_1); %Sum s of subFrames
        s_2(l)=sumsqr(filtered_frame_2);
        
        if (l>1) %Calculate the ds
            sumS_1=sum(s_1(1:l-1));
            sumS_2=sum(s_2(1:l-1));
            ds_1(l)=s_1(l)/((1/(l-1))*sumS_1);
            ds_2(l)=s_2(l)/((1/(l-1))*sumS_2);
        else
            ds_1(l)=1;
            ds_2(l)=1;
        end
        if (s_1(l)>10^(-3) && ds_1(l)>10) %if true for at least once per channel, frame is ESH
            esh_1=true;
            if (strcmp(prevFrameType,'OLS'))
                frameType_1='LSS';
                
            else
                frameType_1='ESH';
            end
        end
        if (s_2(l)>10^(-3) && ds_2(l)>10) %Same for 2nd channel
            esh_2=true;
            if (strcmp(prevFrameType,'OLS'))
                frameType_2='LSS';
                
            else
                frameType_2='ESH';
            end
        end
        index=index+hop;
    end
    if (~esh_1) %True if we did not find any ESH subFrame
        if (strcmp(prevFrameType,'OLS'))
             frameType_1='OLS';
        else
             frameType_1='LPS';
        end
    end
    if (~esh_2)
        if (strcmp(prevFrameType,'OLS'))
             frameType_2='OLS';
        else
             frameType_2='LPS';
        end
    end
    %Now according to given table for the two chanels, define the frameType
    %of the whole Framef
    if (strcmp(frameType_1,'OLS') && strcmp(frameType_2,'OLS'))
        frameType='OLS';
        return;
    elseif (strcmp(frameType_1,'OLS') && strcmp(frameType_2,'LSS'))
        frameType='LSS';
        return;
    elseif (strcmp(frameType_1,'OLS') && strcmp(frameType_2,'ESH'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'OLS') && strcmp(frameType_2,'LPS'))
        frameType='LPS';
        return;
    elseif (strcmp(frameType_1,'LSS') && strcmp(frameType_2,'OLS'))
        frameType='LSS';
        return;
    elseif (strcmp(frameType_1,'LSS') && strcmp(frameType_2,'LSS'))
        frameType='LSS';
        return;
    elseif (strcmp(frameType_1,'LSS') && strcmp(frameType_2,'ESH'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'LSS') && strcmp(frameType_2,'LPS'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'ESH') && strcmp(frameType_2,'OLS'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'ESH') && strcmp(frameType_2,'LSS'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'ESH') && strcmp(frameType_2,'ESH'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'ESH') && strcmp(frameType_2,'LPS'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'LPS') && strcmp(frameType_2,'OLS'))
        frameType='LPS';
        return;
    elseif (strcmp(frameType_1,'LPS') && strcmp(frameType_2,'LSS'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'LPS') && strcmp(frameType_2,'ESH'))
        frameType='ESH';
        return;
    elseif (strcmp(frameType_1,'LPS') && strcmp(frameType_2,'LPS'))
        frameType='LPS';
        return;
    end
    
end

end

