function frameFout = iTNS(frameFin, frameType, TNScoeffs)
%Inverse Function for Temporal Noise Shaping
if (strcmp(frameType,'ESH'))
    frameFout=zeros(128,8);
    for i=1:8
        fil=[1;-TNScoeffs(:,i)];
        frameFout(:,i)=filter(1,fil,frameFin(:,i));
    end
else
    fil=[1;-TNScoeffs];
    frameFout=filter(1,fil,frameFin);
end
end

