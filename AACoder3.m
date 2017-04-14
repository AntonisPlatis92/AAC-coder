function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
%Coding of Level 3
[y,fs]=audioread(fNameIn);
hop=1024;
index=1;

lastFrameSamples=mod(size(y,1),1024);
size(y);
y=padarray(y,[1024 0],'both');
y=padarray(y,[1024-lastFrameSamples 0],'post');
size(y);
totalFrames=ceil(size(y,1)/1024)-1;
winType='SIN';
allFramesT=zeros(2048,2,totalFrames);

for i=1:totalFrames
    frameT=y(index:index+2047,:);
    allFramesT(:,:,i)=frameT;
    if i==1
        prevFrameType='OLS';
        nextFrameT=y(index+hop:index+hop+2047,:);
        frameType=SSC(frameT,nextFrameT,prevFrameType,i);
        prevFrameType=frameType;
        index=index+hop;
    elseif i<totalFrames
        nextFrameT=y(index+hop:index+hop+2047,:);
        frameType=SSC(frameT,nextFrameT,prevFrameType,i);
        index=index+hop;
        prevFrameType=frameType;
    else
        frameType='OLS';
    end

    frameF=filterbank(frameT,frameType,winType);
    
    AACSeq3(i,1).frameType=frameType;
    AACSeq3(i,1).winType=winType;
    if (strcmp(frameType,'ESH'))
        AACSeq1(i,1).chl.frameF=frameF(:,1,:);
        AACSeq1(i,1).chr.frameF=frameF(:,2,:);
    else
        AACSeq1(i,1).chl.frameF=frameF(:,1);
        AACSeq1(i,1).chr.frameF=frameF(:,2);
    end

end


for i=1:totalFrames
    if (strcmp(AACSeq3(i).frameType,'ESH'))
        frameFout=zeros(128,2,8);
        a=zeros(4,2,8);
        [frameFout(:,1,:),a(:,1,:)] = TNS(AACSeq1(i).chl.frameF, AACSeq3(i).frameType);
        [frameFout(:,2,:),a(:,2,:)] = TNS(AACSeq1(i).chr.frameF, AACSeq3(i).frameType);
        AACSeq3(i,1).chl.TNScoeffs=a(:,1,:);
        AACSeq3(i,1).chr.TNScoeffs=a(:,2,:);
        AACSeq2(i,1).chl.frameF=frameFout(:,1,:);
        AACSeq2(i,1).chr.frameF=frameFout(:,2,:);
    else
        frameFout=zeros(1024,2);
        a=zeros(4,2);
        [frameFout(:,1),a(:,1)] = TNS(AACSeq1(i).chl.frameF, AACSeq3(i).frameType);
        [frameFout(:,2),a(:,2)] = TNS(AACSeq1(i).chr.frameF, AACSeq3(i).frameType);
        AACSeq3(i,1).chl.TNScoeffs=a(:,1);
        AACSeq3(i,1).chr.TNScoeffs=a(:,2);
        AACSeq2(i,1).chl.frameF=frameFout(:,1);
        AACSeq2(i,1).chr.frameF=frameFout(:,2);
    end
end

%%%   PART 3   %%%%
for k=1:totalFrames
    if (k==1)
        SMR1=psycho(allFramesT(:,1,k),AACSeq3(k).frameType,zeros(2048,1),zeros(2048,1));
        SMR2=psycho(allFramesT(:,2,k),AACSeq3(k).frameType,zeros(2048,1),zeros(2048,1));
    elseif (k==2)
        SMR1=psycho(allFramesT(:,1,k),AACSeq3(k).frameType,allFramesT(:,1,k-1),zeros(2048,1));
        SMR2=psycho(allFramesT(:,2,k),AACSeq3(k).frameType,allFramesT(:,2,k-1),zeros(2048,1));
    else
        SMR1=psycho(allFramesT(:,1,k),AACSeq3(k).frameType,allFramesT(:,1,k-1),allFramesT(:,1,k-2));
        SMR2=psycho(allFramesT(:,2,k),AACSeq3(k).frameType,allFramesT(:,2,k-1),allFramesT(:,2,k-2));
    end

    [S, sfc, G] = AACquantizer(AACSeq2(k).chl.frameF, AACSeq3(k).frameType, SMR1);
    if (strcmp(AACSeq3(k).frameType,'ESH'))
        stream=cell(1,8);
        codebook=zeros(1,8);
        encodeSfc=cell(1,8);
        for i=1:8
            [stream{i} ,codebook(1,i) ]=encodeHuff(S(:,i) ,loadLUT());
            encodeSfc{i}=encodeHuff(sfc(:,i),loadLUT(),12);
        end
        AACSeq3(k,1).chl.stream=stream;
        AACSeq3(k,1).chl.codebook=codebook;
        AACSeq3(k,1).chl.sfc=encodeSfc;
        AACSeq3(k,1).chl.G=G;
    else
        
        [stream,codebook]=encodeHuff(S,loadLUT());
        AACSeq3(k,1).chl.stream=stream;
        AACSeq3(k,1).chl.codebook=codebook;
        AACSeq3(k,1).chl.sfc=encodeHuff(sfc,loadLUT(),12);
        AACSeq3(k,1).chl.G=G;
    end

    [S, sfc, G] = AACquantizer(AACSeq2(k).chr.frameF, AACSeq3(k).frameType, SMR2);
    
    if (strcmp(AACSeq3(k).frameType,'ESH'))
        stream=cell(1,8);
        codebook=zeros(1,8);
        encodeSfc=cell(1,8);
        for i=1:8
            [stream{i} ,codebook(1,i) ]=encodeHuff(S(:,i) ,loadLUT());
            encodeSfc{i}=encodeHuff(sfc(:,i),loadLUT(),12);
        end
        AACSeq3(k,1).chr.stream=stream;
        AACSeq3(k,1).chr.codebook=codebook;
        AACSeq3(k,1).chr.sfc=encodeSfc;
        AACSeq3(k,1).chr.G=G;
    else
        
        [stream,codebook]=encodeHuff(S,loadLUT());
        AACSeq3(k,1).chr.stream=stream;
        AACSeq3(k,1).chr.codebook=codebook;
        AACSeq3(k,1).chr.sfc=encodeHuff(sfc,loadLUT(),12);
        AACSeq3(k,1).chr.G=G;
    end
    
end  
save(fnameAACoded,'AACSeq3');
end

