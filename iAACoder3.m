function x = iAACoder3(AACSeq3, fNameOut)
%Inverse Coding of Level 1
x=zeros((size(AACSeq3,1)+1)*1024,2);
index=1;

for k=1:size(AACSeq3,1)
    if (strcmp(AACSeq3(k).frameType,'ESH'))
        decFrameF1=zeros(128,8);
        decFrameF2=zeros(128,8);
        decsfc1=zeros(42,8);
        decsfc2=zeros(42,8);
        for i=1:8
            decFrameF1(:,i)=decodeHuff(AACSeq3(k).chl.stream{i},AACSeq3(k).chl.codebook(1,i),loadLUT());
            decFrameF2(:,i)=decodeHuff(AACSeq3(k).chr.stream{i},AACSeq3(k).chr.codebook(1,i),loadLUT());
            decsfc1(:,i)=decodeHuff(AACSeq3(k).chl.sfc{i},12,loadLUT());
            decsfc2(:,i)=decodeHuff(AACSeq3(k).chr.sfc{i},12,loadLUT());
        end
        iQuanFrameF1=iAACquantizer(decFrameF1,decsfc1,AACSeq3(k).chl.G,AACSeq3(k).frameType);
        iQuanFrameF2=iAACquantizer(decFrameF2,decsfc2,AACSeq3(k).chr.G,AACSeq3(k).frameType);
        iTNSframeF1 = iTNS(iQuanFrameF1, AACSeq3(k).frameType, AACSeq3(k).chl.TNScoeffs);
        iTNSframeF2 = iTNS(iQuanFrameF2, AACSeq3(k).frameType, AACSeq3(k).chr.TNScoeffs);
        frameF=zeros(128,2,8);
        frameF(:,1,:)=iTNSframeF1;
        frameF(:,2,:)=iTNSframeF2;
        
    else
        decFrameF1=decodeHuff(AACSeq3(k).chl.stream,AACSeq3(k).chl.codebook,loadLUT());
        decFrameF2=decodeHuff(AACSeq3(k).chr.stream,AACSeq3(k).chr.codebook,loadLUT());
        decsfc1=decodeHuff(AACSeq3(k).chl.sfc,12,loadLUT());
        decsfc2=decodeHuff(AACSeq3(k).chr.sfc,12,loadLUT());
        iQuanFrameF1=iAACquantizer(decFrameF1,decsfc1,AACSeq3(k).chl.G,AACSeq3(k).frameType);
        iQuanFrameF2=iAACquantizer(decFrameF2,decsfc2,AACSeq3(k).chr.G,AACSeq3(k).frameType);
        iTNSframeF1 = iTNS(iQuanFrameF1, AACSeq3(k).frameType, AACSeq3(k).chl.TNScoeffs);
        iTNSframeF2 = iTNS(iQuanFrameF2, AACSeq3(k).frameType, AACSeq3(k).chr.TNScoeffs);
        frameF=zeros(1024,2);
        frameF(:,1)=iTNSframeF1;
        frameF(:,2)=iTNSframeF2;
        
    end

    frameT=iFilterbank(frameF,AACSeq3(k).frameType,AACSeq3(k).winType);
%     if (k==2)
%         AACSeq3(k).frameType
%         frameT(500:600,:)
%     end
    x(index:index+2047,:)=x(index:index+2047,:)+frameT;
    index=index+1024;
    
end
x=x(1025:284002,:);
audiowrite(fNameOut,x,48000);

