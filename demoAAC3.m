function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fnameAACoded)
%function to test Level 3
%Outputs:
%SNR of level 3
%bitrate of Level 3
%compression of level 3
[y,fs]=audioread(fNameIn);
tic;
seqIn=AACoder3(fNameIn,fnameAACoded);
toc;
tic;
x=iAACoder3(seqIn,fNameOut);
toc;
SNR=zeros(2,1);

bitsSum=0;
for k=1:size(seqIn,1)
    if (strcmp(seqIn(k).frameType,'ESH'))
        for i=1:8
            bitsSum=bitsSum+size(seqIn(k).chl.stream{i},2);
            bitsSum=bitsSum+size(seqIn(k).chl.sfc{i},2);
            bitsSum=bitsSum+size(seqIn(k).chr.stream{i},2);
            bitsSum=bitsSum+size(seqIn(k).chr.sfc{i},2);
        end
        bitsSum=bitsSum+2*4*8*64;
        bitsSum=bitsSum+2*8*64; %G 
        bitsSum=bitsSum+2*8*16;       
    else
        bitsSum=bitsSum+size(seqIn(k).chl.stream,2);
        bitsSum=bitsSum+size(seqIn(k).chl.sfc,2);
        bitsSum=bitsSum+size(seqIn(k).chr.stream,2);
        bitsSum=bitsSum+size(seqIn(k).chr.sfc,2);
        bitsSum=bitsSum+2*4*64; %TNSCoeffs
        bitsSum=bitsSum+2*64; %G 
        bitsSum=bitsSum+2*16; %Codebook
    end
    bitsSum=bitsSum+2*3*8; %Window + FrameType
end
KBytesCompressed=bitsSum/(1024*8);
KBytesUncompressed=1117;
Ratio=KBytesUncompressed/KBytesCompressed;
compression=100/Ratio;
z=y-x;
plot(z);
SNR(1)=snr(y(:,1),z(:,1));
SNR(2)=snr(y(:,2),z(:,2));
bitrate=0;

end

