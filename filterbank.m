%Function to covert a frame from time space to MDCT coeffs
function frameF = filterbank(frameT, frameType, winType)
%%Create the window depending of the window type
if (strcmp(frameType,'OLS')) %OLS case
    a=6;
    window=zeros(2048,2);
    if (strcmp(winType,'KDB')) %%KDB window
      num = kaiser(1025, a);
      denum = sum(num(1:1025));
      window = zeros(2048, 2);
      for i = 1:1024
         window(i, :) = sqrt(sum(num(1:i))/denum);                 
      end
      window(1025:2048, :)= window(1024:-1:1,:);
    
    elseif (strcmp(winType,'SIN')) %SIN window
      for n=1:2048
          window(n,:)=sin((pi/2048)*(n-1+1/2));
      end
   end

elseif (strcmp(frameType,'LSS')) %LSS case
    window=zeros(2048,2);
    if (strcmp(winType,'KDB'))
        a=6;
        num = kaiser(1025, a);
        denum = sum(num(1:1025));
        window = zeros(2048, 2);
        for i = 1:1024
           window(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        for n=1025:1472
            window(n,:)=1;
        end
        a=4;
        num = kaiser(129, a);
        denum = sum(num(1:129));
        winTest=zeros(256,2);
        for i = 1:128
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(1473:1600,:)=winTest(128:-1:1,:);
        for n=1601:2048
            window(n,:)=0;
        end
    elseif (strcmp(winType,'SIN'))
        for n=1:1024
            window(n,:)=sin((pi/2048)*(n-1+1/2));
        end
        for n=1025:1472
            window(n,:)=1;
        end
        for n=1473:1600
            window(n,:)=sin((pi/256)*(n+128-1473+1/2));
        end
        for n=1601:2048
            window(n,:)=0;
        end
    end
elseif (strcmp(frameType,'LPS')) %LPS case
    window=zeros(2048,2);
    if (strcmp(winType,'KDB'))
        for n=1:448
            window(n,:)=0;
        end
        a=4;
        num = kaiser(129, a);
        denum = sum(num(1:129));
        winTest=zeros(256,2);
        for i = 1:128
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(449:576,:)=winTest(1:128,:);
        for n=577:1024
            window(n,:)=1;
        end
        a=6;
        num = kaiser(1025, a);
        denum = sum(num(1:1025));
        winTest = zeros(2048, 2);
        for i = 1:1024
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(1025:2048, :)= winTest(1024:-1:1,:);
   
    elseif (strcmp(winType,'SIN')) 
        for n=1:448
            window(n,:)=0;
        end
        for n=449:576
            window(n,:)=sin((pi/256)*(n-449+1/2));
        end
        for n=577:1024
            window(n,:)=1;
        end
        for n=1025:2048
            window(n,:)=sin((pi/2048)*(n-1+1/2));
        end
     end
elseif (strcmp(frameType,'ESH'))  %ESH case
    window=zeros(256,2);
    if (strcmp(winType,'KDB'))
      a=4;
      num = kaiser(129, a);
      denum = sum(num(1:129));
      window = zeros(256, 2);
      for i = 1:128
         window(i, :) = sqrt(sum(num(1:i))/denum);                 
      end
      window(129:256, :)= window(128:-1:1,:);
     
    elseif (strcmp(winType,'SIN'))
      for n=1:256
          window(n,:)=sin((pi/256)*(n-1+1/2));
      end
   end
end

% Now apply the window in the frameT and then use mdct4() to calculate MDCT
% coeffs
if (strcmp(frameType,'ESH')) 
    index=1+448;
    for j=1:8
        windowedFrameT=zeros(256,2);
        windowedFrameT=frameT(index:index+255,:).*window;
        index=index+128;
        frameF(:,1,j)=mdct4(windowedFrameT(:,1));
        frameF(:,2,j)=mdct4(windowedFrameT(:,2));
    end
    
else
    windowedFrameT=frameT.*window;
    frameF(:,1)=mdct4(windowedFrameT(:,1));
    frameF(:,2)=mdct4(windowedFrameT(:,2));
end
end

