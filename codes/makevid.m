writerObj = VideoWriter('test_coupled4.mp4','MPEG-4');
writerObj.FrameRate=20;
open(writerObj);
for K = 1:24900
  if mod(K,20)==0
    K
    filename = sprintf('images/time%i.png', K);
    thisimage = imread(filename);
    writeVideo(writerObj, thisimage);
  end
end
close(writerObj);